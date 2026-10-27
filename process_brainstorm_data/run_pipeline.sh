#!/bin/bash
# Waits for the currently-running batch_cat12.m (driven by run_batch_cat12.sh)
# to finish, then runs the rest of the pipeline headlessly, in order:
#   preprocess_meg_notch_highpass.m
#   build_duneuro_forward_models.m
#   process_meg_compute_mn_source_fem.m
#   calculate_dipole_p_inputs.m
#   parcellate_cortical_column_res_vol.m   (only if its required data files are present)
#   cell_res_avg.m                         (only if its required data files are present)
#
# Each stage uses the same restart-on-stall pattern as run_batch_cat12.sh: if
# a stage's log stops growing for STALL_SECS, MATLAB is killed and restarted
# for that stage (every stage's .m script is idempotent -- it skips subjects
# already done -- so a restart just resumes instead of redoing work).
#
# Runs stages sequentially, never two Brainstorm/MATLAB sessions at once
# against the same protocol, to avoid the database-corruption risk that
# comes with concurrent writers to protocol.mat.
#
# Usage: nohup ./run_pipeline.sh > log/run_pipeline.out 2>&1 &

set -u
set -m
cd "$(dirname "$0")"

LOGDIR="log"
mkdir -p "$LOGDIR"
STALL_SECS=$((10 * 60))
POLL_SECS=15
MAX_ATTEMPTS=50

# run_stage <script.m>
# Runs one .m script headlessly, restarting on a stall/crash, up to
# MAX_ATTEMPTS times. Returns 0 on a clean MATLAB exit, 1 if attempts run out.
run_stage() {
    local script="$1"
    local stage="${script%.m}"

    for ((attempt=1; attempt<=MAX_ATTEMPTS; attempt++)); do
        local logfile="$LOGDIR/${stage}_$(date +%Y%m%d_%H%M%S).log"
        echo "[$(date)] $stage: attempt $attempt/$MAX_ATTEMPTS -- launching MATLAB (log: $logfile)"

        # -batch (not -r): guarantees a non-zero exit code if the script
        # errors. -r would still run a trailing `exit;` after an uncaught
        # error, making every failed stage look like a false-positive
        # success (this happened for real: all 5 downstream stages "passed"
        # in ~2 minutes total on 2026-07-29 while actually erroring out
        # immediately on missing bst_headless_init).
        matlab21b -batch "run('$script')" < /dev/null > "$logfile" 2>&1 &
        local mpid=$!

        tail -n +1 -f "$logfile" --pid="$mpid" &
        local tailpid=$!

        local status=-1
        while true; do
            if ! kill -0 "$mpid" 2>/dev/null; then
                wait "$mpid"
                status=$?
                break
            fi
            sleep "$POLL_SECS"
            local now last_mod idle
            now=$(date +%s)
            last_mod=$(stat -c %Y "$logfile" 2>/dev/null) || last_mod=$now
            idle=$((now - last_mod))
            if [ "$idle" -ge "$STALL_SECS" ]; then
                echo "[$(date)] $stage: no new log output for ${STALL_SECS}s (stuck) -- killing and retrying."
                kill -TERM -- -"$mpid" 2>/dev/null
                sleep 5
                kill -KILL -- -"$mpid" 2>/dev/null
                wait "$mpid" 2>/dev/null
                status=124
                break
            fi
        done
        kill "$tailpid" 2>/dev/null
        wait "$tailpid" 2>/dev/null

        if [ "$status" -eq 0 ]; then
            echo "[$(date)] $stage: completed normally."
            return 0
        fi
        echo "[$(date)] $stage: MATLAB exited with status $status -- restarting."
        sleep 10
    done

    echo "[$(date)] $stage: reached MAX_ATTEMPTS without a clean exit. Stopping pipeline here."
    return 1
}

echo "[$(date)] Waiting for run_batch_cat12.sh to finish..."
# Wait on the wrapper script itself, not the inner MATLAB process -- the
# wrapper kills and restarts MATLAB on stalls, which would otherwise create
# gaps where no MATLAB process exists yet the batch isn't actually done.
while pgrep -f "run_batch_cat12\.sh" > /dev/null; do
    sleep "$POLL_SECS"
done
echo "[$(date)] run_batch_cat12.sh is no longer running -- starting downstream pipeline."

run_stage preprocess_meg_notch_highpass.m || exit 1
run_stage build_duneuro_forward_models.m || exit 1
run_stage process_meg_compute_mn_source_fem.m || exit 1
run_stage calculate_dipole_p_inputs.m || exit 1

# parcellate_cortical_column_res_vol.m needs s600_17_to_7.mat / s600_region_names.mat.
# generate_schaefer_17_to_7_mapping.m builds these from the first subject
# that has both the 17-network and 7-network Schaefer-600 atlases (the
# 7-network one was only just added to cat_defaults.m, so this may still
# come up empty if batch_cat12.m finished before that took effect for any
# subject -- in that case it errors out and this stage is skipped).
PBD="$(pwd)"
OUTPUTS_DIR="/export02/data/vikramn/hbm_manuscript_code/outputs"
if [ ! -f "$OUTPUTS_DIR/s600_17_to_7.mat" ] || [ ! -f "$OUTPUTS_DIR/s600_region_names.mat" ]; then
    run_stage generate_schaefer_17_to_7_mapping.m
fi
if [ -f "$OUTPUTS_DIR/s600_17_to_7.mat" ] && [ -f "$OUTPUTS_DIR/s600_region_names.mat" ]; then
    run_stage parcellate_cortical_column_res_vol.m || exit 1
else
    echo "[$(date)] SKIPPING parcellate_cortical_column_res_vol.m: no subject with both Schaefer-600 atlases yet."
fi

RESDIR="$OUTPUTS_DIR/resistances"
if [ -f "$PBD/cell_props_filtered_redo.csv" ] && ls "$RESDIR"/*_vertex_res_gm.mat >/dev/null 2>&1; then
    run_stage cell_res_avg.m || exit 1
else
    echo "[$(date)] SKIPPING cell_res_avg.m: missing cell_props_filtered_redo.csv and/or no *_vertex_res_gm.mat from parcellate_cortical_column_res_vol.m in $RESDIR"
fi

echo "[$(date)] Pipeline finished."
