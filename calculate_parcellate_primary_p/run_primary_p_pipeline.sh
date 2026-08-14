#!/bin/bash
# Runs the calculate_parcellate_primary_p pipeline headlessly, in order:
#   rms_idip_rest.m
#   calculate_fem_column_power.m
#   parcellate_central_surface_s600.m
#   build_pri_power_met_csv.py
#   schaeffer_rest_power_map.R
#   snr_sensitivity_supp_fig.R (supplementary, non-fatal)
#
# Mirrors process_brainstorm_data/run_pipeline.sh's conventions: same
# 32-subject list, the same restart-on-stall pattern for MATLAB stages (each
# .m script is idempotent -- it skips subjects/SnrFixed combinations already
# done -- so a restart just resumes instead of redoing work), and the same
# rule of never running a Brainstorm/MATLAB session here at the same time as
# process_brainstorm_data/run_pipeline.sh: both would open the OmegaSubset
# protocol, and concurrent writers risk database corruption (this happened
# for real once already -- see process_brainstorm_data/bst_headless_init.m).
#
# Usage: nohup ./run_primary_p_pipeline.sh > log/run_primary_p_pipeline.out 2>&1 &

set -u
set -m
cd "$(dirname "$0")"

OUTPUTS_DIR="/export02/data/vikramn/hbm_manuscript_code/outputs"

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
        # errors, so a failed stage doesn't look like a false-positive success.
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

# run_r_stage <script.R>
# R stages here don't touch Brainstorm and don't hang the way a headless
# Brainstorm session can, so just run once and report status -- no
# stall-restart loop needed.
run_r_stage() {
    local script="$1"
    local stage="${script%.R}"
    local logfile="$LOGDIR/${stage}_$(date +%Y%m%d_%H%M%S).log"
    echo "[$(date)] $stage: launching Rscript (log: $logfile)"
    Rscript "$script" > "$logfile" 2>&1
    local status=$?
    if [ "$status" -eq 0 ]; then
        echo "[$(date)] $stage: completed normally."
        return 0
    fi
    echo "[$(date)] $stage: Rscript exited with status $status."
    return 1
}

echo "[$(date)] Waiting for process_brainstorm_data/run_pipeline.sh to finish (never run a Brainstorm session concurrently with it)..."
while pgrep -f "process_brainstorm_data/run_pipeline\.sh" > /dev/null; do
    sleep "$POLL_SECS"
done
echo "[$(date)] process_brainstorm_data/run_pipeline.sh is not running -- proceeding."

# rms_idip_rest.m needs process_meg_compute_mn_source_fem.m's kernel manifest
# and preprocess_meg_notch_highpass.m's condition manifest.
if [ ! -f "$OUTPUTS_DIR/mn_kernel_manifest.mat" ] || [ ! -f "$OUTPUTS_DIR/preproc_condition_manifest.mat" ]; then
    echo "[$(date)] Missing mn_kernel_manifest.mat and/or preproc_condition_manifest.mat in $OUTPUTS_DIR."
    echo "[$(date)] Run process_brainstorm_data/run_pipeline.sh through process_meg_compute_mn_source_fem.m first. Stopping."
    exit 1
fi
run_stage rms_idip_rest.m || exit 1

# calculate_fem_column_power.m needs thickness_areas_omega_fem.mat and each
# subject's uniform-grey-matter-resistivity column resistance file.
if [ ! -f "$OUTPUTS_DIR/thickness_areas_omega_fem.mat" ]; then
    echo "[$(date)] Missing thickness_areas_omega_fem.mat in $OUTPUTS_DIR."
    echo "[$(date)] Run process_brainstorm_data/calculate_dipole_p_inputs.m first. Stopping."
    exit 1
fi
if ! ls "$OUTPUTS_DIR"/resistances/*_vertex_res_gm.mat >/dev/null 2>&1; then
    echo "[$(date)] No *_vertex_res_gm.mat files in $OUTPUTS_DIR/resistances."
    echo "[$(date)] Run process_brainstorm_data/parcellate_cortical_column_res_vol.m first. Stopping."
    exit 1
fi
run_stage calculate_fem_column_power.m || exit 1

# parcellate_central_surface_s600.m needs the Schaefer-600 17-to-7 mapping.
if [ ! -f "$OUTPUTS_DIR/s600_17_to_7.mat" ] || [ ! -f "$OUTPUTS_DIR/s600_region_names.mat" ]; then
    echo "[$(date)] Missing s600_17_to_7.mat and/or s600_region_names.mat in $OUTPUTS_DIR."
    echo "[$(date)] Run process_brainstorm_data/generate_schaefer_17_to_7_mapping.m first. Stopping."
    exit 1
fi
run_stage parcellate_central_surface_s600.m || exit 1

# build_pri_power_met_csv.m needs a neuromaps CMRO2/CMRGlu reference CSV
# (see its module docstring) -- not produced by any stage in this pipeline.
if [ ! -f "$OUTPUTS_DIR/raichle_metabolism_s600_7.csv" ]; then
    echo "[$(date)] Missing raichle_metabolism_s600_7.csv in $OUTPUTS_DIR."
    echo "[$(date)] See build_pri_power_met_csv.py's docstring (neuromaps_rest_final.ipynb's fetch_annotation cells produce this). Stopping."
    exit 1
fi
echo "[$(date)] build_pri_power_met_csv: launching Python (log: $LOGDIR/build_pri_power_met_csv_$(date +%Y%m%d_%H%M%S).log)"
python3 build_pri_power_met_csv.py > "$LOGDIR/build_pri_power_met_csv_$(date +%Y%m%d_%H%M%S).log" 2>&1 || {
    echo "[$(date)] build_pri_power_met_csv.py failed. Stopping."
    exit 1
}

run_r_stage schaeffer_rest_power_map.R || exit 1

# snr_sensitivity_supp_fig.R is a supplementary analysis (regularization
# parameter sensitivity), not part of the core Fig. 1 pipeline -- a failure
# here shouldn't block the main outputs above, so don't exit on failure.
run_r_stage snr_sensitivity_supp_fig.R || echo "[$(date)] snr_sensitivity_supp_fig.R failed (non-fatal, supplementary figure only)."

echo "[$(date)] Pipeline finished."
