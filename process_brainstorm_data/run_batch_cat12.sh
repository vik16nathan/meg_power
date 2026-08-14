#!/bin/bash
# Runs batch_cat12.m headlessly and auto-restarts MATLAB if it gets killed or
# hangs mid-run. Known causes so far: an OOM kill during CAT12 segmentation,
# a CAT12 telemetry call with a broken timeout/retry (fixed directly, see
# cat_defaults.m send_info=0), and at least one more hang whose exact call
# site wasn't pinned down (candidates: a Brainstorm network check, or NFS
# I/O stalls -- /home/bic is a hard-mounted NFS4 share, and spm_jobman
# ('initcfg') does a huge recursive file scan across it).
#
# 
# To check for timeouts, this script watches the log file: every real processing phase in this pipeline prints
# something, so no new log output for a long stretch means it's actually
# stuck, not just busy. 
#
# batch_cat12.m skips subjects that already have completed CAT12 results, so
# each restart just resumes with the next unfinished subject instead of
# losing the whole overnight batch to one bad subject.
#
# Usage: nohup ./run_batch_cat12.sh > log/run_batch_cat12.out 2>&1 &

set -u
# Job control (monitor mode) so each backgrounded MATLAB run gets its own
# process group -- required for the `kill -- -$mpid` process-group kill below
# to work instead of accidentally targeting this script's own group.
set -m
cd "$(dirname "$0")"

MAX_ATTEMPTS=50
STALL_SECS=$((10 * 60))   # Kill only if the log hasn't grown at all for this long.
POLL_SECS=15
LOGDIR="log"
mkdir -p "$LOGDIR"

for ((attempt=1; attempt<=MAX_ATTEMPTS; attempt++)); do
    LOGFILE="$LOGDIR/batch_cat12_$(date +%Y%m%d_%H%M%S).log"
    echo "[$(date)] Attempt $attempt/$MAX_ATTEMPTS -- launching MATLAB (log: $LOGFILE)"

    # stdin must be detached (< /dev/null): with `set -m` job control active,
    # a background job that still has the controlling terminal as stdin gets
    # SIGTTIN-stopped the moment it touches stdin, freezing it silently
    # before it can print anything -- this caused the empty-log hangs.
    matlab21b -nodisplay -nosplash -r "run('batch_cat12.m'); exit;" < /dev/null > "$LOGFILE" 2>&1 &
    mpid=$!

    # Mirror the log to the terminal live (for foreground/interactive runs).
    tail -n +1 -f "$LOGFILE" --pid="$mpid" &
    tailpid=$!

    status=-1
    while true; do
        if ! kill -0 "$mpid" 2>/dev/null; then
            wait "$mpid"
            status=$?
            break
        fi
        sleep "$POLL_SECS"
        now=$(date +%s)
        last_mod=$(stat -c %Y "$LOGFILE" 2>/dev/null) || last_mod=$now
        idle=$((now - last_mod))
        if [ "$idle" -ge "$STALL_SECS" ]; then
            echo "[$(date)] No new log output for ${STALL_SECS}s (stuck, not just slow) -- killing and moving to the next subject."
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
        echo "[$(date)] batch_cat12.m exited normally -- all subjects processed or skipped."
        exit 0
    fi

    echo "[$(date)] MATLAB exited with status $status -- restarting to resume remaining subjects."
    sleep 10
done

echo "[$(date)] Reached MAX_ATTEMPTS ($MAX_ATTEMPTS) without a clean exit. Check $LOGDIR for the last log and investigate."
exit 1
