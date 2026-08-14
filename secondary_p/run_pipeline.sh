#!/bin/bash
# Runs the secondary_p pipeline end to end, in order:
#   1. convert_fem_mat_to_vc.py   (any python3 -- fast, no duneuropy needed)
#   2. eeg_dipole_power_computation.py   (needs python3.8 + duneuropy; slow --
#      ~1 hour/subject with N_WORKERS=20 parallel dipole solves)
#   3. eeg_dipole_vsum_sec_curr.py       (same python3.8 requirement; also
#      ~1 hour/subject, run only for the exemplar/available subjects used
#      by the tissue-breakdown pie charts)
#   4. generate_figures.py       (any python3 -- fast; produces everything
#      under figures/, degrading gracefully for any subject/figure whose
#      upstream data isn't ready yet)
#
# Each stage is idempotent (skips subjects/files already done), so this can
# be safely re-run to pick up wherever a prior run left off.
#
# Usage: nohup ./run_pipeline.sh > log/run_pipeline_$(date +%Y%m%d_%H%M%S).log 2>&1 &

set -eu
cd "$(dirname "$0")"

PY38=/usr/bin/python3.8
LOGDIR="log"
mkdir -p "$LOGDIR"

run_stage() {
    local label="$1"; shift
    local logfile="$LOGDIR/${label}_$(date +%Y%m%d_%H%M%S).log"
    echo "[$(date)] === $label: starting (log: $logfile) ==="
    if "$@" > "$logfile" 2>&1; then
        echo "[$(date)] === $label: completed ==="
    else
        echo "[$(date)] === $label: FAILED (see $logfile) ==="
        exit 1
    fi
}

run_stage convert                      python3 convert_fem_mat_to_vc.py
run_stage eeg_dipole_power_computation "$PY38" eeg_dipole_power_computation.py
run_stage eeg_dipole_vsum_sec_curr     "$PY38" eeg_dipole_vsum_sec_curr.py
run_stage generate_figures             python3 generate_figures.py

echo "[$(date)] Pipeline finished."
