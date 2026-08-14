#!/bin/bash
# Waits for the sub-0002 validation run (run_pipeline.sh, PID passed as $1)
# to exit, sanity-checks its vsum/pie-chart output, and if that looks sane,
# runs the remaining subjects (everything with full inputs except sub-0002)
# through convert -> eeg_dipole_power_computation -> generate_figures only
# -- deliberately skipping eeg_dipole_vsum_sec_curr.py, which is scoped to
# just the 4 exemplar subjects and run separately/later, per instruction.
#
# Usage: nohup ./run_remaining_subjects.sh <pid-to-wait-for> > log/run_remaining_subjects.out 2>&1 &

set -u
cd "$(dirname "$0")"

WAIT_PID="$1"
LOGDIR="log"
mkdir -p "$LOGDIR"
PY38=/usr/bin/python3.8
REMAINING_SUBJECTS="sub-0008,sub-0019,sub-0022,sub-0029,sub-0032"

echo "[$(date)] Waiting for PID $WAIT_PID (sub-0002 validation run) to exit..."
while kill -0 "$WAIT_PID" 2>/dev/null; do
    sleep 30
done
echo "[$(date)] PID $WAIT_PID has exited."

if ! grep -q "Pipeline finished." log/run_pipeline_test_sub0002.out; then
    echo "[$(date)] sub-0002 run did not report 'Pipeline finished.' -- not proceeding. Check log/run_pipeline_test_sub0002.out."
    exit 1
fi

echo "[$(date)] Sanity-checking sub-0002 vsum/pie-chart output..."
if ! "$PY38" - <<'PYEOF'
import sys
import numpy as np

p = np.load('../outputs/secondary_p/vsum_sec_curr/sub-0002_p_output.npz')
v_els, i_sq_els, p_els = p['v_els'], p['i_sq_els'], p['p_els']

problems = []
for name, arr in [('v_els', v_els), ('i_sq_els', i_sq_els), ('p_els', p_els)]:
    if np.isnan(arr).any():
        problems.append(f'{name} has NaNs')
    if (arr < 0).any():
        problems.append(f'{name} has negative values')
if not (0.001 < v_els.sum() < 0.02):  # plausible whole-head volume in m^3
    problems.append(f'v_els.sum()={v_els.sum():.6f} m^3 outside plausible whole-head range')

if problems:
    print('SANITY CHECK FAILED:', '; '.join(problems))
    sys.exit(1)
print(f'Sanity check passed: v_els.sum()={v_els.sum():.6f} m^3, '
      f'p_els.sum()={p_els.sum():.3e}, no NaN/negative values.')
PYEOF
then
    echo "[$(date)] Sanity check failed -- NOT launching remaining subjects. Investigate sub-0002's vsum output first."
    exit 1
fi

echo "[$(date)] Sanity check passed. Launching remaining subjects: $REMAINING_SUBJECTS"

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

run_stage convert python3 convert_fem_mat_to_vc.py
SECONDARY_P_SUBJECTS="$REMAINING_SUBJECTS" run_stage eeg_dipole_power_computation "$PY38" eeg_dipole_power_computation.py
run_stage generate_figures python3 generate_figures.py

echo "[$(date)] Remaining-subjects run finished."
