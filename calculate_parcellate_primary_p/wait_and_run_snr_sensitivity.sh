#!/bin/bash
# Polls until rms_idip_rest.m has produced all 31 subjects x 3 SnrFixed
# RMS I_dip files, then runs the SNR sensitivity analysis on them.
set -u

RMS_DIR=/export02/data/vikramn/hbm_manuscript_code/outputs/primary_p/idip_rms
SUBJECTS=(sub-0002 sub-0008 sub-0009 sub-0010 sub-0011 sub-0012 sub-0014 sub-0015 \
          sub-0016 sub-0018 sub-0019 sub-0020 sub-0021 sub-0022 sub-0023 sub-0024 sub-0025 \
          sub-0026 sub-0028 sub-0029 sub-0030 sub-0031 sub-0032 sub-0033 sub-0034 sub-0035 sub-0036 \
          sub-0037 sub-0039 sub-0040 sub-0041)
SNRS=(1 3 5)

while true; do
  missing=0
  for s in "${SUBJECTS[@]}"; do
    for snr in "${SNRS[@]}"; do
      f="$RMS_DIR/${s}_snr${snr}_idip_rms.csv"
      [ -f "$f" ] || missing=$((missing+1))
    done
  done
  echo "$(date): missing $missing/93 RMS I_dip files"
  if [ "$missing" -eq 0 ]; then
    break
  fi
  sleep 60
done

echo "$(date): all 93 RMS I_dip files present -- running SNR sensitivity analysis"
cd /export02/data/vikramn/hbm_manuscript_code/calculate_parcellate_primary_p
Rscript snr_sensitivity_idip_rms.R
echo "EXIT: $?"
