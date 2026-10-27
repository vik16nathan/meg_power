#!/usr/bin/env python3
"""Histogram of whole-brain total primary dipole P (one point per subject),
at a single SnrFixed value.

Reads parcellate_central_surface_s600.m's per-subject, per-region CSVs
(outputs/primary_p/cent_surf_fem_fwd/s600_<subject>_snr<SNR>_p_pri_2min_rest.csv,
columns: region, value) and sums each subject's regions to get that
subject's whole-brain total P, matching the summary-statistics method used
throughout this pipeline's reporting.

Discovers whichever subjects currently have output for SNR (via glob)
rather than assuming a fixed subject count, so this can be run at any point
-- including while rms_idip_rest.m is still filling in the remaining
subjects -- and will simply show however many are done so far.
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # headless: no display, just save to file
import matplotlib.pyplot as plt
import pandas as pd

SNR = 3  # SnrFixed value to plot (change here to re-target 1 or 5)

# Confirmed outlier (~11x the excluding-outlier mean whole-brain current
# dipole moment, R and cortical thickness both normal -- flagged but not
# root-caused, likely subject-specific data quality rather than a pipeline
# bug). Reported separately below rather than silently dropped.
OUTLIER_SUBJECT = 'sub-0011'

OUTPUTS_DIR = Path('/export02/data/vikramn/hbm_manuscript_code/outputs')
POWER_DIR = OUTPUTS_DIR / 'primary_p' / 'cent_surf_fem_fwd'  # parcellate_central_surface_s600.m output
FIG_DIR = OUTPUTS_DIR / 'primary_p' / 'figures'


def load_subject_totals(snr):
    files = sorted(POWER_DIR.glob(f's600_sub-*_snr{snr}_p_pri_2min_rest.csv'))
    if not files:
        raise FileNotFoundError(
            f'No s600_sub-*_snr{snr}_p_pri_2min_rest.csv files in {POWER_DIR} -- '
            f'run parcellate_central_surface_s600.m first.')

    rows = []
    for f in files:
        subject = f.stem.split('_')[1]  # s600_<subject>_snr<N>_... -> <subject>
        df = pd.read_csv(f)
        rows.append({'subject': subject, 'total_p_W': df['value'].sum()})
    return pd.DataFrame(rows).sort_values('subject').reset_index(drop=True)


def print_summary(label, values):
    print(f'{label} (n={len(values)}):')
    print(f'  mean   = {values.mean():.2e} W')
    print(f'  median = {values.median():.2e} W')
    print(f'  min    = {values.min():.2e} W')
    print(f'  max    = {values.max():.2e} W')


def main():
    df = load_subject_totals(SNR)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(df['total_p_W'], bins=15, color='#4C72B0', edgecolor='white')
    ax.set_xlabel('Whole-brain total primary dipole P (W)')
    ax.set_ylabel('Number of subjects')
    ax.set_title(f'Total primary dipole P across subjects (SnrFixed={SNR}, n={len(df)})')

    out_file = FIG_DIR / f'total_p_histogram_snr{SNR}.png'
    fig.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f'n = {len(df)} subjects')
    print(df.to_string(index=False))
    print(f'\nSaved {out_file}')

    # --- Main-text summary statistics -----------------------------------
    print('\n=== Main-text summary statistics (whole-brain total primary dipole power) ===')
    print('sub-0013 excluded from the cohort throughout (DUNEuro FEM forward model')
    print('failed to converge even after widening the BEM brain-inner margin to')
    print('3/6/9 mm; see SUPPLEMENTARY_METHODS.md Section 1) -- not a statistical')
    print('exclusion, so it is not present in the totals below at all.\n')
    print_summary(f'All subjects', df['total_p_W'])
    if OUTLIER_SUBJECT in df['subject'].values:
        df_excl = df[df['subject'] != OUTLIER_SUBJECT]
        print(f'\n{OUTLIER_SUBJECT} flagged as an outlier (elevated whole-brain current')
        print('dipole moment with normal R/thickness; excluded from the statistics below,')
        print('shown separately here):')
        print(f'  {OUTLIER_SUBJECT} total_p_W = {df.loc[df.subject == OUTLIER_SUBJECT, "total_p_W"].iloc[0]:.2e} W')
        print()
        print_summary(f'Excluding {OUTLIER_SUBJECT}', df_excl['total_p_W'])


if __name__ == '__main__':
    main()
