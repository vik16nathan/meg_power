#!/usr/bin/env python3
"""Builds rms_2min_cent_pri_power_met_broad.csv: primary dipole P/I
(parcellate_central_surface_s600.m's Schaeffer-600 subject-average output)
joined with neuromaps CMRO2/CMRGlu region metabolism, for
schaeffer_rest_power_map.R.

Replaces the equivalent join/rename/to_csv logic that used to live at the
end of metabolism_correlations/neuromaps_rest_final.ipynb (that notebook now
just reads this script's output CSV back in for its own downstream
analysis).

Prerequisites:
- parcellate_central_surface_s600.m must have written
  s600_p_pri_2min_rest_subavg_snr<SNR>.csv and
  s600_i_pri_2min_rest_subavg_snr<SNR>.csv (columns: region, value) to
  OUTPUTS_DIR/primary_p/cent_surf_fem_fwd/.
- A neuromaps CMRO2/CMRGlu reference CSV (columns: region, cmro2, cmrglu),
  matching the same Schaefer-600 7-network parcellation. This is a one-time,
  subject-independent annotation fetch (see
  metabolism_correlations/neuromaps_rest_final.ipynb's neuromaps
  fetch_annotation cells, still in the notebook) -- NOT regenerated here.
  Region-name conventions between the neuromaps fetch and this pipeline's
  Brainstorm/CAT12-derived region labels (e.g. "Cont_Cing_1 L") are not
  guaranteed to match; harmonize_region_name() below is a no-op hook to
  patch that if/when regions fail to match (see the printed warning).
"""
import argparse
from pathlib import Path

import pandas as pd

OUTPUTS_DIR = Path('/export02/data/vikramn/hbm_manuscript_code/outputs')
PRIMARY_P_DIR = OUTPUTS_DIR / 'primary_p' / 'cent_surf_fem_fwd'


def harmonize_region_name(name):
    """No-op hook: patch region-name mismatches between the MEG pipeline's
    Schaefer-600 labels and the metabolism CSV's labels here if needed."""
    return name


def build(snr, metabolism_csv, output_csv):
    p_file = PRIMARY_P_DIR / f's600_p_pri_2min_rest_subavg_snr{snr}.csv'
    i_file = PRIMARY_P_DIR / f's600_i_pri_2min_rest_subavg_snr{snr}.csv'
    for f in (p_file, i_file, metabolism_csv):
        if not f.exists():
            raise FileNotFoundError(
                f'{f} not found -- run parcellate_central_surface_s600.m (for the P/I '
                f'files) or see this script\'s module docstring (for the metabolism '
                f'reference file) first.')

    meg_p = pd.read_csv(p_file).rename(columns={'value': 'mean_p'})
    meg_i = pd.read_csv(i_file).rename(columns={'value': 'mean_i'})
    metabolism = pd.read_csv(metabolism_csv)[['region', 'cmro2', 'cmrglu']]

    for df in (meg_p, meg_i, metabolism):
        df['region'] = df['region'].map(harmonize_region_name)

    # A duplicated 'region' key in any input silently turns merge() into a
    # cartesian join for that key (e.g. 1 metabolism duplicate + otherwise
    # clean 600-region P/I files still produces 601+ output rows, with the
    # duplicated region's P/I value paired against multiple different
    # metabolism values) -- fail loudly instead of writing a corrupted CSV.
    for name, df in (('meg_p', meg_p), ('meg_i', meg_i), ('metabolism', metabolism)):
        dupes = df.loc[df['region'].duplicated(keep=False), 'region'].unique().tolist()
        if dupes:
            raise ValueError(f'{name} has duplicate region name(s), fix the source file before merging: {dupes}')

    merged = meg_p.merge(metabolism, on='region', how='inner').merge(meg_i, on='region', how='inner')

    missing = sorted(set(meg_p['region']) - set(merged['region']))
    if missing:
        print(f'WARNING: {len(missing)}/{len(meg_p)} region(s) in {p_file.name} had no '
              f'metabolism match (check harmonize_region_name() above): {missing[:5]}'
              f'{"..." if len(missing) > 5 else ""}')

    merged = merged[['region', 'mean_p', 'cmro2', 'cmrglu', 'mean_i']]
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_csv, index=False)
    print(f'Saved {output_csv} ({len(merged)} regions).')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--snr', type=int, default=3,
                         help='SnrFixed value of the primary P/I files to use (default: 3).')
    parser.add_argument('--metabolism-csv', type=Path,
                         default=OUTPUTS_DIR / 'raichle_metabolism_s600_7.csv',
                         help='CMRO2/CMRGlu neuromaps reference CSV (region, cmro2, cmrglu).')
    parser.add_argument('--output-csv', type=Path,
                         default=PRIMARY_P_DIR / 'rms_2min_cent_pri_power_met_broad.csv',
                         help='Output path for the merged CSV.')
    args = parser.parse_args()
    build(args.snr, args.metabolism_csv, args.output_csv)
