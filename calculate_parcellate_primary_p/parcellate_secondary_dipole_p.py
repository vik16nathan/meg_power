#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parcellates secondary dipole P (secondary_p/eeg_dipole_power_computation.py's
fem_p_rms.csv) into Schaefer-600 7-network regions, for the 4 exemplar
subjects that have it, and averages primary + secondary dipole P over just
those same 4 subjects -- feeding schaeffer_rest_power_map.R's Figure 3D
section (primary vs. secondary dipole P, parcellated, same 4 subjects Fig.
3B already uses).

Mirrors parcellate_central_surface_s600.m's region-assignment logic (region
name -> vertex indices comes from each subject's own tess_cortex_central_low.mat
Schaefer_600_17net atlas, mapped 17-network -> 7-network via
s600_17_to_7_mapping.csv, region value = SUM of the per-vertex quantity over
that region's vertices) but reads the vertex-to-region mapping directly in
Python instead of needing a MATLAB/Brainstorm session, since Atlas/Scouts is
a plain struct array (unlike s600_17_to_7.mat's containers.Map, which is why
that one still needed a one-off MATLAB export to CSV -- see
outputs/s600_17_to_7_mapping.csv).
"""
import os

import numpy as np
import pandas as pd
import scipy.io as sio

EXEMPLAR_SUBJECTS = ['sub-0019', 'sub-0023', 'sub-0031', 'sub-0025']
SNR = 3

OUTPUTS_DIR = '/export02/data/vikramn/hbm_manuscript_code/outputs'
OMEGA_ANAT_DIR = '/export02/data/vikramn/OmegaSubset/anat'
PRIMARY_PARCELLATED_DIR = os.path.join(OUTPUTS_DIR, 'primary_p', 'cent_surf_fem_fwd')
FEM_P_RMS_DIR = os.path.join(OUTPUTS_DIR, 'secondary_p', 'fem_p_rms')
OUTPUT_DIR = os.path.join(OUTPUTS_DIR, 'secondary_p', 'cent_surf_fem_fwd')

MEDIAL_WALL_NAMES = {
    'Background+FreeSurfer_Defined_Medial_Wall L',
    'Background+FreeSurfer_Defined_Medial_Wall R',
}


def load_s17_to_s7_mapping():
    df = pd.read_csv(os.path.join(OUTPUTS_DIR, 's600_17_to_7_mapping.csv'))
    return dict(zip(df['s17_name'], df['s7_name']))


def region_vertices_by_s7_name(subject, s17_to_s7):
    """{7-network region name: 0-based vertex indices} for one subject,
    from that subject's own Schaefer_600_17net atlas -- same source
    parcellate_central_surface_s600.m uses, read directly instead of via
    Brainstorm/MATLAB."""
    anat_file = os.path.join(OMEGA_ANAT_DIR, subject, 'tess_cortex_central_low.mat')
    surf = sio.loadmat(anat_file, simplify_cells=True)
    atlas = next(
        a for a in surf['Atlas']
        if 'schaefer' in a['Name'].lower() and '600' in a['Name'] and '17' in a['Name'].lower()
    )
    region_vertices = {}
    for scout in atlas['Scouts']:
        s17_name = scout['Label']
        if s17_name in MEDIAL_WALL_NAMES or s17_name not in s17_to_s7:
            continue  # same skip-list as parcellate_central_surface_s600.m
        s7_name = s17_to_s7[s17_name]
        # MATLAB 1-based vertex indices -> 0-based, matching fem_p_rms.csv's
        # row order (both come from this same subject's low-res cortex surface)
        region_vertices[s7_name] = np.atleast_1d(scout['Vertices']).astype(np.int64) - 1
    return region_vertices


def parcellate_subject(subject, per_vertex_values, region_vertices):
    """Region value = SUM of the per-vertex quantity over that region's
    vertices, matching parcellate_central_surface_s600.m (physically
    additive: power in a region is the sum of its vertices' power, not a
    mean)."""
    return {name: per_vertex_values[idx].sum() for name, idx in region_vertices.items()}


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    s17_to_s7 = load_s17_to_s7_mapping()

    per_subject_secondary = {}
    per_subject_primary = {}
    for subject in EXEMPLAR_SUBJECTS:
        fem_p_rms_file = os.path.join(FEM_P_RMS_DIR, f'{subject}_snr{SNR}_fem_p_rms.csv')
        primary_file = os.path.join(PRIMARY_PARCELLATED_DIR, f's600_{subject}_snr{SNR}_p_pri_2min_rest.csv')
        if not (os.path.exists(fem_p_rms_file) and os.path.exists(primary_file)):
            print(f'Skipping {subject}: missing {fem_p_rms_file} and/or {primary_file}.')
            continue

        secondary_dipole_p = pd.read_csv(fem_p_rms_file, header=None).to_numpy().flatten()
        region_vertices = region_vertices_by_s7_name(subject, s17_to_s7)
        region_secondary = parcellate_subject(subject, secondary_dipole_p, region_vertices)

        out_file = os.path.join(OUTPUT_DIR, f's600_{subject}_snr{SNR}_p_sec_dipole_2min_rest.csv')
        pd.DataFrame({'region': list(region_secondary.keys()), 'value': list(region_secondary.values())}) \
            .to_csv(out_file, index=False)
        print(f'{subject}: wrote {out_file} ({len(region_secondary)} regions)')

        per_subject_secondary[subject] = region_secondary
        per_subject_primary[subject] = dict(zip(
            pd.read_csv(primary_file)['region'], pd.read_csv(primary_file)['value']))

    if not per_subject_secondary:
        raise FileNotFoundError('No exemplar subject has both fem_p_rms.csv and parcellated primary P.')

    # Average over exactly these subjects (not the full 30-subject cohort's
    # existing subavg file), so both quantities in Figure 3D are averaged
    # over the identical n=4 -- a fair, matched comparison.
    subjects_done = sorted(per_subject_secondary)
    all_regions = sorted(set.union(*(set(d) for d in per_subject_secondary.values())))

    def average_over(per_subject_dict):
        rows = []
        for region in all_regions:
            values = [per_subject_dict[s][region] for s in subjects_done if region in per_subject_dict[s]]
            if values:
                rows.append({'region': region, 'value': np.mean(values)})
        return pd.DataFrame(rows)

    sec_avg = average_over(per_subject_secondary)
    pri_avg = average_over(per_subject_primary)

    sec_out = os.path.join(OUTPUT_DIR, f's600_p_sec_dipole_2min_rest_subavg_snr{SNR}_exemplars.csv')
    pri_out = os.path.join(OUTPUT_DIR, f's600_p_pri_2min_rest_subavg_snr{SNR}_exemplars.csv')
    sec_avg.to_csv(sec_out, index=False)
    pri_avg.to_csv(pri_out, index=False)
    print(f'\nAveraged over {subjects_done}:')
    print(f'  wrote {sec_out} ({len(sec_avg)} regions)')
    print(f'  wrote {pri_out} ({len(pri_avg)} regions)')


if __name__ == '__main__':
    main()
