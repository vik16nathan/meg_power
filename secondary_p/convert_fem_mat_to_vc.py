#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Converts each subject's Brainstorm/DUNEuro 3-layer FEM head mesh into the
.npz volume-conductor format consumed by eeg_dipole_power_computation.py.

Reads directly from the live Brainstorm protocol database
(OMEGA_ANAT_DIR/<subject>/ for the FEM mesh, OMEGA_DATA_DIR/<subject>/<cond>/
for headmodel_surf_duneuro.mat, wherever <cond> is -- found by search, since
which raw recording's condition folder ended up holding the head model
varies per subject) rather than a hand-curated backup, so this picks up
every subject build_duneuro_forward_models.m has finished, not just the
original 8 whose mesh files were separately copied into
outputs/mergemesh_fem_backup/.

Two DUNEuro FEM meshing methods were used across the cohort (see
process_brainstorm_data/build_duneuro_forward_models.m for why): meshes
built with 'iso2mesh-2021' already carry named TissueLabels
(['brain','skull','scalp']); meshes built with 'iso2mesh-2026' carry raw
numeric tetgen region IDs (['1','2','3']) whose ordering has no consistent
meaning across subjects, and must be reclassified by testing which of that
subject's own scalp/outer-skull/inner-skull BEM surfaces each tissue
region's boundary vertices coincide with -- ported from that same MATLAB
script's ClassifyFemTissueLabels().
"""
import glob
import os

import numpy as np
import scipy.io as sio
from scipy.spatial import cKDTree

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OMEGA_ANAT_DIR = '/export02/data/vikramn/OmegaSubset/anat'
OMEGA_DATA_DIR = '/export02/data/vikramn/OmegaSubset/data'
OUTPUT_DIR = os.path.join(SCRIPT_DIR, '..', 'outputs', 'secondary_p', 'vc_3layer')
# OmegaSubset/anat holds the whole OMEGA dataset, not just this manuscript's
# cohort -- outputs/primary_p/idip_rms/ (rms_idip_rest.m's output) is the
# authoritative list of subjects actually analyzed (n=31, sub-0013 already
# excluded there), so use it to scope which subjects to even attempt.
IDIP_RMS_DIR = os.path.join(SCRIPT_DIR, '..', 'outputs', 'primary_p', 'idip_rms')

# sub-0013: DUNEuro forward-model computation failed outright (see
# SUPPLEMENTARY_METHODS.md) -- excluded from primary P too, final cohort is
# n=31. sub-0009: DUNEuro succeeded for EEG (has idip_rms/primary P) but its
# only head model on disk is headmodel_surf_duneuro_OPENMEEG.mat -- a
# boundary-element (surface-only) fallback with no volumetric FEM mesh/
# per-element conductivity, which this pipeline's power integral needs.
# Secondary P is therefore only computable for 30 of the 31 subjects.
EXCLUDED_SUBJECTS = {'sub-0009', 'sub-0013'}

# Matches what Brainstorm's own DUNEuro forward modeling actually used for
# these subjects' head models (build_duneuro_forward_models.m's
# GetDefaultConductivitySimbio 'simbio' reference values) -- confirmed
# directly from each subject's headmodel_surf_duneuro.mat History field
# ("FemCond: isotropic, ..."), consistently {0.33, 0.008, 0.43} (just
# permuted per subject to match that subject's raw tissue-label order)
# across every subject checked. NOT the manuscript text's stated 0.33/
# 0.0042/0.33 (citing Oostenveld et al. 2011) -- that combination appears in
# neither Brainstorm's actual output nor here; using Brainstorm's real
# values instead keeps this pipeline consistent with the primary-P Gain
# matrix computed from the exact same head models.
CONDUCTIVITIES = [0.33, 0.008, 0.43]
LABEL_INFO = [(0, 'brain'), (1, 'skull'), (2, 'scalp')]
EXPECTED_TISSUE_ORDER = ['brain', 'skull', 'scalp']

# 10 micrometers: FEM mesh boundary nodes coincide near-exactly with the
# input BEM surface nodes it was built from.
BEM_MATCH_TOL = 1e-5


def find_subjects():
    cohort = sorted(
        os.path.basename(f).split('_snr')[0]
        for f in glob.glob(os.path.join(IDIP_RMS_DIR, 'sub-*_snr3_idip_rms.csv'))
    )
    return [s for s in cohort if s not in EXCLUDED_SUBJECTS]


def find_headmodel_file(subject):
    # The condition folder holding the head model varies per subject
    # (whichever preprocessed resting recording build_duneuro_forward_
    # models.m happened to attach it to), so search rather than guess the
    # name. A subject with only headmodel_surf_duneuro_openmeeg.mat (no
    # plain headmodel_surf_duneuro.mat) has no usable volumetric FEM
    # solution -- see EXCLUDED_SUBJECTS.
    matches = glob.glob(os.path.join(OMEGA_DATA_DIR, subject, '*', 'headmodel_surf_duneuro.mat'))
    if not matches:
        return None
    if len(matches) > 1:
        raise ValueError(f'{subject}: multiple headmodel_surf_duneuro.mat found: {matches}')
    return matches[0]


def classify_numeric_tissue_labels(subject, fem):
    """Reclassify iso2mesh-2026's raw numeric TissueLabels (['1','2','3'])
    as brain/skull/scalp by proximity to that subject's own BEM surfaces.
    Mirrors build_duneuro_forward_models.m's ClassifyFemTissueLabels()."""
    def load_bem(pattern):
        matches = glob.glob(os.path.join(OMEGA_ANAT_DIR, subject, pattern))
        if len(matches) != 1:
            raise ValueError(f'{subject}: expected exactly one {pattern}, found {matches}')
        return sio.loadmat(matches[0], simplify_cells=True)['Vertices']

    tree_scalp = cKDTree(load_bem('tess_head_bem_*V.mat'))
    tree_outer = cKDTree(load_bem('tess_outerskull_bem_*V.mat'))
    tree_inner = cKDTree(load_bem('tess_innerskull_bem_*V.mat'))

    vertices = fem['Vertices']
    elements = fem['Elements'].astype(np.int64) - 1
    tissue = fem['Tissue'].astype(np.int64)
    unique_labels = np.unique(tissue)
    if not np.array_equal(unique_labels, np.arange(1, len(unique_labels) + 1)):
        raise ValueError(f'{subject}: unexpected FEM Tissue values {unique_labels} (expected 1..N, no gaps).')

    names = {}
    for lbl in unique_labels:
        pts = vertices[np.unique(elements[tissue == lbl].flatten())]
        touches_scalp = np.mean(tree_scalp.query(pts)[0] < BEM_MATCH_TOL) > 0.01
        touches_outer = np.mean(tree_outer.query(pts)[0] < BEM_MATCH_TOL) > 0.01
        touches_inner = np.mean(tree_inner.query(pts)[0] < BEM_MATCH_TOL) > 0.01
        if touches_inner and not touches_outer and not touches_scalp:
            names[lbl] = 'brain'
        elif touches_inner and touches_outer and not touches_scalp:
            names[lbl] = 'skull'
        elif touches_outer and touches_scalp and not touches_inner:
            names[lbl] = 'scalp'
        else:
            raise ValueError(
                f'{subject}: could not classify tissue label {lbl} (touches: '
                f'scalp={touches_scalp} outer={touches_outer} inner={touches_inner}).')
    # order by the numeric label so the result lines up with `tissue` as-is
    return [names[lbl] for lbl in sorted(names)]


def convert_subject(subject, output_dir):
    headmodel_file = find_headmodel_file(subject)
    if headmodel_file is None:
        return None, 'no headmodel_surf_duneuro.mat (only an OpenMEEG fallback, if anything)'
    headmodel = sio.loadmat(headmodel_file, simplify_cells=True)

    fem_mesh_path = os.path.join(OMEGA_ANAT_DIR, headmodel['FemMeshFile'])
    fem = sio.loadmat(fem_mesh_path, simplify_cells=True)

    tissue_labels = list(fem['TissueLabels'])
    if all(label.isdigit() for label in tissue_labels):
        tissue_labels = classify_numeric_tissue_labels(subject, fem)
    if sorted(tissue_labels) != sorted(EXPECTED_TISSUE_ORDER):
        raise ValueError(
            f"{subject}: FEM mesh TissueLabels {tissue_labels} don't match the expected "
            f"set {EXPECTED_TISSUE_ORDER}.")

    nodes = np.asarray(fem['Vertices'], dtype=np.float64)  # already in meters
    elements = np.asarray(fem['Elements'], dtype=np.int64) - 1  # MATLAB 1-based -> 0-based

    # tissue_labels[i] names raw (1-based) Tissue value i+1 -- this ORDER
    # varies per subject (confirmed: e.g. sub-0010 is ['skull','scalp','brain']),
    # so raw Tissue values can't be assumed to already be brain=1/skull=2/
    # scalp=3. Remap explicitly to the canonical 0=brain/1=skull/2=scalp
    # index LABEL_INFO/CONDUCTIVITIES below assume, rather than blindly
    # doing `Tissue - 1`.
    canonical_index = {name: i for i, name in LABEL_INFO}
    raw_to_canonical = np.array([canonical_index[name] for name in tissue_labels])
    tissue_raw = np.asarray(fem['Tissue'], dtype=np.int64).flatten()  # 1-based
    labels = raw_to_canonical[tissue_raw - 1]

    out_file = os.path.join(output_dir, f'{subject}_vc_3layer.npz')
    np.savez(out_file, nodes=nodes, elements=elements, labels=labels,
             label_info=LABEL_INFO, conductivities=CONDUCTIVITIES)
    return out_file, None


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    subject_list = find_subjects()
    if not subject_list:
        raise FileNotFoundError(f'No subject folders found in {OMEGA_ANAT_DIR}')

    print(f'Found {len(subject_list)} candidate subject(s) (excluding {sorted(EXCLUDED_SUBJECTS)}): {subject_list}')
    skipped = []
    for subject in subject_list:
        out_file = os.path.join(OUTPUT_DIR, f'{subject}_vc_3layer.npz')
        if os.path.exists(out_file):
            print(f'Skipping {subject}: already converted.')
            continue
        print(f'Converting {subject}...')
        try:
            result, reason = convert_subject(subject, OUTPUT_DIR)
            if result is None:
                print(f'  SKIPPED {subject}: {reason}')
                skipped.append(subject)
            else:
                print(f'  wrote {result}')
        except Exception as exc:
            print(f'  FAILED {subject}: {exc}')
            skipped.append(subject)

    if skipped:
        print(f'\n{len(skipped)} subject(s) not converted: {skipped}')


if __name__ == '__main__':
    main()
