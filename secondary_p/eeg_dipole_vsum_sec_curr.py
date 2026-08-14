#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GOAL: For each of a few exemplar subjects, compute the vector sum (across
all ~15k cortical primary-current dipoles) of the secondary current density
at every FEM element, then use it to get power/current/volume broken down
by tissue (brain/skull/scalp) -- the inputs generate_figures.py needs for
the tissue-proportion pie charts.

Must be run with a Python 3.8 interpreter matching the duneuropy.so build
(e.g. /usr/bin/python3.8), same as eeg_dipole_power_computation.py.
"""
import glob
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

DUNEUROPY_PATH = '/home/bic/vikramn/build-release/duneuro-py/src/'
sys.path.append(DUNEUROPY_PATH)
import duneuropy as dp  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eeg_dipole_power_computation import (  # noqa: E402
    SNR, VC_DIR, VERTICES_DIR, IDIP_RMS_DIR, N_WORKERS, build_driver_cfg,
)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, '..', 'outputs', 'secondary_p', 'vsum_sec_curr')

# Exemplar subjects for the individual/supplemental pie-chart panels --
# generate_figures.py averages over every subject this script actually
# produces output for, not just these four, but these are the ones picked
# out for per-subject display alongside the average.
# sub-0022 -> sub-0023: swapped for representativeness (sub-0022's total
# primary P/I_dip were ~208%/~53% above the n=26-remaining-subject mean --
# see exemplar_sensitivity_supp_fig.R) and sex balance (sub-0023 is F,
# making this 2M/2F instead of 3M/1F). sub-0023's total P sits almost
# exactly at the cohort median (0.0% deviation) and I_dip -2.7%.
EXEMPLAR_SUBJECTS = ['sub-0019', 'sub-0023', 'sub-0031', 'sub-0025']

_worker_driver = None
_worker_driver_cfg = None


def _init_worker(driver_cfg):
    global _worker_driver, _worker_driver_cfg
    _worker_driver_cfg = driver_cfg
    _worker_driver = dp.MEEGDriver3d(driver_cfg)


def _sec_curr_density(position_and_moment):
    position, moment = position_and_moment
    dipole = dp.Dipole3d(position, moment)
    solution_storage = _worker_driver.makeDomainFunction()
    _worker_driver.solveEEGForward(dipole, solution_storage, _worker_driver_cfg)
    result = _worker_driver.exportVolumeConductorAndFunction(solution_storage)
    return np.asarray(result['currentDensityAtElementCenters'])


def calc_tetr_vol_from_vert(v_list):
    """Volume of a tetrahedral element from its 4 ordered 3D vertex coordinates."""
    a = v_list.T
    b = np.row_stack((a, [1, 1, 1, 1]))
    return 1.0 / 6.0 * np.abs(np.linalg.det(b))


def calculate_fem_p(nodes, elements, labels, conductivities, current_density):
    all_el_xyz = np.asarray(nodes)[elements, :]
    v_els = np.array([calc_tetr_vol_from_vert(all_el_xyz[i]) for i in range(elements.shape[0])])
    conds = np.asarray(conductivities)[labels]
    i_sq_els = np.linalg.norm(current_density, axis=1) ** 2
    p_els = i_sq_els * v_els / conds
    return v_els, i_sq_els, p_els


def find_subjects():
    # see eeg_dipole_power_computation.py's find_subjects() for what this does
    only = os.environ.get('SECONDARY_P_SUBJECTS')
    only = set(only.split(',')) if only else None

    subjects = []
    for vc_file in sorted(glob.glob(os.path.join(VC_DIR, 'sub-*_vc_3layer.npz'))):
        subject = os.path.basename(vc_file)[:-len('_vc_3layer.npz')]
        if only is not None and subject not in only:
            continue
        vertices_file = os.path.join(VERTICES_DIR, f'{subject}_vertices.csv')
        idip_file = os.path.join(IDIP_RMS_DIR, f'{subject}_snr{SNR}_idip_rms.csv')
        if os.path.exists(vertices_file) and os.path.exists(idip_file):
            subjects.append(subject)
    return subjects


def process_subject(subject):
    vsum_file = os.path.join(OUTPUT_DIR, f'{subject}_summed_sec_curr.npz')
    p_file = os.path.join(OUTPUT_DIR, f'{subject}_p_output.npz')
    if os.path.exists(vsum_file) and os.path.exists(p_file):
        print(f'Skipping {subject}: outputs already exist.')
        return

    print(f'=== {subject}: loading inputs ===')
    vc_data = np.load(os.path.join(VC_DIR, f'{subject}_vc_3layer.npz'))
    nodes = vc_data['nodes']
    elements = vc_data['elements']
    labels = vc_data['labels']
    conductivities = vc_data['conductivities']
    label_info = vc_data['label_info']

    dipole_positions = np.array(pd.read_csv(
        os.path.join(VERTICES_DIR, f'{subject}_vertices.csv'), header=None))
    dipole_moments = np.array(pd.read_csv(
        os.path.join(IDIP_RMS_DIR, f'{subject}_snr{SNR}_idip_rms.csv'), header=None))
    if len(dipole_positions) != len(dipole_moments):
        raise ValueError(
            f'{subject}: {len(dipole_positions)} dipole positions vs. '
            f'{len(dipole_moments)} dipole moments -- must match.')

    driver_cfg = build_driver_cfg(vc_data, 'partial_integration')  # deliberate choice for secondary-P, distinct from Brainstorm's own venant-based forward modeling -- see eeg_dipole_power_computation.py's SOURCE_MODEL comment
    args = list(zip(dipole_positions.tolist(), dipole_moments.tolist()))

    print(f'{subject}: summing secondary current density over {len(args)} dipoles '
          f'with {N_WORKERS} workers...')
    start = time.perf_counter()
    summed = np.zeros((elements.shape[0], 3))
    with Pool(N_WORKERS, initializer=_init_worker, initargs=(driver_cfg,)) as pool:
        for i, current_density in enumerate(pool.imap_unordered(_sec_curr_density, args, chunksize=4)):
            summed += current_density
            if (i + 1) % 1000 == 0:
                elapsed = time.perf_counter() - start
                print(f'  {subject}: {i + 1}/{len(args)} dipoles done ({elapsed:.0f}s elapsed)')
    elapsed = time.perf_counter() - start
    print(f'{subject}: done in {elapsed:.1f}s ({elapsed / len(args) * 1000:.1f} ms/dipole)')

    np.savez(vsum_file, nodes=nodes, elements=elements, labels=labels,
             conductivities=conductivities, label_info=label_info, current_density=summed)
    print(f'  wrote {vsum_file}')

    v_els, i_sq_els, p_els = calculate_fem_p(nodes, elements, labels, conductivities, summed)
    np.savez(p_file, v_els=v_els, i_sq_els=i_sq_els, p_els=p_els)
    print(f'  wrote {p_file}')


def main():
    # This computation (per-dipole full-mesh current-density export, not just
    # a scalar power) is markedly more expensive than eeg_dipole_power_
    # computation.py, so it's deliberately scoped to just the exemplar
    # subjects generate_figures.py shows individually -- not run for every
    # subject the way eeg_dipole_power_computation.py is. Restrict further
    # to one subject via SECONDARY_P_SUBJECTS to run exemplars concurrently
    # as separate processes (with SECONDARY_P_N_WORKERS lowered per-process
    # so they share the core count instead of each grabbing 20).
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    available = find_subjects()
    subject_list = [s for s in EXEMPLAR_SUBJECTS if s in available]
    missing = [s for s in EXEMPLAR_SUBJECTS if s not in available]
    if missing:
        print(f'Note: exemplar subject(s) missing required inputs, skipping: {missing}')

    if not subject_list:
        raise FileNotFoundError('No subject has the inputs needed for vector-sum secondary current computation.')

    print(f'Processing {len(subject_list)} subject(s): {subject_list}')
    for subject in subject_list:
        process_subject(subject)


if __name__ == '__main__':
    main()
