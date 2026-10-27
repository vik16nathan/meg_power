#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Single-dipole SEP validation (three-layer FEM): places one simulated 20 nA.m
current dipole at the postcentral gyrus (somatosensory cortex) of sub-0002's
three-layer FEM mesh, computes the resulting secondary current density with
DUNEuro, and breaks its power/volume/current down by tissue (brain/skull/
scalp) -- the three-layer counterpart to the manuscript's six-layer SEP
comparison (Figure 4/5): "simulating the currents induced by a single
somatosensory evoked potential (SEP) of 20 nA.m located on the postcentral
gyrus."

Dipole location/moment: taken directly from data_from_malte/
dipole_simulation_results_sep.npz's own dipole_position/dipole_moment --
the original, manually-selected-within-Brainstorm postcentral-gyrus
position used for this same validation previously (NOT re-derived from an
atlas scout's seed vertex, which only approximates the intended location).

Must be run with a Python 3.8 interpreter matching the duneuropy.so build
(e.g. /usr/bin/python3.8), same as eeg_dipole_power_computation.py.
"""
import os
import sys

import numpy as np

DUNEUROPY_PATH = '/home/bic/vikramn/build-release/duneuro-py/src/'
sys.path.append(DUNEUROPY_PATH)
import duneuropy as dp  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eeg_dipole_power_computation import build_driver_cfg  # noqa: E402

SUBJECT = 'sub-0002'
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
VC_DIR = os.path.join(SCRIPT_DIR, '..', 'outputs', 'secondary_p', 'vc_3layer')
# Position/moment precomputed by a separate one-off (see git history/README):
# system python3.8's scipy (1.3.3) predates loadmat's simplify_cells kwarg,
# so the Desikan-Killiany scout lookup runs under the regular python3 (modern
# scipy) instead and hands off just the tiny (position, moment) result here.
DIPOLE_FILE = os.path.join(SCRIPT_DIR, '..', 'outputs', 'secondary_p',
                            f'{SUBJECT}_postcentral_dipole.npz')
OUTPUT_FILE = os.path.join(SCRIPT_DIR, '..', 'outputs', 'secondary_p',
                            f'{SUBJECT}_sep_postcentral_20nAm_3layer.npz')


def get_postcentral_dipole():
    d = np.load(DIPOLE_FILE)
    return d['position'], d['moment'], int(d['vertex_idx'])


def calc_tetr_vol_from_vert(v_list):
    a = v_list.T
    b = np.row_stack((a, [1, 1, 1, 1]))
    return 1.0 / 6.0 * np.abs(np.linalg.det(b))


def main():
    vc_data = np.load(os.path.join(VC_DIR, f'{SUBJECT}_vc_3layer.npz'))
    nodes, elements, labels, conductivities = (
        vc_data['nodes'], vc_data['elements'], vc_data['labels'], vc_data['conductivities'])

    position, moment, _ = get_postcentral_dipole()
    print(f'{SUBJECT}: position={position} (from dipole_simulation_results_sep.npz), '
          f'moment={moment} (|moment|={np.linalg.norm(moment):.3e} A.m)')

    driver_cfg = build_driver_cfg(vc_data, 'partial_integration')  # deliberate choice for secondary-P, distinct from Brainstorm's own venant-based forward modeling -- see eeg_dipole_power_computation.py's SOURCE_MODEL comment
    driver = dp.MEEGDriver3d(driver_cfg)
    dipole = dp.Dipole3d(position.tolist(), moment.tolist())
    solution_storage = driver.makeDomainFunction()
    driver.solveEEGForward(dipole, solution_storage, driver_cfg)
    result = driver.exportVolumeConductorAndFunction(solution_storage)
    current_density = np.asarray(result['currentDensityAtElementCenters'])

    el_xyz = nodes[elements]
    v_els = np.array([calc_tetr_vol_from_vert(el_xyz[i]) for i in range(elements.shape[0])])
    conds = conductivities[labels]
    i_sq_els = np.linalg.norm(current_density, axis=1) ** 2
    p_els = i_sq_els * v_els / conds

    np.savez(OUTPUT_FILE, v_els=v_els, i_sq_els=i_sq_els, p_els=p_els, labels=labels,
              dipole_position=position, dipole_moment=moment)
    print(f'wrote {OUTPUT_FILE}')
    print(f'total power = {p_els.sum():.3e} W, total volume = {v_els.sum():.3e} m^3')
    for t, name in [(0, 'brain'), (1, 'skull'), (2, 'scalp')]:
        mask = labels == t
        print(f'  {name}: {p_els[mask].sum()/p_els.sum()*100:.1f}% power, '
              f'{v_els[mask].sum()/v_els.sum()*100:.1f}% volume')


if __name__ == '__main__':
    main()
