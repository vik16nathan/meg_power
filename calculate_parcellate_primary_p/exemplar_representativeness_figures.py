#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supplementary Figure 5: do the 4 secondary-P exemplar subjects' primary
current dipoles (5A) and primary dipole P (5B) look like the remaining 26
subjects', in both magnitude and spatial distribution? Each panel is a
histogram of per-vertex magnitude (pooled within each group, density-
normalized since group sizes differ ~6.5x) + a scatter of the two groups'
parcellated (Schaefer-600) average maps against each other, with Spearman r.

Also prints the whole-brain-total summary stats (mean AND median -- the
distribution across subjects is right-skewed by a few outliers, so these can
diverge substantially; median is arguably the fairer comparison given the
exemplars were selected to be near the cohort MEDIAN of primary dipole P)
that feed the Results-text sentence quoting these numbers.
"""
import glob
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.stats import spearmanr

OUTPUTS_DIR = '/export02/data/vikramn/hbm_manuscript_code/outputs'
PRIMARY_DIR = os.path.join(OUTPUTS_DIR, 'primary_p', 'fem_column_power')
PARCEL_DIR = os.path.join(OUTPUTS_DIR, 'primary_p', 'cent_surf_fem_fwd')
VC_DIR = os.path.join(OUTPUTS_DIR, 'secondary_p', 'vc_3layer')
FIG_DIR = os.path.join(OUTPUTS_DIR, 'primary_p', 'figures')

SNR = 3
EXEMPLARS = ['sub-0019', 'sub-0023', 'sub-0031', 'sub-0025']
EXEMPLAR_COLOR = '#4C72B0'
REMAINING_COLOR = '#888888'


def find_eligible_subjects():
    # same cohort scoping as secondary_p/convert_fem_mat_to_vc.py: only
    # subjects with a usable FEM head model are ever candidates for the
    # secondary-P exemplar pool this figure is validating
    return sorted(
        os.path.basename(f)[:-len('_vc_3layer.npz')]
        for f in glob.glob(os.path.join(VC_DIR, 'sub-*_vc_3layer.npz'))
    )


def load_vertexwise(subject):
    d = sio.loadmat(os.path.join(PRIMARY_DIR, f'{subject}_snr{SNR}_dip_p_i_rms.mat'), simplify_cells=True)
    return d['dip_p'], d['dip_i']


def load_parcellated(subject, kind):
    f = os.path.join(PARCEL_DIR, f's600_{subject}_snr{SNR}_{kind}_pri_2min_rest.csv')
    d = pd.read_csv(f)
    return dict(zip(d['region'], d['value']))


def group_average_map(subjects, kind):
    maps = [load_parcellated(s, kind) for s in subjects]
    common_regions = sorted(set.intersection(*(set(m) for m in maps)))
    return {r: np.mean([m[r] for m in maps]) for r in common_regions}


def print_summary(kind, label):
    eligible = find_eligible_subjects()
    remaining = [s for s in eligible if s not in EXEMPLARS]
    rows = []
    for s in eligible:
        p, i = load_vertexwise(s)
        value = p.sum() if kind == 'p' else i.sum()
        rows.append({'subject': s, 'total': value, 'group': 'exemplar' if s in EXEMPLARS else 'remaining'})
    df = pd.DataFrame(rows)
    ex, rem = df[df.group == 'exemplar']['total'], df[df.group == 'remaining']['total']
    print(f'{label} whole-brain total (n_exemplar=4, n_remaining={len(remaining)}):')
    print(f'  MEAN:   exemplar={ex.mean():.3e}  remaining={rem.mean():.3e}  ratio={rem.mean()/ex.mean():.2f}x')
    print(f'  MEDIAN: exemplar={ex.median():.3e}  remaining={rem.median():.3e}  ratio={rem.median()/ex.median():.2f}x')
    return df


def make_panel(kind, label, unit, out_name, panel_letter):
    eligible = find_eligible_subjects()
    remaining = [s for s in eligible if s not in EXEMPLARS]

    ex_vertex = np.concatenate([load_vertexwise(s)[0 if kind == 'p' else 1] for s in EXEMPLARS])
    rem_vertex = np.concatenate([load_vertexwise(s)[0 if kind == 'p' else 1] for s in remaining])

    ex_map = group_average_map(EXEMPLARS, kind)
    rem_map = group_average_map(remaining, kind)
    common = sorted(set(ex_map) & set(rem_map))
    ex_parcel = np.array([ex_map[r] for r in common])
    rem_parcel = np.array([rem_map[r] for r in common])
    r, p_val = spearmanr(rem_parcel, ex_parcel)

    fig, (ax_hist, ax_scatter) = plt.subplots(1, 2, figsize=(13, 5.5))

    bins = np.logspace(np.log10(min(ex_vertex.min(), rem_vertex.min())),
                        np.log10(max(ex_vertex.max(), rem_vertex.max())), 60)
    ax_hist.hist(rem_vertex, bins=bins, color=REMAINING_COLOR, alpha=0.6, density=True,
                 label=f'Remaining (n={len(remaining)})')
    ax_hist.hist(ex_vertex, bins=bins, color=EXEMPLAR_COLOR, alpha=0.6, density=True,
                 label='Exemplars (n=4)')
    ax_hist.set_xscale('log')
    ax_hist.set_xlabel(f'{label} ({unit})')
    ax_hist.set_ylabel('Density (per vertex, pooled within group)')
    ax_hist.set_title(f'{label} magnitude distribution')
    ax_hist.legend(fontsize=10)
    for spine in ('top', 'right'):
        ax_hist.spines[spine].set_visible(False)

    ax_scatter.scatter(rem_parcel, ex_parcel, s=14, color=EXEMPLAR_COLOR, alpha=0.6, linewidths=0)
    lims = [min(rem_parcel.min(), ex_parcel.min()), max(rem_parcel.max(), ex_parcel.max())]
    ax_scatter.plot(lims, lims, color='#888888', lw=1, ls='--', zorder=0, label='y = x')
    ax_scatter.set_xscale('log')
    ax_scatter.set_yscale('log')
    ax_scatter.set_xlabel(f'Remaining subjects, parcellated avg. (n={len(remaining)})')
    ax_scatter.set_ylabel('Exemplars, parcellated avg. (n=4)')
    ax_scatter.set_title(f'Spatial distribution (Spearman r={r:.3f}, n={len(common)} regions)')
    ax_scatter.legend(fontsize=10, loc='upper left')
    for spine in ('top', 'right'):
        ax_scatter.spines[spine].set_visible(False)

    fig.suptitle(f'Supplementary Figure {panel_letter}: {label}, exemplars vs. remaining subjects', fontsize=14)
    fig.tight_layout()
    out_file = os.path.join(FIG_DIR, out_name)
    fig.savefig(out_file, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {out_file}')
    return r, p_val


def main():
    os.makedirs(FIG_DIR, exist_ok=True)

    print('=== Primary current dipole (Supp. Fig. 5A) ===')
    print_summary('i', 'Primary current dipole')
    r_i, p_i = make_panel('i', 'Primary current dipole', 'A·m', 'supp_fig5a_current_dipole_representativeness.png', '5A')
    print(f'  Spearman r={r_i:.4f}  p={p_i:.3e}')

    print()
    print('=== Primary dipole P (Supp. Fig. 5B) ===')
    print_summary('p', 'Primary dipole P')
    r_p, p_p = make_panel('p', 'Primary dipole P', 'W', 'supp_fig5b_primary_p_representativeness.png', '5B')
    print(f'  Spearman r={r_p:.4f}  p={p_p:.3e}')


if __name__ == '__main__':
    main()
