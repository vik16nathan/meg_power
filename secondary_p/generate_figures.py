#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates the secondary-current-power figures for whichever subjects have
been through the pipeline so far (convert_fem_mat_to_vc.py ->
eeg_dipole_power_computation.py -> eeg_dipole_vsum_sec_curr.py), saving
everything under figures/. Replaces SecondaryCurrentsV3.ipynb, which mixed
this real per-subject analysis in with one-off exploratory cells against a
collaborator's toy sphere/realistic head models on a local Windows path
that no longer exists in this environment.

Two families of figures:
  - Tissue proportions (brain/skull/scalp) of FEM element Volume, and of
    Power/Current from the vector-summed secondary current density
    (eeg_dipole_vsum_sec_curr.py's output): one averaged-across-subjects
    pie chart (with SD/SEM reported) as the main result, plus a 2x2 grid
    of individual exemplar subjects as a supplement.
  - Whole-brain total secondary power per subject (eeg_dipole_power_
    computation.py's per-dipole totals, summed), matching the primary-P
    histogram convention already used in calculate_parcellate_primary_p/
    plot_total_p_histogram.py.
"""
import glob
import os
from functools import lru_cache

import matplotlib
matplotlib.use('Agg')  # headless: no display, just save to file
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUTS_DIR = os.path.join(SCRIPT_DIR, '..', 'outputs')
VC_DIR = os.path.join(OUTPUTS_DIR, 'secondary_p', 'vc_3layer')
FEM_P_RMS_DIR = os.path.join(OUTPUTS_DIR, 'secondary_p', 'fem_p_rms')
VSUM_DIR = os.path.join(OUTPUTS_DIR, 'secondary_p', 'vsum_sec_curr')
PRIMARY_P_DIR = os.path.join(OUTPUTS_DIR, 'primary_p', 'fem_column_power')
FIG_DIR = os.path.join(SCRIPT_DIR, 'figures')

# Two distinct quantities per the paper (Eq. 6 vs Eq. 7) -- NOT interchangeable:
#   secondary dipole P  = sum of each dipole's OWN power, computed independently
#                          (eeg_dipole_power_computation.py's fem_p_rms.csv, summed)
#   total secondary P   = power from the VECTOR SUM of all dipoles' fields first
#                          (eeg_dipole_vsum_sec_curr.py's p_els, summed)
# total secondary P is generally larger (captures constructive interference
# between spatially-correlated dipoles' fields) and is only computed for the
# 4 exemplar subjects (expensive); secondary dipole P scales to the full cohort.

SNR = 3
EXEMPLAR_SUBJECTS = ['sub-0019', 'sub-0023', 'sub-0031', 'sub-0025']

# Same subject flagged in calculate_parcellate_primary_p/plot_total_p_histogram.py
# (elevated whole-brain MNE current dipole moment, R/thickness normal --
# flagged, not root-caused) -- consistent with that, its secondary dipole P
# is ~20-30x every other subject too. Excluded from summary figures/stats
# below, but never silently dropped from per-subject listings.
OUTLIER_SUBJECT = 'sub-0011'

TISSUE_NAMES = {0: 'brain', 1: 'skull', 2: 'scalp'}
TISSUE_COLOR = {'brain': 'gray', 'skull': 'lavender', 'scalp': 'brown'}
TISSUE_ORDER = ['brain', 'skull', 'scalp']


def tetrahedron_volumes(el_xyz):
    """Volume of every tetrahedral element at once, given (n_elements, 4, 3) vertex coordinates."""
    a = el_xyz[:, 0] - el_xyz[:, 3]
    b = el_xyz[:, 1] - el_xyz[:, 3]
    c = el_xyz[:, 2] - el_xyz[:, 3]
    return np.abs(np.einsum('ij,ij->i', a, np.cross(b, c))) / 6.0


def tissue_proportions(labels, values):
    """Fraction of `values` (one entry per FEM element) attributable to each tissue."""
    sums = {name: 0.0 for name in TISSUE_ORDER}
    for label_idx, value in zip(labels, values):
        sums[TISSUE_NAMES[label_idx]] += value
    total = sum(sums.values())
    return {name: sums[name] / total * 100 for name in TISSUE_ORDER}


def plot_pie(ax, percentages, title):
    labels = TISSUE_ORDER
    values = [percentages[label] for label in labels]
    colors = [TISSUE_COLOR[label] for label in labels]
    explode = [0.1 if v < 5 else 0 for v in values]
    wedges, _, autotexts = ax.pie(
        values, labels=None, autopct='%1.1f%%', startangle=140, colors=colors,
        textprops={'fontsize': 14}, explode=explode, pctdistance=0.6, radius=1.0,
    )
    for autotext in autotexts:
        autotext.set_fontsize(16)
        autotext.set_bbox(dict(facecolor='white', edgecolor='none', alpha=0.6, pad=2))
    ax.set_title(title, fontsize=18, pad=10)
    return wedges


@lru_cache(maxsize=None)
def load_subject_volume_conductor(subject):
    return np.load(os.path.join(VC_DIR, f'{subject}_vc_3layer.npz'))


@lru_cache(maxsize=None)
def load_subject_vsum_power(subject):
    p_file = os.path.join(VSUM_DIR, f'{subject}_p_output.npz')
    if not os.path.exists(p_file):
        return None
    return np.load(p_file)


@lru_cache(maxsize=None)
def volume_values(subject):
    vc = load_subject_volume_conductor(subject)
    nodes, elements, labels = vc['nodes'], vc['elements'], vc['labels']
    v_els = tetrahedron_volumes(nodes[elements])
    return labels, v_els


def find_subjects_with_vc():
    return sorted(
        os.path.basename(f)[:-len('_vc_3layer.npz')]
        for f in glob.glob(os.path.join(VC_DIR, 'sub-*_vc_3layer.npz'))
    )


def find_subjects_with_vsum():
    return sorted(
        os.path.basename(f)[:-len('_p_output.npz')]
        for f in glob.glob(os.path.join(VSUM_DIR, 'sub-*_p_output.npz'))
    )


def make_average_pie_figure(subjects, variable, load_values_fn, out_name):
    """One pie chart of mean tissue proportions across `subjects`, with the
    across-subject SD reported in the legend for each slice."""
    per_subject = []
    for subject in subjects:
        labels, values = load_values_fn(subject)
        if labels is None:
            continue
        per_subject.append(tissue_proportions(labels, values))
    if not per_subject:
        print(f'  no subjects with data for "{variable}" -- skipping {out_name}')
        return

    means = {t: np.mean([p[t] for p in per_subject]) for t in TISSUE_ORDER}
    sds = {t: np.std([p[t] for p in per_subject], ddof=1) if len(per_subject) > 1 else 0.0
           for t in TISSUE_ORDER}

    fig, ax = plt.subplots(figsize=(8, 8))
    plot_pie(ax, means, f'{variable} by Tissue\n(mean across n={len(per_subject)} subjects)')
    legend_labels = [f'{t}: {means[t]:.1f}% ± {sds[t]:.1f}% (SD)' for t in TISSUE_ORDER]
    ax.legend(legend_labels, title='Tissue', loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)
    fig.tight_layout()
    out_file = os.path.join(FIG_DIR, out_name)
    fig.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {out_file} (n={len(per_subject)} subjects)')


def make_individual_pie_grid(subjects, variable, load_values_fn, out_name):
    subjects = [s for s in subjects if load_values_fn(s)[0] is not None]
    if not subjects:
        print(f'  no exemplar subjects with data for "{variable}" -- skipping {out_name}')
        return

    ncols = 2
    nrows = int(np.ceil(len(subjects) / ncols))
    fig, axs = plt.subplots(nrows, ncols, figsize=(7 * ncols, 7 * nrows))
    axs = np.atleast_1d(axs).flatten()

    wedges_for_legend = None
    for ax, subject in zip(axs, subjects):
        labels, values = load_values_fn(subject)
        percentages = tissue_proportions(labels, values)
        wedges = plot_pie(ax, percentages, subject)
        wedges_for_legend = wedges_for_legend or wedges
    for ax in axs[len(subjects):]:
        ax.axis('off')

    fig.suptitle(f'{variable} by Tissue (exemplar subjects)', fontsize=20)
    if wedges_for_legend is not None:
        fig.legend(wedges_for_legend, TISSUE_ORDER, title='Tissue', loc='lower center',
                   ncol=len(TISSUE_ORDER), fontsize=14, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    out_file = os.path.join(FIG_DIR, out_name)
    fig.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {out_file} (subjects={subjects})')


def power_values(subject):
    vc = load_subject_volume_conductor(subject)
    p_output = load_subject_vsum_power(subject)
    if p_output is None:
        return None, None
    return vc['labels'], p_output['p_els']


def current_values(subject):
    vc = load_subject_volume_conductor(subject)
    p_output = load_subject_vsum_power(subject)
    if p_output is None:
        return None, None
    return vc['labels'], p_output['i_sq_els']


# single_dipole_sep_validation.py's output: ONE simulated 20 nA.m dipole at
# sub-0002's postcentral gyrus (three-layer FEM), not real per-subject data
# -- the three-layer counterpart to the manuscript's six-layer single-dipole
# SEP comparison (Figure 4/5), for a like-for-like check (the existing text
# instead compares the six-layer single-dipole result against the three-
# layer n=4 REAL-PARTICIPANT average, which mixes a simulated dipole against
# real MEG data).
SEP_VALIDATION_FILE = os.path.join(OUTPUTS_DIR, 'secondary_p', 'sub-0002_sep_postcentral_20nAm_3layer.npz')


def make_sep_validation_pie_charts():
    if not os.path.exists(SEP_VALIDATION_FILE):
        print(f'  {SEP_VALIDATION_FILE} not found -- run single_dipole_sep_validation.py first, skipping.')
        return
    sep = np.load(SEP_VALIDATION_FILE)
    labels = sep['labels']
    panels = [('Volume', sep['v_els']), ('Power Dissipation', sep['p_els']),
              ('Squared Secondary Current', sep['i_sq_els'])]

    fig, axs = plt.subplots(1, 3, figsize=(21, 7))
    wedges_for_legend = None
    for ax, (variable, values) in zip(axs, panels):
        percentages = tissue_proportions(labels, values)
        wedges = plot_pie(ax, percentages, variable)
        wedges_for_legend = wedges_for_legend or wedges

    fig.suptitle('Single-dipole SEP validation (20 nA.m, postcentral gyrus, sub-0002 three-layer FEM)', fontsize=16)
    fig.legend(wedges_for_legend, TISSUE_ORDER, title='Tissue', loc='lower center',
               ncol=len(TISSUE_ORDER), fontsize=14, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    out_file = os.path.join(FIG_DIR, 'fig_sep_validation_3layer_by_tissue.png')
    fig.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {out_file}')


def make_secondary_dipole_p_histogram():
    """Whole-brain secondary dipole P per subject (Eq. 7: sum of each dipole's
    own power) -- NOT total secondary P (Eq. 6), which needs the vsum pipeline."""
    files = sorted(glob.glob(os.path.join(FEM_P_RMS_DIR, f'sub-*_snr{SNR}_fem_p_rms.csv')))
    if not files:
        print(f'  no {FEM_P_RMS_DIR}/sub-*_snr{SNR}_fem_p_rms.csv files -- '
              f'run eeg_dipole_power_computation.py first, skipping histogram.')
        return

    rows = []
    for f in files:
        subject = os.path.basename(f).split('_')[0]
        total_p = pd.read_csv(f, header=None).to_numpy().sum()
        rows.append({'subject': subject, 'secondary_dipole_p_W': total_p})
    df = pd.DataFrame(rows).sort_values('subject').reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(df['secondary_dipole_p_W'], bins=max(5, len(df) // 2), color='#4C72B0', edgecolor='white')
    ax.set_xlabel('Whole-brain secondary dipole P (W)')
    ax.set_ylabel('Number of subjects')
    ax.set_title(f'Secondary dipole P across subjects (SnrFixed={SNR}, n={len(df)})')

    out_file = os.path.join(FIG_DIR, f'secondary_dipole_p_histogram_snr{SNR}.png')
    fig.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {out_file}')
    print(df.to_string(index=False))
    print_summary_excluding_outlier('secondary_dipole_p_W', df)


def print_summary_excluding_outlier(column, df):
    """Matches calculate_parcellate_primary_p/plot_total_p_histogram.py's
    convention: OUTLIER_SUBJECT stays in the plot/listing above (never
    silently dropped), but is excluded from the reported summary stats."""
    def stats(label, values):
        print(f'{label} (n={len(values)}): mean={values.mean():.2e}  median={values.median():.2e}  '
              f'min={values.min():.2e}  max={values.max():.2e}')

    stats('All subjects', df[column])
    if OUTLIER_SUBJECT in df['subject'].values:
        print(f'{OUTLIER_SUBJECT} flagged as an outlier (see OUTLIER_SUBJECT comment above) -- '
              f'value = {df.loc[df.subject == OUTLIER_SUBJECT, column].iloc[0]:.2e}, excluded below:')
        stats(f'Excluding {OUTLIER_SUBJECT}', df.loc[df['subject'] != OUTLIER_SUBJECT, column])


def print_exemplar_validation_stats():
    """Checks two manuscript claims against the 4 exemplar subjects
    (EXEMPLAR_SUBJECTS -- same set as Figure 3B/3D):
      1. Brain tissue's share of total secondary power ('the vast majority
         of the power ... is dissipated in the gray matter').
      2. Average secondary power across the 4 subjects, reported both ways
         (Eq. 6 total secondary P and Eq. 7 secondary dipole P) since the
         manuscript text doesn't specify which -- compare each against the
         stated 1.73e-9 W before deciding which one (if either) it means."""
    brain_pcts, total_secondary_p, dipole_secondary_p = [], [], []
    for subject in EXEMPLAR_SUBJECTS:
        labels, values = power_values(subject)
        p_output = load_subject_vsum_power(subject)
        _, sec = load_primary_secondary_pair(subject)
        if labels is None or p_output is None or sec is None:
            print(f'  {subject}: missing vsum and/or secondary dipole P data -- skipping')
            continue
        brain_pcts.append(tissue_proportions(labels, values)['brain'])
        total_secondary_p.append(float(p_output['p_els'].sum()))
        dipole_secondary_p.append(float(sec.sum()))

    n = len(brain_pcts)
    print(f'Brain tissue share of secondary power dissipation (n={n} exemplars):')
    for subject, pct in zip(EXEMPLAR_SUBJECTS, brain_pcts):
        print(f'  {subject}: {pct:.1f}%')
    print(f'  mean = {np.mean(brain_pcts):.1f}%  (manuscript states 93.6%)')
    print()
    print(f'Average secondary power across exemplars (n={n}), both definitions:')
    print(f'  mean total secondary P (Eq. 6)  = {np.mean(total_secondary_p):.3e} W')
    print(f'  mean secondary dipole P (Eq. 7) = {np.mean(dipole_secondary_p):.3e} W')
    print(f'  (manuscript states 1.73e-09 W -- neither matches exactly; '
          f'Eq. 7 is the closer of the two)')


def print_total_secondary_p_ratio(scope='exemplars'):
    """Total secondary P (Eq. 6, vsum pipeline) as a percentage of total
    primary P, per subject and averaged. `scope`: 'exemplars' (the 4
    representative subjects) or 'all' (every subject with total secondary P
    computed so far -- currently the same 4, since that's the pipeline's
    real bottleneck, but this will pick up more automatically if/when total
    secondary P is ever computed beyond the exemplars)."""
    if scope not in ('exemplars', 'all'):
        raise ValueError(f"scope must be 'exemplars' or 'all', got {scope!r}")

    vsum_subjects = find_subjects_with_vsum()
    subjects = [s for s in EXEMPLAR_SUBJECTS if s in vsum_subjects] if scope == 'exemplars' else vsum_subjects
    subjects = [s for s in subjects if s != OUTLIER_SUBJECT]

    rows = []
    for subject in subjects:
        p_output = load_subject_vsum_power(subject)
        prim_file = os.path.join(PRIMARY_P_DIR, f'{subject}_snr{SNR}_dip_p_i_rms.mat')
        if p_output is None or not os.path.exists(prim_file):
            continue
        total_secondary_p = float(p_output['p_els'].sum())
        total_primary_p = float(sio.loadmat(prim_file, simplify_cells=True)['dip_p'].sum())
        rows.append({'subject': subject, 'total_secondary_p_W': total_secondary_p,
                     'total_primary_p_W': total_primary_p,
                     'pct_secondary_of_primary': total_secondary_p / total_primary_p * 100})

    label = 'representative subjects' if scope == 'exemplars' else 'all subjects with total secondary P'
    if not rows:
        print(f'  no subjects with both total secondary P and primary P yet ({label}) -- skipping')
        return
    df = pd.DataFrame(rows)
    print(f'Total secondary P as % of total primary P ({label}, n={len(df)}):')
    print(df.to_string(index=False))
    print(f'  mean = {df["pct_secondary_of_primary"].mean():.1f}%  '
          f'range = [{df["pct_secondary_of_primary"].min():.1f}%, {df["pct_secondary_of_primary"].max():.1f}%]')


def load_primary_secondary_pair(subject):
    """Vertex-wise (primary dipole P, secondary dipole P) for one subject, or
    (None, None) if either input is missing or their vertex counts disagree."""
    prim_file = os.path.join(PRIMARY_P_DIR, f'{subject}_snr{SNR}_dip_p_i_rms.mat')
    sec_file = os.path.join(FEM_P_RMS_DIR, f'{subject}_snr{SNR}_fem_p_rms.csv')
    if not (os.path.exists(prim_file) and os.path.exists(sec_file)):
        return None, None
    prim = np.asarray(sio.loadmat(prim_file, simplify_cells=True)['dip_p'])
    sec = pd.read_csv(sec_file, header=None).to_numpy().flatten()
    if len(prim) != len(sec):
        print(f'  WARNING {subject}: {len(prim)} primary vertices vs. {len(sec)} '
              f'secondary vertices -- skipping (should always match).')
        return None, None
    return prim, sec


def _fig3_scope_subjects(scope):
    available = sorted(
        os.path.basename(f).split('_')[0]
        for f in glob.glob(os.path.join(FEM_P_RMS_DIR, f'sub-*_snr{SNR}_fem_p_rms.csv'))
    )
    # OUTLIER_SUBJECT excluded either way: it isn't one of the exemplars, and
    # in 'all' scope its ~20-30x-elevated magnitude would stretch the log-log
    # axes and correlation color normalization for every other subject.
    if scope == 'exemplars':
        return [s for s in EXEMPLAR_SUBJECTS if s in available and s != OUTLIER_SUBJECT]
    return [s for s in available if s != OUTLIER_SUBJECT]


def make_figure_3c_scatter(scope='exemplars', scatter_subsample=1500, seed=0):
    """Figure 3C: primary vs. secondary dipole P scatter (log-log) with one
    regression line per subject, colored by that subject's own Pearson r --
    the standalone scatter half of what Figure 3B's combined barplot+scatter
    layout shows, split out since it only needs secondary dipole P (Eq. 7),
    not total secondary P (Eq. 6/vsum), so it's ready independent of the
    (slower) vsum pipeline."""
    if scope not in ('exemplars', 'all'):
        raise ValueError(f"scope must be 'exemplars' or 'all', got {scope!r}")

    subjects = _fig3_scope_subjects(scope)
    per_subject = []
    for subject in subjects:
        prim, sec = load_primary_secondary_pair(subject)
        if prim is None:
            continue
        log_prim, log_sec = np.log10(prim), np.log10(sec)
        r = np.corrcoef(log_prim, log_sec)[0, 1]
        slope, intercept = np.polyfit(log_prim, log_sec, 1)
        per_subject.append(dict(subject=subject, prim=prim, sec=sec, r=r, slope=slope, intercept=intercept))
    if not per_subject:
        print('  no subjects with both primary and secondary dipole P yet -- skipping fig 3c')
        return
    per_subject.sort(key=lambda d: d['r'])

    cmap = plt.get_cmap('viridis')
    norm = plt.Normalize(vmin=per_subject[0]['r'], vmax=per_subject[-1]['r'])

    fig, ax = plt.subplots(figsize=(8, 6.5))
    rng = np.random.default_rng(seed)
    all_prim = np.concatenate([d['prim'] for d in per_subject])
    x_range = np.array([all_prim.min(), all_prim.max()])
    for d in per_subject:
        n = len(d['prim'])
        idx = rng.choice(n, size=min(scatter_subsample, n), replace=False)
        ax.scatter(d['prim'][idx], d['sec'][idx], s=3, color=cmap(norm(d['r'])),
                   alpha=0.3, linewidths=0, rasterized=True)
    for d in per_subject:
        y_fit = 10 ** (d['slope'] * np.log10(x_range) + d['intercept'])
        ax.plot(x_range, y_fit, color=cmap(norm(d['r'])), lw=1.3, alpha=0.9, label=d['subject'])

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Primary dipole P (W)')
    ax.set_ylabel('Secondary dipole P (W)')
    scope_label = 'representative subjects' if scope == 'exemplars' else 'all subjects'
    ax.set_title(f'Figure 3C: Primary vs. secondary dipole P ({scope_label}, n={len(per_subject)})', fontsize=12)
    ax.grid(True, which='major', alpha=0.15, lw=0.5)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    ax.legend(fontsize=8, loc='upper left', title='Subject (line)')

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.05, pad=0.02)
    cbar.set_label('Per-subject Pearson r', fontsize=10)

    fig.tight_layout()
    out_file = os.path.join(FIG_DIR, f'fig3c_scatter_primary_vs_secondary_dipole_p_snr{SNR}_{scope}.png')
    fig.savefig(out_file, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {out_file} (n={len(per_subject)} subjects, '
          f'r range [{per_subject[0]["r"]:.3f}, {per_subject[-1]["r"]:.3f}])')


def make_figure_3b(scope='exemplars', scatter_subsample=1500, seed=0):
    """Figure 3B: per-subject correlation (barplot) and primary-vs-secondary
    dipole P (scatter, log-log). `scope` toggles which subjects are shown:
    'exemplars' -- the 4 representative subjects, matching the other Figure 3
    panels' exemplar-based reporting; 'all' -- every subject with secondary
    dipole P computed so far (grows automatically as the cohort job
    finishes). ~15k points/subject is too dense for a legible per-individual
    color scheme, so: scatter points are subsampled and drawn at very low
    alpha (density only, not identity), and each subject's regression line
    is colored by ITS correlation strength on a single sequential colormap
    shared with the barplot -- encoding the one thing that's actually being
    compared across subjects (fit quality), rather than uninterpretable
    per-subject categorical hues.
    """
    if scope not in ('exemplars', 'all'):
        raise ValueError(f"scope must be 'exemplars' or 'all', got {scope!r}")

    available = sorted(
        os.path.basename(f).split('_')[0]
        for f in glob.glob(os.path.join(FEM_P_RMS_DIR, f'sub-*_snr{SNR}_fem_p_rms.csv'))
    )
    # OUTLIER_SUBJECT excluded either way: it isn't one of the exemplars, and
    # in 'all' scope its ~20-30x-elevated magnitude would stretch the log-log
    # axes and correlation color normalization for every other subject.
    if scope == 'exemplars':
        subjects = [s for s in EXEMPLAR_SUBJECTS if s in available and s != OUTLIER_SUBJECT]
    else:
        subjects = [s for s in available if s != OUTLIER_SUBJECT]
    per_subject = []
    for subject in subjects:
        prim, sec = load_primary_secondary_pair(subject)
        if prim is None:
            continue
        log_prim, log_sec = np.log10(prim), np.log10(sec)
        r = np.corrcoef(log_prim, log_sec)[0, 1]
        slope, intercept = np.polyfit(log_prim, log_sec, 1)
        per_subject.append(dict(subject=subject, prim=prim, sec=sec, r=r, slope=slope, intercept=intercept))
    if not per_subject:
        print('  no subjects with both primary and secondary dipole P yet -- skipping fig 3b')
        return
    per_subject.sort(key=lambda d: d['r'])

    # Total secondary P (Eq. 6, vsum pipeline) only exists for subjects that
    # have been through eeg_dipole_vsum_sec_curr.py -- absent for the rest,
    # not an error, just a bar that isn't drawn for that subject.
    for d in per_subject:
        p_output = load_subject_vsum_power(d['subject'])
        d['total_sec_p'] = float(p_output['p_els'].sum()) if p_output is not None else None

    cmap = plt.get_cmap('viridis')
    norm = plt.Normalize(vmin=per_subject[0]['r'], vmax=per_subject[-1]['r'])

    fig = plt.figure(figsize=(15, max(6.5, 0.32 * len(per_subject))))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.7], wspace=0.35)
    ax_bar = fig.add_subplot(gs[0])
    ax_scatter = fig.add_subplot(gs[1])

    # Panel A: per-subject magnitude of primary dipole P, total secondary P
    # (Eq. 6, where computed), and secondary dipole P (Eq. 7) -- collapses
    # the 3 magnitude columns of the original Fig. 3B's per-subject table
    # into one grouped bar chart. Correlation itself isn't repeated here --
    # it's already the scatter panel's regression-line color.
    bar_series = [
        ('Primary dipole P', lambda d: d['prim'].sum(), '#4C72B0'),
        ('Total secondary P (Eq. 6)', lambda d: d['total_sec_p'], '#DD8452'),
        ('Secondary dipole P (Eq. 7)', lambda d: d['sec'].sum(), '#55A868'),
    ]
    n_series = len(bar_series)
    bar_height = 0.8 / n_series
    ys = np.arange(len(per_subject))
    for i, (label, value_fn, color) in enumerate(bar_series):
        rows, values = [], []
        for j, d in enumerate(per_subject):
            v = value_fn(d)
            if v is not None:
                rows.append(ys[j] + (i - (n_series - 1) / 2) * bar_height)
                values.append(v)
        if values:
            ax_bar.barh(rows, values, height=bar_height, color=color, label=label)
    ax_bar.set_yticks(ys)
    ax_bar.set_yticklabels([d['subject'] for d in per_subject], fontsize=11)
    ax_bar.set_xscale('log')
    ax_bar.set_xlabel('Power (W)', fontsize=14)
    ax_bar.set_title(f'Per-subject P magnitude (n={len(per_subject)})', fontsize=16)
    ax_bar.legend(fontsize=12, loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=3)
    ax_bar.tick_params(axis='x', labelsize=12)
    for spine in ('top', 'right'):
        ax_bar.spines[spine].set_visible(False)

    # Panel B: subsampled low-alpha scatter (density signal only) + one
    # regression line per subject (identity + fit-quality signal)
    rng = np.random.default_rng(seed)
    all_prim = np.concatenate([d['prim'] for d in per_subject])
    x_range = np.array([all_prim.min(), all_prim.max()])
    for d in per_subject:
        n = len(d['prim'])
        idx = rng.choice(n, size=min(scatter_subsample, n), replace=False)
        ax_scatter.scatter(d['prim'][idx], d['sec'][idx], s=3, color=cmap(norm(d['r'])),
                            alpha=0.3, linewidths=0, rasterized=True)
    for d in per_subject:
        y_fit = 10 ** (d['slope'] * np.log10(x_range) + d['intercept'])
        ax_scatter.plot(x_range, y_fit, color=cmap(norm(d['r'])), lw=1.1, alpha=0.9)

    ax_scatter.set_xscale('log')
    ax_scatter.set_yscale('log')
    ax_scatter.set_xlabel('Primary dipole P (W)', fontsize=14)
    ax_scatter.set_ylabel('Secondary dipole P (W)', fontsize=14)
    scope_label = 'representative subjects' if scope == 'exemplars' else 'all subjects'
    ax_scatter.set_title(f'Primary vs. secondary dipole P, {scope_label} (n={len(per_subject)})', fontsize=16)
    ax_scatter.grid(True, which='major', alpha=0.15, lw=0.5)
    ax_scatter.tick_params(axis='both', labelsize=12)
    for spine in ('top', 'right'):
        ax_scatter.spines[spine].set_visible(False)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax_scatter, fraction=0.05, pad=0.02)
    cbar.set_label('Per-subject Pearson r', fontsize=13)
    cbar.ax.tick_params(labelsize=11)

    fig.suptitle(f'Figure 3B: Primary vs. secondary dipole P ({scope_label})', fontsize=18)
    out_file = os.path.join(FIG_DIR, f'fig3b_primary_vs_secondary_dipole_p_snr{SNR}_{scope}.png')
    fig.savefig(out_file, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {out_file} (n={len(per_subject)} subjects, '
          f'r range [{per_subject[0]["r"]:.3f}, {per_subject[-1]["r"]:.3f}])')


def main():
    os.makedirs(FIG_DIR, exist_ok=True)

    vc_subjects = find_subjects_with_vc()
    vsum_subjects = find_subjects_with_vsum()
    exemplars = [s for s in EXEMPLAR_SUBJECTS if s in vsum_subjects]
    # OUTLIER_SUBJECT stays in the individual/exemplar panels (it's never one
    # of the exemplars anyway) and in Figure 3B's per-subject listing -- only
    # dropped from cross-subject averages/summary stats, never silently
    # absent from a per-subject figure.
    vc_subjects_for_summary = [s for s in vc_subjects if s != OUTLIER_SUBJECT]
    vsum_subjects_for_summary = [s for s in vsum_subjects if s != OUTLIER_SUBJECT]

    print('Volume-by-tissue pie charts (needs only convert_fem_mat_to_vc.py output)...')
    make_average_pie_figure(vc_subjects_for_summary, 'Volume', volume_values, f'fig_volume_by_tissue_average.png')
    make_average_pie_figure(
        [s for s in EXEMPLAR_SUBJECTS if s in vc_subjects and s != OUTLIER_SUBJECT],
        'Volume', volume_values, 'fig_volume_by_tissue_average_exemplars.png')
    make_individual_pie_grid(
        [s for s in EXEMPLAR_SUBJECTS if s in vc_subjects], 'Volume', volume_values,
        'fig_volume_by_tissue_individual.png')

    print('Power-dissipation-by-tissue pie charts (needs eeg_dipole_vsum_sec_curr.py output)...')
    make_average_pie_figure(vsum_subjects_for_summary, 'Power Dissipation', power_values,
                             'fig_power_by_tissue_average.png')
    make_individual_pie_grid(exemplars, 'Power Dissipation', power_values,
                              'fig_power_by_tissue_individual.png')

    print('Single-dipole SEP validation pie charts (needs single_dipole_sep_validation.py output)...')
    make_sep_validation_pie_charts()

    print('Secondary-current-by-tissue pie charts...')
    make_average_pie_figure(vsum_subjects_for_summary, 'Squared Secondary Current', current_values,
                             'fig_current_by_tissue_average.png')
    make_individual_pie_grid(exemplars, 'Squared Secondary Current', current_values,
                              'fig_current_by_tissue_individual.png')

    print('Secondary dipole P histogram (needs eeg_dipole_power_computation.py output)...')
    make_secondary_dipole_p_histogram()

    print('Figure 3B: primary vs. secondary dipole P (needs primary_p and eeg_dipole_power_computation.py output)...')
    make_figure_3b(scope='exemplars')
    make_figure_3b(scope='all')

    print('Total secondary P as % of total primary P (needs eeg_dipole_vsum_sec_curr.py output)...')
    print_total_secondary_p_ratio(scope='exemplars')
    print_total_secondary_p_ratio(scope='all')

    print('Exemplar validation stats (brain tissue %, average secondary power)...')
    print_exemplar_validation_stats()


if __name__ == '__main__':
    main()
