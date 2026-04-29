#!/usr/bin/env python3
"""
generate_difficulty_supps.py — Red Analysis Pipeline

Difficulty-stratified supplementary figures.

  figS6: AMBER pre-Rosetta dual-effect by BM5.5 difficulty tier (3 panels)
  figS7: Rosetta protocol MolProbity by difficulty tier (3 panels)
  figS8: TM/MP tradeoff by difficulty tier (3 panels), highlighting dualspace_beta
         and cartesian_beta cells

Targets per tier: rigid n=162, medium n=60, difficult n=35.
"""

import os
import sqlite3
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

DB = "/data/p_csb_meiler/agarwm5/red_analysis/db/bm55.sqlite"
METRICS_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"
FIGDIR = "/data/p_csb_meiler/agarwm5/red_analysis/figures"

DIFF_ORDER = ['rigid', 'medium', 'difficult']
DIFF_COLORS = {
    'rigid':     '#56B4E9',
    'medium':    '#E69F00',
    'difficult': '#D55E00',
}
DIFF_N = {'rigid': 162, 'medium': 60, 'difficult': 35}

PROTOCOL_ORDER = ['cartesian_beta', 'cartesian_ref15',
                  'dualspace_beta', 'dualspace_ref15',
                  'normal_beta',    'normal_ref15']
PROTO_COLORS = {
    'cartesian_beta':  '#0072B2',
    'cartesian_ref15': '#56B4E9',
    'dualspace_beta':  '#009E73',
    'dualspace_ref15': '#CC79A7',
    'normal_beta':     '#E69F00',
    'normal_ref15':    '#D55E00',
}
PROTO_LABELS = {
    'cartesian_beta':  'Cart-beta',
    'cartesian_ref15': 'Cart-REF15',
    'dualspace_beta':  'Dual-beta',
    'dualspace_ref15': 'Dual-REF15',
    'normal_beta':     'Norm-beta',
    'normal_ref15':    'Norm-REF15',
}

SOURCE_ORDER = ['crystal', 'af_relaxed', 'af_unrelaxed', 'boltz',
                'amber_af', 'amber_boltz']
SOURCE_LABELS = {
    'crystal':      'Crystal',
    'af_relaxed':   'AF rlx',
    'af_unrelaxed': 'AF unrlx',
    'boltz':        'Boltz',
    'amber_af':     'AMBER(AF)',
    'amber_boltz':  'AMBER(Boltz)',
}


def load_targets():
    con = sqlite3.connect(DB)
    df = pd.read_sql("SELECT target_id, difficulty FROM targets", con)
    con.close()
    return df


def classify_src(s):
    s = s.lower()
    if s.startswith('crystal'): return 'crystal'
    if s.startswith('amber_af'): return 'amber_af'
    if s.startswith('amber_boltz'): return 'amber_boltz'
    if s.startswith('amber_crystal'): return 'amber_crystal'
    if s.startswith('af_relaxed'): return 'af_relaxed'
    if s.startswith('af_unrelaxed'): return 'af_unrelaxed'
    if s.startswith('boltz'): return 'boltz'
    return 'unknown'


def load_pre_mp(targets):
    """Pre-Rosetta MolProbity per target/source/pipeline."""
    f = os.path.join(METRICS_DIR, 'combined_molprobity.tsv')
    df = pd.read_csv(f, sep='\t')
    for c in ['molprobity_score', 'clashscore']:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    df['source'] = df['source'].apply(classify_src)
    by = df.groupby(['target', 'pipeline', 'source'], as_index=False).agg(
        mp=('molprobity_score', 'mean'),
        clash=('clashscore', 'mean'))
    by = by.merge(targets, left_on='target', right_on='target_id', how='left')
    return by


def load_pre_tm(targets):
    """Pre-Rosetta TM per target/source/pipeline."""
    f = os.path.join(METRICS_DIR, 'combined_tmscore.tsv')
    df = pd.read_csv(f, sep='\t')
    df['tmscore'] = pd.to_numeric(df['tmscore'], errors='coerce')
    by = df.groupby(['target', 'pipeline', 'source'], as_index=False).agg(
        tm=('tmscore', 'mean'))
    by = by.merge(targets, left_on='target', right_on='target_id', how='left')
    return by


def load_rosetta_mp_tm(targets):
    """Rosetta MP and TM per target/source/protocol/pipeline."""
    mp_f = os.path.join(METRICS_DIR, 'combined_rosetta_molprobity.tsv')
    tm_f = os.path.join(METRICS_DIR, 'combined_rosetta_tmscore.tsv')
    mp = pd.read_csv(mp_f, sep='\t')
    tm = pd.read_csv(tm_f, sep='\t')
    for c in ['molprobity_score', 'clashscore']:
        mp[c] = pd.to_numeric(mp[c], errors='coerce')
    tm['tmscore'] = pd.to_numeric(tm['tmscore'], errors='coerce')
    mp['source'] = mp['source'].apply(classify_src)
    tm['source'] = tm['source'].apply(classify_src)

    mp_by = mp.groupby(['target', 'pipeline', 'source', 'protocol'],
                        as_index=False).agg(mp=('molprobity_score', 'mean'),
                                            clash=('clashscore', 'mean'))
    tm_by = tm.groupby(['target', 'pipeline', 'source', 'protocol'],
                        as_index=False).agg(tm=('tmscore', 'mean'))
    mp_by = mp_by.merge(targets, left_on='target', right_on='target_id', how='left')
    tm_by = tm_by.merge(targets, left_on='target', right_on='target_id', how='left')
    return mp_by, tm_by


def save_fig(fig, name):
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(FIGDIR, f"{name}.{fmt}"),
                    dpi=300 if fmt == 'png' else None,
                    bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {name}")


# ----------------- figS6 -----------------------

def figS6_amber_dual_effect_by_difficulty(pre_mp, pre_tm):
    """Three-panel AMBER pre-Rosetta dual-effect, one per difficulty tier.

    Each panel: per-target Delta TM (x) vs Delta MolProbity score (y), where the
    pair is (af_unrelaxed -> amber_af) or (boltz -> amber_boltz). Pooled across
    pipelines (blue + green) for visual density. The d=-0.99 effect lives in the
    pre-Rosetta MP delta direction.
    """
    pairs = [('af_unrelaxed', 'amber_af',    '#1f6feb', 'AF unrlx -> AMBER(AF)'),
             ('boltz',        'amber_boltz', '#2ea043', 'Boltz -> AMBER(Boltz)')]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), sharey=True)

    for ax, diff in zip(axes, DIFF_ORDER):
        n_targets = DIFF_N[diff]
        for src_b, src_a, col, lbl in pairs:
            sub_b_mp = pre_mp[(pre_mp.source == src_b) & (pre_mp.difficulty == diff)]
            sub_a_mp = pre_mp[(pre_mp.source == src_a) & (pre_mp.difficulty == diff)]
            sub_b_tm = pre_tm[(pre_tm.source == src_b) & (pre_tm.difficulty == diff)]
            sub_a_tm = pre_tm[(pre_tm.source == src_a) & (pre_tm.difficulty == diff)]

            # Per (target, pipeline) pairing
            b_mp = sub_b_mp.groupby(['target', 'pipeline'])['mp'].mean().rename('b_mp')
            a_mp = sub_a_mp.groupby(['target', 'pipeline'])['mp'].mean().rename('a_mp')
            b_tm = sub_b_tm.groupby(['target', 'pipeline'])['tm'].mean().rename('b_tm')
            a_tm = sub_a_tm.groupby(['target', 'pipeline'])['tm'].mean().rename('a_tm')

            joined = pd.concat([b_mp, a_mp, b_tm, a_tm], axis=1).dropna()
            if len(joined) == 0:
                continue
            d_tm = joined['a_tm'] - joined['b_tm']
            d_mp = joined['a_mp'] - joined['b_mp']

            ax.scatter(d_tm, d_mp, alpha=0.55, s=22, color=col,
                       edgecolors='black', linewidths=0.25, label=lbl)

            # Overlay group mean
            ax.scatter(d_tm.mean(), d_mp.mean(), marker='*', s=240,
                       color=col, edgecolors='black', linewidths=1.0, zorder=10)

        ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.6, alpha=0.6)
        ax.axvline(x=0, color='gray', linestyle='-', linewidth=0.6, alpha=0.6)
        ax.set_xlabel('Delta TM-score (after - before AMBER)', fontsize=10)
        if ax is axes[0]:
            ax.set_ylabel('Delta MolProbity score (after - before AMBER)\n'
                          'lower = better', fontsize=10)
        ax.set_title(f'{diff.capitalize()}  (n={n_targets} targets)',
                     fontsize=11, fontweight='bold')
        ax.set_xlim(-0.04, 0.04)
        ax.legend(fontsize=8, loc='upper left')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle('AMBER pre-Rosetta dual-effect by BM5.5 difficulty tier',
                 fontsize=13, fontweight='bold')
    fig.text(0.5, 0.005,
             'TM-score is unchanged; MolProbity collapses by ~0.5-1.0 across all '
             'tiers. Stars mark per-pair group means.',
             fontsize=8.5, ha='center', style='italic', color='#333')
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    save_fig(fig, 'figS6_amber_dual_effect_by_difficulty')


# ----------------- figS7 -----------------------

def figS7_rosetta_protocol_mp_by_difficulty(ros_mp):
    """Three-panel Rosetta protocol MP score by difficulty tier.

    Each panel: 6 protocols x 6 sources, box plots, MP score y-axis. Difficulty
    fixed per panel. Visualizes whether protocol ranking is stable across tiers.
    """
    fig, axes = plt.subplots(1, 3, figsize=(20, 6), sharey=True)

    pooled = ros_mp.groupby(['target', 'difficulty', 'source', 'protocol'],
                            as_index=False).agg(mp=('mp', 'mean'))

    for ax, diff in zip(axes, DIFF_ORDER):
        sub = pooled[pooled.difficulty == diff]
        # Build box-plot data: x positions = source x protocol (grouped)
        positions = []
        widths = []
        plot_data = []
        plot_colors = []
        labels = []
        x = 0
        group_centers = []
        for s_idx, src in enumerate(SOURCE_ORDER):
            group_start = x
            for p_idx, proto in enumerate(PROTOCOL_ORDER):
                vals = sub[(sub.source == src) &
                           (sub.protocol == proto)]['mp'].dropna().values
                if len(vals) >= 2:
                    plot_data.append(vals)
                    positions.append(x)
                    widths.append(0.7)
                    plot_colors.append(PROTO_COLORS[proto])
                    labels.append(PROTO_LABELS[proto])
                x += 1
            x += 1
            group_centers.append((group_start + x - 2) / 2.0)

        if plot_data:
            bp = ax.boxplot(plot_data, positions=positions, widths=widths,
                            patch_artist=True, showfliers=False,
                            medianprops={'color': 'black', 'linewidth': 1.0})
            for patch, col in zip(bp['boxes'], plot_colors):
                patch.set_facecolor(col)
                patch.set_alpha(0.75)

        # Group separator labels (sources at the bottom)
        ax.set_xticks(group_centers)
        ax.set_xticklabels([SOURCE_LABELS[s] for s in SOURCE_ORDER],
                           rotation=30, fontsize=8)
        ax.set_title(f'{diff.capitalize()}  (n={DIFF_N[diff]} targets)',
                     fontsize=11, fontweight='bold')
        if ax is axes[0]:
            ax.set_ylabel('MolProbity score (lower = better)', fontsize=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Legend (single, top-right)
    from matplotlib.patches import Patch
    legend_handles = [Patch(facecolor=PROTO_COLORS[p], alpha=0.75,
                            label=PROTO_LABELS[p])
                      for p in PROTOCOL_ORDER]
    axes[2].legend(handles=legend_handles, loc='upper right',
                   fontsize=8, ncol=2, frameon=True, title='Protocol',
                   title_fontsize=8)

    fig.suptitle('Rosetta protocol MolProbity by BM5.5 difficulty tier',
                 fontsize=13, fontweight='bold')
    plt.tight_layout(rect=[0, 0.0, 1, 0.96])
    save_fig(fig, 'figS7_rosetta_protocol_mp_by_difficulty')


# ----------------- figS8 -----------------------

def figS8_tm_mp_tradeoff_by_difficulty(pre_mp, pre_tm, ros_mp, ros_tm):
    """Three-panel TM/MP tradeoff by difficulty tier.

    Per panel: per-target Delta TM (x) vs Delta MolProbity (y) where the delta is
    (post-Rosetta - pre-Rosetta) for the pre-Rosetta source af_relaxed pooled
    across both pipelines. Highlight dualspace_beta and cartesian_beta cells.
    """
    headline_protos = ['dualspace_beta', 'cartesian_beta']
    other_protos = [p for p in PROTOCOL_ORDER if p not in headline_protos]
    src = 'af_relaxed'  # main source for the tradeoff

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), sharey=True)

    # Pre values per (target, pipeline)
    pre_mp_idx = pre_mp[pre_mp.source == src].groupby(
        ['target', 'pipeline', 'difficulty'], as_index=False).agg(
        pre_mp=('mp', 'mean'))
    pre_tm_idx = pre_tm[pre_tm.source == src].groupby(
        ['target', 'pipeline', 'difficulty'], as_index=False).agg(
        pre_tm=('tm', 'mean'))
    pre = pre_mp_idx.merge(pre_tm_idx,
                           on=['target', 'pipeline', 'difficulty'])

    ros_mp_src = ros_mp[ros_mp.source == src].groupby(
        ['target', 'pipeline', 'difficulty', 'protocol'], as_index=False).agg(
        post_mp=('mp', 'mean'))
    ros_tm_src = ros_tm[ros_tm.source == src].groupby(
        ['target', 'pipeline', 'difficulty', 'protocol'], as_index=False).agg(
        post_tm=('tm', 'mean'))

    joined = pre.merge(ros_mp_src, on=['target', 'pipeline', 'difficulty']) \
                .merge(ros_tm_src,
                       on=['target', 'pipeline', 'difficulty', 'protocol'])
    joined['d_tm'] = joined['post_tm'] - joined['pre_tm']
    joined['d_mp'] = joined['post_mp'] - joined['pre_mp']

    for ax, diff in zip(axes, DIFF_ORDER):
        sub = joined[joined.difficulty == diff]
        # Plot non-headline protos as small grey dots first
        oth = sub[sub.protocol.isin(other_protos)]
        ax.scatter(oth['d_tm'], oth['d_mp'],
                   color='#bbbbbb', alpha=0.30, s=10, edgecolors='none',
                   label='Other protocols')

        for proto in headline_protos:
            psub = sub[sub.protocol == proto]
            if len(psub) == 0:
                continue
            ax.scatter(psub['d_tm'], psub['d_mp'],
                       color=PROTO_COLORS[proto], alpha=0.75, s=30,
                       edgecolors='black', linewidths=0.3,
                       label=PROTO_LABELS[proto])
            ax.scatter(psub['d_tm'].mean(), psub['d_mp'].mean(),
                       marker='*', s=260, color=PROTO_COLORS[proto],
                       edgecolors='black', linewidths=1.0, zorder=10)

        ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.6, alpha=0.6)
        ax.axvline(x=0, color='gray', linestyle='-', linewidth=0.6, alpha=0.6)
        ax.set_xlabel('Delta TM-score (Rosetta - pre)', fontsize=10)
        if ax is axes[0]:
            ax.set_ylabel('Delta MolProbity (Rosetta - pre)\nlower = MP improved',
                          fontsize=10)
        ax.set_title(f'{diff.capitalize()}  (n={DIFF_N[diff]} targets)',
                     fontsize=11, fontweight='bold')
        ax.legend(fontsize=7.5, loc='best')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle('Rosetta TM/MP tradeoff by BM5.5 difficulty tier '
                 '(source = af_relaxed)', fontsize=12.5, fontweight='bold')
    fig.text(0.5, 0.005,
             'Headline cells: dualspace_beta and cartesian_beta. Stars show '
             'per-protocol means. Bottom-left quadrant = ideal (better MP, '
             'preserved TM).',
             fontsize=8.5, ha='center', style='italic', color='#333')
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    save_fig(fig, 'figS8_tm_mp_tradeoff_by_difficulty')


# ----------------- figS9 (optional outlier) -----------------------

def figS9_boltz_outlier_1Y64(pre_mp):
    """Per-target Delta MP score (boltz -> amber_boltz) with 1Y64 highlighted."""
    boltz = pre_mp[pre_mp.source == 'boltz'].groupby(
        ['target', 'pipeline'], as_index=False).agg(b=('mp', 'mean'))
    amber = pre_mp[pre_mp.source == 'amber_boltz'].groupby(
        ['target', 'pipeline'], as_index=False).agg(a=('mp', 'mean'))
    d = boltz.merge(amber, on=['target', 'pipeline'])
    d['delta'] = d['a'] - d['b']
    # Use blue pipeline for the headline; green is the second bar
    d_blue = d[d.pipeline == 'blue'].sort_values('delta').reset_index(drop=True)
    if len(d_blue) == 0:
        print("  SKIP figS9: no blue pipeline data")
        return

    fig, ax = plt.subplots(figsize=(13, 5.5))
    cols = ['#888888'] * len(d_blue)
    is_1y64 = d_blue['target'] == '1Y64'
    if is_1y64.any():
        i_1y64 = int(np.where(is_1y64)[0][0])
        cols[i_1y64] = '#D55E00'

    ax.bar(np.arange(len(d_blue)), d_blue['delta'].values,
           color=cols, edgecolor='black', linewidth=0.2, width=0.9)
    ax.axhline(y=0, color='black', linewidth=0.6)

    # Highlight 1Y64
    if is_1y64.any():
        ax.annotate(f'1Y64\nDelta MP = {d_blue.loc[i_1y64, "delta"]:+.2f}',
                    xy=(i_1y64, d_blue.loc[i_1y64, 'delta']),
                    xytext=(i_1y64 + 5, d_blue.loc[i_1y64, 'delta']),
                    fontsize=10, color='#D55E00',
                    arrowprops=dict(arrowstyle='->',
                                    color='#D55E00', lw=1.0))

    ax.set_xlabel(f'Target index (sorted by Delta MP, n={len(d_blue)})',
                  fontsize=10)
    ax.set_ylabel('Delta MolProbity score (amber_boltz - boltz)',
                  fontsize=10)
    ax.set_title('Boltz outlier: 1Y64 is the one target where AMBER worsens MP '
                 '(Blue pipeline)', fontsize=12, fontweight='bold')
    # Reference text
    median_d = float(np.median(d_blue['delta']))
    n_worse = int((d_blue['delta'] > 0).sum())
    ax.text(0.02, 0.97,
            f'Distribution stats:\n'
            f'  median Delta MP = {median_d:.3f}\n'
            f'  n worsened (Delta>0) = {n_worse} of {len(d_blue)}',
            transform=ax.transAxes, fontsize=8.5, va='top',
            bbox=dict(boxstyle='round,pad=0.4',
                      facecolor='#f7f7f7', edgecolor='#aaa', linewidth=0.6))

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    save_fig(fig, 'figS9_boltz_outlier_1Y64')


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    targets = load_targets()
    print(f"Loaded {len(targets)} targets, "
          f"{(targets.difficulty == 'rigid').sum()} rigid, "
          f"{(targets.difficulty == 'medium').sum()} medium, "
          f"{(targets.difficulty == 'difficult').sum()} difficult")

    pre_mp = load_pre_mp(targets)
    pre_tm = load_pre_tm(targets)
    ros_mp, ros_tm = load_rosetta_mp_tm(targets)

    print(f"  pre_mp rows={len(pre_mp)}  pre_tm rows={len(pre_tm)}")
    print(f"  ros_mp rows={len(ros_mp)}  ros_tm rows={len(ros_tm)}")

    print("Generating figS6...")
    figS6_amber_dual_effect_by_difficulty(pre_mp, pre_tm)
    print("Generating figS7...")
    figS7_rosetta_protocol_mp_by_difficulty(ros_mp)
    print("Generating figS8...")
    figS8_tm_mp_tradeoff_by_difficulty(pre_mp, pre_tm, ros_mp, ros_tm)
    print("Generating figS9...")
    figS9_boltz_outlier_1Y64(pre_mp)
    print("Done.")


if __name__ == '__main__':
    main()
