#!/usr/bin/env python3
"""
generate_rosetta_figures.py — Red Analysis Pipeline

Publication figures for Rosetta relaxation analysis.

Figures:
  fig7: Protocol comparison (TM-score by protocol, grouped bar)
  fig8: Rosetta effect (pre vs post TM-score, paired scatter)
  fig9: Protocol × Source heatmap
  fig10: Replicate convergence (std dev histogram)
  fig11: Blue-Green Rosetta agreement scatter
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

METRICS_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"
FIGDIR = "/data/p_csb_meiler/agarwm5/red_analysis/figures"

# Colorblind-friendly palette (Wong 2011)
COLORS = {
    'cartesian_beta': '#E69F00',   # Orange
    'cartesian_ref15': '#56B4E9',  # Sky blue
    'dualspace_beta': '#009E73',   # Bluish green
    'dualspace_ref15': '#F0E442',  # Yellow
    'normal_beta': '#0072B2',      # Blue
    'normal_ref15': '#D55E00',     # Vermillion
}

SOURCE_ORDER = ['af_relaxed', 'af_unrelaxed', 'boltz', 'amber_af', 'amber_boltz', 'crystal']
SOURCE_LABELS = {
    'af_relaxed': 'AF2\n(relaxed)',
    'af_unrelaxed': 'AF2\n(unrelaxed)',
    'boltz': 'Boltz-1',
    'amber_af': 'AMBER\n(AF)',
    'amber_boltz': 'AMBER\n(Boltz)',
    'crystal': 'Crystal'
}

PROTOCOL_ORDER = ['cartesian_beta', 'cartesian_ref15', 'dualspace_beta',
                  'dualspace_ref15', 'normal_beta', 'normal_ref15']
PROTOCOL_LABELS = {
    'cartesian_beta': 'Cart β',
    'cartesian_ref15': 'Cart REF15',
    'dualspace_beta': 'Dual β',
    'dualspace_ref15': 'Dual REF15',
    'normal_beta': 'Norm β',
    'normal_ref15': 'Norm REF15'
}


def load_data():
    """Load combined Rosetta TM-score data."""
    path = os.path.join(METRICS_DIR, "combined_rosetta_tmscore.tsv")
    data = pd.read_csv(path, sep='\t')
    for col in ['rmsd', 'tmscore', 'gdtts', 'gdtha']:
        data[col] = pd.to_numeric(data[col], errors='coerce')
    return data


def fig7_protocol_comparison(data):
    """Bar chart of TM-score by protocol for each pipeline."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    for ax, pipe in zip(axes, ['blue', 'green']):
        pipe_data = data[data['pipeline'] == pipe]
        # Per-target mean across reps and source types
        per_target = pipe_data.groupby(['target', 'protocol'])['tmscore'].mean().reset_index()
        by_proto = per_target.groupby('protocol')['tmscore'].agg(['mean', 'std']).reset_index()
        by_proto = by_proto.set_index('protocol').loc[PROTOCOL_ORDER].reset_index()

        bars = ax.bar(range(len(PROTOCOL_ORDER)),
                     by_proto['mean'],
                     yerr=by_proto['std'],
                     color=[COLORS[p] for p in PROTOCOL_ORDER],
                     edgecolor='black', linewidth=0.5,
                     capsize=3, error_kw={'linewidth': 1})

        ax.set_xticks(range(len(PROTOCOL_ORDER)))
        ax.set_xticklabels([PROTOCOL_LABELS[p] for p in PROTOCOL_ORDER], fontsize=9, rotation=30, ha='right')
        ax.set_ylabel('TM-score' if pipe == 'blue' else '')
        ax.set_title(f'{pipe.capitalize()} Pipeline')
        ax.set_ylim(0.88, 0.96)
        ax.axhline(y=0.945, color='gray', linestyle='--', linewidth=0.5, alpha=0.5, label='Pre-Rosetta avg')
        if pipe == 'blue':
            ax.legend(fontsize=8, loc='upper right')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle('Rosetta Protocol Comparison (TM-score to Crystal)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    save_fig(fig, 'fig7_rosetta_protocol_comparison')


def fig8_rosetta_effect(data):
    """Paired scatter: pre-Rosetta vs post-Rosetta TM-score."""
    pre_path = os.path.join(METRICS_DIR, "combined_tmscore.tsv")
    if not os.path.exists(pre_path):
        print("  SKIP fig8: no pre-Rosetta data")
        return

    pre_data = pd.read_csv(pre_path, sep='\t')
    pre_data['tmscore'] = pd.to_numeric(pre_data['tmscore'], errors='coerce')

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    sources = ['af_relaxed', 'af_unrelaxed', 'boltz']

    for ax, src in zip(axes, sources):
        for pipe, marker, color in [('blue', 'o', '#0072B2'), ('green', 's', '#009E73')]:
            pre = pre_data[(pre_data['pipeline'] == pipe) & (pre_data['source'] == src)]
            pre_by_target = pre.groupby('target')['tmscore'].mean()

            post = data[(data['pipeline'] == pipe) & (data['source'] == src)]
            post_by_target = post.groupby('target')['tmscore'].mean()

            common = pre_by_target.index.intersection(post_by_target.index)
            if len(common) < 5:
                continue

            x = pre_by_target.loc[common].values
            y = post_by_target.loc[common].values
            mask = np.isfinite(x) & np.isfinite(y)

            ax.scatter(x[mask], y[mask], alpha=0.5, s=20, marker=marker,
                      color=color, label=f'{pipe.capitalize()}')

        ax.plot([0.5, 1.0], [0.5, 1.0], 'k--', linewidth=0.5, alpha=0.3)
        ax.set_xlabel('Pre-Rosetta TM-score')
        ax.set_ylabel('Post-Rosetta TM-score')
        ax.set_title(SOURCE_LABELS.get(src, src).replace('\n', ' '))
        ax.set_xlim(0.5, 1.0)
        ax.set_ylim(0.5, 1.0)
        ax.set_aspect('equal')
        ax.legend(fontsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle('Rosetta Effect: Pre vs Post TM-score\n(below diagonal = degradation)',
                fontsize=13, fontweight='bold')
    plt.tight_layout()
    save_fig(fig, 'fig8_rosetta_effect')


def fig9_protocol_source_heatmap(data):
    """Heatmap of TM-score by protocol × source."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, pipe in zip(axes, ['blue', 'green']):
        pipe_data = data[data['pipeline'] == pipe]
        per_target = pipe_data.groupby(['target', 'source', 'protocol'])['tmscore'].mean().reset_index()
        pivot = per_target.groupby(['source', 'protocol'])['tmscore'].mean().unstack()

        # Reorder
        valid_sources = [s for s in SOURCE_ORDER if s in pivot.index]
        valid_protos = [p for p in PROTOCOL_ORDER if p in pivot.columns]
        if not valid_sources or not valid_protos:
            continue
        pivot = pivot.loc[valid_sources, valid_protos]
        pivot.index = [SOURCE_LABELS.get(s, s).replace('\n', ' ') for s in valid_sources]
        pivot.columns = [PROTOCOL_LABELS[p] for p in valid_protos]

        sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdYlGn', ax=ax,
                   vmin=0.85, vmax=0.96, linewidths=0.5,
                   cbar_kws={'label': 'TM-score'})
        ax.set_title(f'{pipe.capitalize()} Pipeline')

    fig.suptitle('TM-score by Input Type × Rosetta Protocol', fontsize=14, fontweight='bold')
    plt.tight_layout()
    save_fig(fig, 'fig9_protocol_source_heatmap')


def fig10_replicate_convergence(data):
    """Histogram of replicate std dev."""
    per_rep = data.groupby(['pipeline', 'source', 'target', 'protocol', 'rep'])[
        'tmscore'].mean().reset_index()
    rep_std = per_rep.groupby(['pipeline', 'source', 'target', 'protocol'])[
        'tmscore'].std().reset_index()
    rep_std.columns = ['pipeline', 'source', 'target', 'protocol', 'rep_std']
    rep_std = rep_std.dropna(subset=['rep_std'])

    fig, ax = plt.subplots(figsize=(8, 5))

    for pipe, color in [('blue', '#0072B2'), ('green', '#009E73')]:
        vals = rep_std[rep_std['pipeline'] == pipe]['rep_std'].values
        ax.hist(vals, bins=50, alpha=0.6, color=color, label=f'{pipe.capitalize()} (n={len(vals)})',
               edgecolor='white', linewidth=0.3)

    ax.axvline(x=0.01, color='red', linestyle='--', linewidth=1, label='Convergence threshold (0.01)')
    ax.set_xlabel('Std Dev of TM-score across 5 replicates')
    ax.set_ylabel('Count')
    ax.set_title('Rosetta Replicate Convergence', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    save_fig(fig, 'fig10_replicate_convergence')


def fig11_blue_green_rosetta(data):
    """Scatter: Blue vs Green Rosetta TM-score."""
    from scipy import stats as sp_stats

    per_target = data.groupby(['pipeline', 'target', 'protocol'])['tmscore'].mean().reset_index()

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for idx, proto in enumerate(PROTOCOL_ORDER):
        ax = axes[idx]
        blue = per_target[(per_target['pipeline'] == 'blue') & (per_target['protocol'] == proto)]
        green = per_target[(per_target['pipeline'] == 'green') & (per_target['protocol'] == proto)]

        b_by_t = blue.set_index('target')['tmscore']
        g_by_t = green.set_index('target')['tmscore']
        common = b_by_t.index.intersection(g_by_t.index)

        if len(common) > 5:
            bv = b_by_t.loc[common].values
            gv = g_by_t.loc[common].values
            mask = np.isfinite(bv) & np.isfinite(gv)

            ax.scatter(bv[mask], gv[mask], alpha=0.5, s=15,
                      color=COLORS[proto], edgecolors='black', linewidths=0.3)
            ax.plot([0.5, 1.0], [0.5, 1.0], 'k--', linewidth=0.5, alpha=0.3)

            r, p = sp_stats.pearsonr(bv[mask], gv[mask])
            ax.text(0.05, 0.95, f'r = {r:.3f}', transform=ax.transAxes,
                   fontsize=10, verticalalignment='top')

        ax.set_xlabel('Blue TM-score')
        ax.set_ylabel('Green TM-score')
        ax.set_title(PROTOCOL_LABELS[proto])
        ax.set_xlim(0.5, 1.0)
        ax.set_ylim(0.5, 1.0)
        ax.set_aspect('equal')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle('Blue vs Green Rosetta TM-score Agreement\nby Protocol',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    save_fig(fig, 'fig11_blue_green_rosetta')


def save_fig(fig, name):
    """Save figure as PDF and PNG."""
    for fmt in ['pdf', 'png']:
        path = os.path.join(FIGDIR, f"{name}.{fmt}")
        fig.savefig(path, dpi=300 if fmt == 'png' else None, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {name}")


def main():
    print("Loading Rosetta TM-score data...")
    data = load_data()
    print(f"  {len(data)} rows, {data['target'].nunique()} targets")

    os.makedirs(FIGDIR, exist_ok=True)

    print("\nGenerating figures...")
    fig7_protocol_comparison(data)
    fig8_rosetta_effect(data)
    fig9_protocol_source_heatmap(data)
    fig10_replicate_convergence(data)
    fig11_blue_green_rosetta(data)

    print("\nAll Rosetta figures generated!")


if __name__ == '__main__':
    main()
