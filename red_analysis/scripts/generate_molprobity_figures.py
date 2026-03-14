#!/usr/bin/env python3
"""
generate_molprobity_figures.py — Red Analysis Pipeline

Publication figures for MolProbity validation analysis.

Figures:
  fig12: MolProbity metrics by source (grouped bars)
  fig13: AMBER effect on clashscore (before/after paired plot)
  fig14: TM-score vs MolProbity scatter (the dual-effect plot)
  fig15: Crystal vs predicted MolProbity comparison
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

METRICS_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"
FIGDIR = "/data/p_csb_meiler/agarwm5/red_analysis/figures"

# Source colors (same palette as TM-score figures)
SRC_COLORS = {
    'crystal': '#999999',
    'af_relaxed': '#0072B2',
    'af_unrelaxed': '#56B4E9',
    'boltz': '#009E73',
    'amber_af': '#E69F00',
    'amber_boltz': '#D55E00',
}

SOURCE_ORDER = ['crystal', 'af_relaxed', 'af_unrelaxed', 'boltz', 'amber_af', 'amber_boltz']
SOURCE_LABELS = {
    'crystal': 'Crystal',
    'af_relaxed': 'AF2\n(relaxed)',
    'af_unrelaxed': 'AF2\n(unrelaxed)',
    'boltz': 'Boltz-1',
    'amber_af': 'AMBER\n(AF)',
    'amber_boltz': 'AMBER\n(Boltz)',
}


def load_data():
    mp = pd.read_csv(os.path.join(METRICS_DIR, "combined_molprobity.tsv"), sep='\t')
    for col in ['clashscore', 'rama_outliers', 'rama_favored', 'rota_outliers',
                'molprobity_score', 'cbeta_outliers', 'rms_bonds', 'rms_angles']:
        mp[col] = pd.to_numeric(mp[col], errors='coerce')
    return mp


def save_fig(fig, name):
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(FIGDIR, f"{name}.{fmt}"),
                   dpi=300 if fmt == 'png' else None, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {name}")


def fig12_molprobity_by_source(data):
    """Grouped bar charts for key MolProbity metrics."""
    metrics = [('clashscore', 'Clashscore\n(lower = better)'),
               ('molprobity_score', 'MolProbity Score\n(lower = better)'),
               ('rama_favored', 'Ramachandran Favored %\n(higher = better)'),
               ('rota_outliers', 'Rotamer Outliers %\n(lower = better)')]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    # Per-target mean (average over 5 models)
    per_target = data[data['pipeline'] == 'blue'].groupby(['source', 'target'])[
        [m[0] for m in metrics]].mean().reset_index()

    for ax, (metric, ylabel) in zip(axes, metrics):
        vals = []
        labels = []
        colors = []
        for src in SOURCE_ORDER:
            src_data = per_target[per_target['source'] == src][metric].dropna()
            if len(src_data) > 0:
                vals.append(src_data.values)
                labels.append(SOURCE_LABELS[src])
                colors.append(SRC_COLORS[src])

        bp = ax.boxplot(vals, tick_labels=labels, patch_artist=True, widths=0.6,
                       medianprops={'color': 'black', 'linewidth': 1.5})
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax.set_ylabel(ylabel, fontsize=10)
        ax.tick_params(axis='x', labelsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle('MolProbity Validation Metrics by Structure Source\n(Blue pipeline, per-target averages)',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    save_fig(fig, 'fig12_molprobity_by_source')


def fig13_amber_clashscore_effect(data):
    """Paired before/after plot showing AMBER's effect on clashscore."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    pairs = [
        ('af_unrelaxed', 'amber_af', 'AF unrelaxed → AMBER(AF)', axes[0]),
        ('boltz', 'amber_boltz', 'Boltz-1 → AMBER(Boltz)', axes[1]),
    ]

    for src_before, src_after, title, ax in pairs:
        for pipe, alpha in [('blue', 0.7), ('green', 0.3)]:
            before = data[(data['pipeline'] == pipe) & (data['source'] == src_before)]
            after = data[(data['pipeline'] == pipe) & (data['source'] == src_after)]

            b_mean = before.groupby('target')['clashscore'].mean()
            a_mean = after.groupby('target')['clashscore'].mean()
            common = b_mean.index.intersection(a_mean.index)

            if len(common) < 3:
                continue

            bv = b_mean.loc[common].values
            av = a_mean.loc[common].values

            # Draw lines connecting before and after
            for b, a in zip(bv, av):
                color = 'green' if a < b else 'red'
                ax.plot([0, 1], [b, a], color=color, alpha=alpha*0.3, linewidth=0.5)

            if pipe == 'blue':
                ax.scatter(np.zeros(len(bv)), bv, color=SRC_COLORS[src_before],
                          s=20, alpha=0.6, zorder=5, label='Before')
                ax.scatter(np.ones(len(av)), av, color=SRC_COLORS[src_after],
                          s=20, alpha=0.6, zorder=5, label='After')

        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Before AMBER', 'After AMBER'], fontsize=11)
        ax.set_ylabel('Clashscore')
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.legend(fontsize=9)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle('AMBER Relaxation Effect on Clashscore', fontsize=14, fontweight='bold')
    plt.tight_layout()
    save_fig(fig, 'fig13_amber_clashscore_effect')


def fig14_tm_vs_molprobity(data):
    """THE KEY FIGURE: TM-score (no change) vs MolProbity (big change)."""
    tm_path = os.path.join(METRICS_DIR, "combined_tmscore.tsv")
    if not os.path.exists(tm_path):
        print("  SKIP fig14: no TM-score data")
        return

    tm_data = pd.read_csv(tm_path, sep='\t')
    tm_data['tmscore'] = pd.to_numeric(tm_data['tmscore'], errors='coerce')

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    pipe = 'blue'

    # Panel 1: TM-score change (before vs after AMBER)
    ax = axes[0]
    for src_before, src_after, color, label in [
        ('af_unrelaxed', 'amber_af', '#E69F00', 'AF→AMBER(AF)'),
        ('boltz', 'amber_boltz', '#D55E00', 'Boltz→AMBER(Boltz)'),
    ]:
        b = tm_data[(tm_data['pipeline'] == pipe) & (tm_data['source'] == src_before)]
        a = tm_data[(tm_data['pipeline'] == pipe) & (tm_data['source'] == src_after)]
        bm = b.groupby('target')['tmscore'].mean()
        am = a.groupby('target')['tmscore'].mean()
        common = bm.index.intersection(am.index)
        if len(common) > 5:
            ax.scatter(bm.loc[common], am.loc[common], alpha=0.4, s=15, color=color, label=label)
    ax.plot([0.5, 1], [0.5, 1], 'k--', linewidth=0.5, alpha=0.3)
    ax.set_xlabel('Before AMBER (TM-score)')
    ax.set_ylabel('After AMBER (TM-score)')
    ax.set_title('TM-score: No Change', fontsize=12, fontweight='bold')
    ax.set_xlim(0.5, 1.0)
    ax.set_ylim(0.5, 1.0)
    ax.set_aspect('equal')
    ax.legend(fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Panel 2: MolProbity change (before vs after AMBER)
    ax = axes[1]
    for src_before, src_after, color, label in [
        ('af_unrelaxed', 'amber_af', '#E69F00', 'AF→AMBER(AF)'),
        ('boltz', 'amber_boltz', '#D55E00', 'Boltz→AMBER(Boltz)'),
    ]:
        b = data[(data['pipeline'] == pipe) & (data['source'] == src_before)]
        a = data[(data['pipeline'] == pipe) & (data['source'] == src_after)]
        bm = b.groupby('target')['molprobity_score'].mean()
        am = a.groupby('target')['molprobity_score'].mean()
        common = bm.index.intersection(am.index)
        if len(common) > 5:
            ax.scatter(bm.loc[common], am.loc[common], alpha=0.4, s=15, color=color, label=label)
    ax.plot([0, 4], [0, 4], 'k--', linewidth=0.5, alpha=0.3)
    ax.set_xlabel('Before AMBER (MP Score)')
    ax.set_ylabel('After AMBER (MP Score)')
    ax.set_title('MolProbity Score: Large Improvement', fontsize=12, fontweight='bold')
    ax.set_xlim(0, 3.5)
    ax.set_ylim(0, 3.5)
    ax.set_aspect('equal')
    ax.legend(fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.suptitle('AMBER Dual Effect: Zero TM Change, Large MolProbity Improvement',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    save_fig(fig, 'fig14_amber_dual_effect')


def fig15_crystal_vs_predicted(data):
    """Crystal structures have worse MolProbity than predictions."""
    fig, ax = plt.subplots(figsize=(10, 6))

    pipe = 'blue'
    metrics_to_plot = ['clashscore', 'molprobity_score', 'rota_outliers']
    metric_labels = ['Clashscore', 'MP Score', 'Rota Out%']

    per_target = data[data['pipeline'] == pipe].groupby(['source', 'target'])[metrics_to_plot].mean().reset_index()

    x = np.arange(len(metrics_to_plot))
    width = 0.12
    offsets = np.linspace(-2.5*width, 2.5*width, len(SOURCE_ORDER))

    for i, src in enumerate(SOURCE_ORDER):
        src_data = per_target[per_target['source'] == src]
        if len(src_data) == 0:
            continue
        means = [src_data[m].mean() for m in metrics_to_plot]
        stds = [src_data[m].std() for m in metrics_to_plot]
        ax.bar(x + offsets[i], means, width, yerr=stds, label=SOURCE_LABELS[src].replace('\n', ' '),
               color=SRC_COLORS[src], edgecolor='black', linewidth=0.3, capsize=2,
               error_kw={'linewidth': 0.5})

    ax.set_xticks(x)
    ax.set_xticklabels(metric_labels, fontsize=11)
    ax.set_ylabel('Value (lower = better)')
    ax.set_title('MolProbity: Crystal vs Predicted Structures\n(Crystal has worst rotamers and many clashes)',
                fontsize=13, fontweight='bold')
    ax.legend(fontsize=8, ncol=3, loc='upper right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    save_fig(fig, 'fig15_crystal_vs_predicted')


def main():
    print("Loading MolProbity data...")
    data = load_data()
    print(f"  {len(data)} rows, {data['target'].nunique()} targets")

    os.makedirs(FIGDIR, exist_ok=True)

    print("\nGenerating MolProbity figures...")
    fig12_molprobity_by_source(data)
    fig13_amber_clashscore_effect(data)
    fig14_tm_vs_molprobity(data)
    fig15_crystal_vs_predicted(data)

    print("\nAll MolProbity figures generated!")


if __name__ == '__main__':
    main()
