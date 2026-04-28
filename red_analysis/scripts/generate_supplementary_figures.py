#!/usr/bin/env python3
"""
generate_supplementary_figures.py - Red Analysis Pipeline

Supplementary figures for the paper.

Figures:
  figS1: Per-metric AMBER effect (all 7 MolProbity metrics, before/after)
  figS2: MolProbity distribution by source (violin plots, all metrics)
  figS3: Rosetta TM-score by input type (the protocol × source detail)
  figS4: Correlation matrix of MolProbity metrics
  figS5: Target-level MolProbity scatter (AF vs Boltz)
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
    'af_relaxed': 'AF2 (relaxed)',
    'af_unrelaxed': 'AF2 (unrelaxed)',
    'boltz': 'Boltz-1',
    'amber_af': 'AMBER(AF)',
    'amber_boltz': 'AMBER(Boltz)',
}


def save_fig(fig, name):
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(FIGDIR, f"{name}.{fmt}"),
                   dpi=300 if fmt == 'png' else None, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {name}")


def figS1_amber_all_metrics(data):
    """All 7 MolProbity metrics before/after AMBER. Both pipelines."""
    metrics = [
        ('clashscore', 'Clashscore', True),
        ('molprobity_score', 'MolProbity Score', True),
        ('rama_outliers', 'Rama Outliers %', True),
        ('rama_favored', 'Rama Favored %', False),
        ('rota_outliers', 'Rota Outliers %', True),
        ('rms_bonds', 'RMS Bonds (Å)', True),
        ('rms_angles', 'RMS Angles (°)', True),
    ]

    pairs = [
        ('af_unrelaxed', 'amber_af', 'AF→AMBER(AF)', '#E69F00'),
        ('boltz', 'amber_boltz', 'Boltz→AMBER(Boltz)', '#D55E00'),
    ]

    for pipe in ['blue', 'green']:
        pipe_label = pipe.capitalize()
        pipe_data = data[data['pipeline'] == pipe]

        fig, axes = plt.subplots(2, 4, figsize=(20, 10))
        axes = axes.flatten()

        for idx, (metric, label, lower_better) in enumerate(metrics):
            ax = axes[idx]

            for src_before, src_after, pair_label, color in pairs:
                before = pipe_data[pipe_data['source'] == src_before].groupby('target')[metric].mean()
                after = pipe_data[pipe_data['source'] == src_after].groupby('target')[metric].mean()
                common = before.index.intersection(after.index)

                if len(common) < 5:
                    continue

                delta = after.loc[common] - before.loc[common]

                ax.hist(delta.values, bins=30, alpha=0.5, color=color, label=pair_label,
                       edgecolor='white', linewidth=0.3)

            ax.axvline(x=0, color='black', linewidth=0.5, linestyle='--')
            ax.set_xlabel(f'Δ{label}')
            ax.set_ylabel('Count')
            direction = '←better' if lower_better else 'better→'
            ax.set_title(f'{label}\n({direction})', fontsize=10)
            ax.legend(fontsize=7)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        axes[7].set_visible(False)

        fig.suptitle(f'AMBER Effect on All MolProbity Metrics ({pipe_label})\n(Δ = After - Before)',
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        save_fig(fig, f'figS1_amber_all_metrics_{pipe}')


def figS2_molprobity_violins(data):
    """Violin plots for all MolProbity metrics by source. Both pipelines."""
    metrics = [
        ('clashscore', 'Clashscore'),
        ('molprobity_score', 'MP Score'),
        ('rama_favored', 'Rama Favored %'),
        ('rota_outliers', 'Rota Outliers %'),
        ('rms_bonds', 'RMS Bonds'),
        ('rms_angles', 'RMS Angles'),
    ]

    for pipe in ['blue', 'green']:
        pipe_label = pipe.capitalize()
        per_target = data[data['pipeline'] == pipe].groupby(['source', 'target'])[
            [m[0] for m in metrics]].mean().reset_index()

        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()

        for ax, (metric, label) in zip(axes, metrics):
            plot_data = []
            positions = []
            colors = []
            tick_labels = []

            for i, src in enumerate(SOURCE_ORDER):
                vals = per_target[per_target['source'] == src][metric].dropna().values
                if len(vals) > 0:
                    plot_data.append(vals)
                    positions.append(i)
                    colors.append(SRC_COLORS[src])
                    tick_labels.append(SOURCE_LABELS[src])

            if plot_data:
                vp = ax.violinplot(plot_data, positions=positions, showmedians=True, showextrema=False)
                for body, color in zip(vp['bodies'], colors):
                    body.set_facecolor(color)
                    body.set_alpha(0.7)
                vp['cmedians'].set_color('black')

            ax.set_xticks(positions)
            ax.set_xticklabels(tick_labels, fontsize=8, rotation=30, ha='right')
            ax.set_ylabel(label)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        fig.suptitle(f'MolProbity Metric Distributions by Source\n({pipe_label} pipeline, per-target averages)',
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        save_fig(fig, f'figS2_molprobity_violins_{pipe}')


def figS3_rosetta_by_source(data_rosetta):
    """Rosetta TM-score by input type for each protocol. Both pipelines."""
    if data_rosetta is None:
        print("  SKIP figS3: no Rosetta data")
        return

    PROTOCOL_ORDER = ['cartesian_beta', 'cartesian_ref15', 'dualspace_beta',
                      'dualspace_ref15', 'normal_beta', 'normal_ref15']
    PROTOCOL_LABELS = {
        'cartesian_beta': 'Cart β',
        'cartesian_ref15': 'Cart REF15',
        'dualspace_beta': 'Dual β',
        'dualspace_ref15': 'Dual REF15',
        'normal_beta': 'Norm β',
        'normal_ref15': 'Norm REF15',
    }
    src_order = ['af_relaxed', 'af_unrelaxed', 'boltz', 'amber_af', 'amber_boltz', 'crystal']

    for pipe in ['blue', 'green']:
        pipe_label = pipe.capitalize()
        pipe_data = data_rosetta[data_rosetta['pipeline'] == pipe]

        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()

        for idx, proto in enumerate(PROTOCOL_ORDER):
            ax = axes[idx]
            proto_data = pipe_data[pipe_data['protocol'] == proto]
            per_target = proto_data.groupby(['source', 'target'])['tmscore'].mean().reset_index()

            plot_data = []
            labels = []
            colors = []
            for src in src_order:
                vals = per_target[per_target['source'] == src]['tmscore'].dropna().values
                if len(vals) > 0:
                    plot_data.append(vals)
                    labels.append(SOURCE_LABELS.get(src, src))
                    colors.append(SRC_COLORS.get(src, '#888888'))

            if plot_data:
                bp = ax.boxplot(plot_data, tick_labels=labels, patch_artist=True, widths=0.6,
                               medianprops={'color': 'black', 'linewidth': 1.5})
                for patch, color in zip(bp['boxes'], colors):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)

            ax.set_title(PROTOCOL_LABELS[proto], fontsize=11, fontweight='bold')
            ax.set_ylabel('TM-score')
            ax.tick_params(axis='x', labelsize=7, rotation=30)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        fig.suptitle(f'Rosetta TM-score by Input Type\n({pipe_label} pipeline, per-target averages)',
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        save_fig(fig, f'figS3_rosetta_by_source_{pipe}')


def figS4_mp_correlation_matrix(data):
    """Correlation matrix of MolProbity metrics. Both pipelines."""
    metrics = ['clashscore', 'molprobity_score', 'rama_outliers', 'rama_favored',
               'rota_outliers', 'rms_bonds', 'rms_angles']
    labels = ['Clashscore', 'MP Score', 'Rama Out%', 'Rama Fav%',
              'Rota Out%', 'RMS Bond', 'RMS Angle']

    for pipe in ['blue', 'green']:
        pipe_label = pipe.capitalize()
        per_target = data[data['pipeline'] == pipe].groupby(['source', 'target'])[metrics].mean().reset_index()

        corr = per_target[metrics].corr()
        corr.index = labels
        corr.columns = labels

        fig, ax = plt.subplots(figsize=(8, 7))
        mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
        sns.heatmap(corr, mask=mask, annot=True, fmt='.2f', cmap='RdBu_r',
                   vmin=-1, vmax=1, center=0, ax=ax, linewidths=0.5,
                   square=True, cbar_kws={'shrink': 0.8})
        ax.set_title(f'MolProbity Metric Correlations\n({pipe_label} pipeline)', fontsize=13, fontweight='bold')
        plt.tight_layout()
        save_fig(fig, f'figS4_mp_correlation_{pipe}')


def figS5_af_vs_boltz_mp(data):
    """AF vs Boltz MolProbity scatter (target-level). Both pipelines."""
    metrics = [('clashscore', 'Clashscore'), ('molprobity_score', 'MP Score'),
               ('rota_outliers', 'Rota Outliers %')]

    for pipe in ['blue', 'green']:
        pipe_label = pipe.capitalize()
        per_target = data[data['pipeline'] == pipe].groupby(['source', 'target'])[
            [m[0] for m in metrics]].mean().reset_index()

        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        for ax, (metric, label) in zip(axes, metrics):
            af = per_target[per_target['source'] == 'af_relaxed'].set_index('target')[metric]
            boltz = per_target[per_target['source'] == 'boltz'].set_index('target')[metric]
            common = af.index.intersection(boltz.index)

            if len(common) > 5:
                x = af.loc[common].values
                y = boltz.loc[common].values
                mask = np.isfinite(x) & np.isfinite(y)
                ax.scatter(x[mask], y[mask], alpha=0.5, s=20, color='#0072B2', edgecolors='black', linewidths=0.3)

                max_val = max(x[mask].max(), y[mask].max()) * 1.1
                ax.plot([0, max_val], [0, max_val], 'k--', linewidth=0.5, alpha=0.3)

                r, p = stats.pearsonr(x[mask], y[mask])
                ax.text(0.05, 0.95, f'r = {r:.3f}', transform=ax.transAxes,
                       fontsize=10, verticalalignment='top')

                af_better = (x[mask] < y[mask]).sum()
                boltz_better = (y[mask] < x[mask]).sum()
                ax.text(0.05, 0.85, f'AF better: {af_better}\nBoltz better: {boltz_better}',
                       transform=ax.transAxes, fontsize=8, verticalalignment='top')

            ax.set_xlabel(f'AF2 {label}')
            ax.set_ylabel(f'Boltz {label}')
            ax.set_title(label)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        fig.suptitle(f'AF2 vs Boltz MolProbity Comparison\n(per-target, {pipe_label} pipeline)',
                    fontsize=13, fontweight='bold')
        plt.tight_layout()
        save_fig(fig, f'figS5_af_vs_boltz_mp_{pipe}')


def main():
    os.makedirs(FIGDIR, exist_ok=True)

    print("Loading data...")

    # MolProbity
    mp_path = os.path.join(METRICS_DIR, "combined_molprobity.tsv")
    mp = pd.read_csv(mp_path, sep='\t')
    for col in ['clashscore', 'rama_outliers', 'rama_favored', 'rota_outliers',
                'molprobity_score', 'cbeta_outliers', 'rms_bonds', 'rms_angles']:
        mp[col] = pd.to_numeric(mp[col], errors='coerce')
    print(f"  MolProbity: {len(mp)} rows, {mp['target'].nunique()} targets")

    # Rosetta
    rt_path = os.path.join(METRICS_DIR, "combined_rosetta_tmscore.tsv")
    rt = None
    if os.path.exists(rt_path):
        rt = pd.read_csv(rt_path, sep='\t')
        for col in ['tmscore', 'rmsd', 'gdtts', 'gdtha']:
            rt[col] = pd.to_numeric(rt[col], errors='coerce')
        print(f"  Rosetta: {len(rt)} rows, {rt['target'].nunique()} targets")

    print("\nGenerating supplementary figures...")
    figS1_amber_all_metrics(mp)
    figS2_molprobity_violins(mp)
    figS3_rosetta_by_source(rt)
    figS4_mp_correlation_matrix(mp)
    figS5_af_vs_boltz_mp(mp)

    print("\nAll supplementary figures generated!")


if __name__ == '__main__':
    main()
