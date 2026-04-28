#!/usr/bin/env python3
"""
generate_figures.py - Red Analysis Pipeline

Publication-quality figures for the BM5.5 protein relaxation benchmark.

Figures generated:
  1. Violin/box plot: TM-score by source type (AF vs Boltz vs AMBER)
  2. Violin/box plot: RMSD by source type
  3. Scatter: Blue vs Green TM-score (reproducibility)
  4. Heatmap: per-target TM-score across source types
  5. AMBER effect: paired difference plots (before vs after relaxation)
  6. Outlier analysis: targets with TM-score < 0.8

COMMENT: Using matplotlib + seaborn for publication figures. We could
use plotly for interactive exploration but publication needs static PDFs.
Nature/Science journals prefer vector graphics (PDF/SVG).

COMMENT: Color scheme follows colorblind-friendly palette (Wong 2011).
"""

import os
import sys
import numpy as np

# Try importing full stack, fall back to basic
try:
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend for cluster
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_FULL_STACK = True
except ImportError:
    HAS_FULL_STACK = False
    print("WARNING: pandas/matplotlib/seaborn not available. Using basic mode.")

METRICS_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"
FIGURES_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/figures"
TABLES_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/tables"

# Colorblind-friendly palette (Wong 2011, Nature Methods)
COLORS = {
    'af_relaxed': '#0072B2',      # Blue
    'af_unrelaxed': '#56B4E9',    # Sky blue
    'boltz': '#E69F00',           # Orange
    'amber_af': '#009E73',        # Green
    'amber_boltz': '#D55E00',     # Vermillion
    'crystal': '#CC79A7',         # Pink
}

# Human-readable labels
LABELS = {
    'af_relaxed': 'AF2 (AMBER)',
    'af_unrelaxed': 'AF2 (unrelaxed)',
    'boltz': 'Boltz-1',
    'amber_af': 'AMBER(AF2)',
    'amber_boltz': 'AMBER(Boltz)',
    'crystal': 'Crystal',
}

SOURCE_ORDER = ['af_unrelaxed', 'af_relaxed', 'amber_af', 'boltz', 'amber_boltz']


def load_data():
    """Load combined TM-score data."""
    combined = os.path.join(METRICS_DIR, "combined_tmscore.tsv")
    if not os.path.exists(combined):
        print("ERROR: Run aggregate_tmscore.py first")
        sys.exit(1)

    df = pd.read_csv(combined, sep='\t')
    # Convert numeric columns
    for col in ['rmsd', 'tmscore', 'gdtts', 'gdtha', 'aligned_len', 'seq_len']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    return df


def fig1_tmscore_by_source(df):
    """
    Figure 1: TM-score distribution by prediction source.

    This is the main result figure. Shows how well each prediction method
    reconstructs the crystal structure.

    COMMENT: Using violin plots rather than box plots because they show
    the distribution shape. Many of our TM-scores cluster near 1.0, so
    the violin shape will reveal the heavy tail of difficult targets.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, pipeline in zip(axes, ['blue', 'green']):
        sub = df[df['pipeline'] == pipeline].copy()

        # Mean of 5 models per target per source
        means = sub.groupby(['target', 'source'])['tmscore'].mean().reset_index()

        # Filter to standard sources and order
        means = means[means['source'].isin(SOURCE_ORDER)]
        means['source'] = pd.Categorical(means['source'], categories=SOURCE_ORDER, ordered=True)

        sns.violinplot(
            data=means, x='source', y='tmscore',
            palette=[COLORS[s] for s in SOURCE_ORDER],
            inner='box', cut=0, ax=ax
        )

        ax.set_xticklabels([LABELS[s] for s in SOURCE_ORDER], rotation=30, ha='right')
        ax.set_ylabel('TM-score to Crystal')
        ax.set_xlabel('')
        ax.set_title(f'{"Blue" if pipeline == "blue" else "Green"} Pipeline')
        ax.set_ylim(0.3, 1.02)
        ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, label='Random threshold')

        # Add mean markers
        for i, src in enumerate(SOURCE_ORDER):
            mean_val = means[means['source'] == src]['tmscore'].mean()
            ax.scatter(i, mean_val, color='black', marker='D', s=30, zorder=5)

    plt.suptitle('Structural Similarity to Crystal (TM-score)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    outf = os.path.join(FIGURES_DIR, "fig1_tmscore_by_source.pdf")
    plt.savefig(outf, dpi=300, bbox_inches='tight')
    plt.savefig(outf.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Figure 1 -> {outf}")


def fig2_rmsd_by_source(df):
    """
    Figure 2: RMSD distribution by prediction source.

    COMMENT: RMSD is more intuitive than TM-score for most audiences.
    But it's sensitive to a few badly placed residues (one disordered loop
    can dominate). TM-score is more robust.

    QUESTION: Should we use log scale for RMSD y-axis? Some targets have
    RMSD > 15Å which compresses the interesting range. Using log for now.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, pipeline in zip(axes, ['blue', 'green']):
        sub = df[df['pipeline'] == pipeline].copy()
        means = sub.groupby(['target', 'source'])['rmsd'].mean().reset_index()
        means = means[means['source'].isin(SOURCE_ORDER)]
        means['source'] = pd.Categorical(means['source'], categories=SOURCE_ORDER, ordered=True)

        sns.boxplot(
            data=means, x='source', y='rmsd',
            palette=[COLORS[s] for s in SOURCE_ORDER],
            showfliers=True, flierprops={'marker': 'o', 'markersize': 3, 'alpha': 0.5},
            ax=ax
        )

        ax.set_xticklabels([LABELS[s] for s in SOURCE_ORDER], rotation=30, ha='right')
        ax.set_ylabel('RMSD to Crystal (Å)')
        ax.set_xlabel('')
        ax.set_title(f'{"Blue" if pipeline == "blue" else "Green"} Pipeline')
        ax.set_yscale('log')
        ax.axhline(y=2.0, color='gray', linestyle='--', alpha=0.5)

    plt.suptitle('RMSD to Crystal Reference', fontsize=14, fontweight='bold')
    plt.tight_layout()
    outf = os.path.join(FIGURES_DIR, "fig2_rmsd_by_source.pdf")
    plt.savefig(outf, dpi=300, bbox_inches='tight')
    plt.savefig(outf.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Figure 2 -> {outf}")


def fig3_blue_vs_green(df):
    """
    Figure 3: Blue vs Green reproducibility scatter.

    This is the key reproducibility figure. Each point is one (target, source)
    pair. If pipelines are reproducible, points cluster on the diagonal.

    COMMENT: We expect systematic scatter because AF/Boltz are stochastic.
    The question is HOW MUCH scatter. If mean-of-5 is on the diagonal,
    that validates the benchmark design.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, metric, label in zip(axes, ['tmscore', 'rmsd'], ['TM-score', 'RMSD (Å)']):
        # Get mean-of-5 per target per source per pipeline
        means = df.groupby(['target', 'source', 'pipeline'])[metric].mean().reset_index()

        blue = means[means['pipeline'] == 'blue'].drop(columns='pipeline')
        green = means[means['pipeline'] == 'green'].drop(columns='pipeline')

        merged = blue.merge(green, on=['target', 'source'], suffixes=('_blue', '_green'))

        for src in SOURCE_ORDER:
            sub = merged[merged['source'] == src]
            ax.scatter(
                sub[f'{metric}_blue'], sub[f'{metric}_green'],
                c=COLORS[src], label=LABELS[src], alpha=0.4, s=15, edgecolors='none'
            )

        # Diagonal
        lims = [ax.get_xlim(), ax.get_ylim()]
        lo = min(lims[0][0], lims[1][0])
        hi = max(lims[0][1], lims[1][1])
        ax.plot([lo, hi], [lo, hi], 'k--', alpha=0.5, linewidth=1)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect('equal')

        ax.set_xlabel(f'Blue {label}')
        ax.set_ylabel(f'Green {label}')
        ax.legend(fontsize=8, loc='lower right')

        # Compute Pearson correlation
        from scipy import stats
        r, p = stats.pearsonr(merged[f'{metric}_blue'].dropna(), merged[f'{metric}_green'].dropna())
        ax.text(0.05, 0.95, f'r = {r:.4f}', transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.suptitle('Blue vs Green: Pipeline Reproducibility', fontsize=14, fontweight='bold')
    plt.tight_layout()
    outf = os.path.join(FIGURES_DIR, "fig3_blue_vs_green.pdf")
    plt.savefig(outf, dpi=300, bbox_inches='tight')
    plt.savefig(outf.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Figure 3 -> {outf}")


def fig4_amber_effect(df):
    """
    Figure 4: Effect of AMBER relaxation on structural similarity.

    Paired comparison: for each target, what happens to TM-score when you
    apply AMBER relaxation?

    AF unrelaxed → AF relaxed (built-in AMBER)
    AF unrelaxed → standalone AMBER(AF)
    Boltz → AMBER(Boltz)

    Generates separate figures for each pipeline (Blue + Green).
    """
    for pipeline in ['blue', 'green']:
        pipe_label = pipeline.capitalize()
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        pairs = [
            ('af_unrelaxed', 'af_relaxed', 'AF2: Built-in AMBER'),
            ('af_unrelaxed', 'amber_af', 'AF2: Standalone AMBER'),
            ('boltz', 'amber_boltz', 'Boltz: Standalone AMBER'),
        ]

        for ax, (before_src, after_src, title) in zip(axes, pairs):
            sub = df[df['pipeline'] == pipeline].copy()

            before = sub[sub['source'] == before_src].groupby('target')['tmscore'].mean()
            after = sub[sub['source'] == after_src].groupby('target')['tmscore'].mean()

            common = before.index.intersection(after.index)
            diff = after.loc[common] - before.loc[common]

            ax.hist(diff, bins=50, color='steelblue', edgecolor='black', linewidth=0.5, alpha=0.7)
            ax.axvline(x=0, color='red', linestyle='--', linewidth=2)
            ax.axvline(x=diff.mean(), color='orange', linestyle='-', linewidth=2, label=f'Mean: {diff.mean():.4f}')
            ax.axvline(x=diff.median(), color='green', linestyle='-', linewidth=2, label=f'Median: {diff.median():.4f}')

            ax.set_xlabel('ΔTM-score (after - before)')
            ax.set_ylabel('Number of targets')
            ax.set_title(title)
            ax.legend(fontsize=8)

            n_better = (diff > 0.001).sum()
            n_worse = (diff < -0.001).sum()
            n_same = len(diff) - n_better - n_worse
            ax.text(0.05, 0.85, f'Better: {n_better}\nWorse: {n_worse}\nSame: {n_same}',
                    transform=ax.transAxes, fontsize=9, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

        plt.suptitle(f'Effect of AMBER Relaxation on TM-score ({pipe_label} pipeline)',
                     fontsize=14, fontweight='bold')
        plt.tight_layout()
        outf = os.path.join(FIGURES_DIR, f"fig4_amber_effect_{pipeline}.pdf")
        plt.savefig(outf, dpi=300, bbox_inches='tight')
        plt.savefig(outf.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Figure 4 ({pipe_label}) -> {outf}")


def fig5_af_vs_boltz(df):
    """
    Figure 5: AF2 vs Boltz-1 head-to-head comparison.

    Scatter plot: each target is a point. x = AF2 TM-score, y = Boltz TM-score.
    Points above diagonal = Boltz better. Below = AF2 better.

    Generates separate figures for each pipeline (Blue + Green).
    """
    for pipeline in ['blue', 'green']:
        pipe_label = pipeline.capitalize()
        fig, ax = plt.subplots(figsize=(8, 8))

        sub = df[df['pipeline'] == pipeline].copy()

        af = sub[sub['source'] == 'af_unrelaxed'].groupby('target')['tmscore'].mean()
        boltz = sub[sub['source'] == 'boltz'].groupby('target')['tmscore'].mean()

        common = af.index.intersection(boltz.index)

        ax.scatter(af.loc[common], boltz.loc[common], c='steelblue', alpha=0.5, s=20, edgecolors='none')
        ax.plot([0, 1], [0, 1], 'k--', alpha=0.5)

        ax.set_xlabel('AF2 TM-score to Crystal', fontsize=12)
        ax.set_ylabel('Boltz-1 TM-score to Crystal', fontsize=12)
        ax.set_title(f'AF2 vs Boltz-1: Per-Target Comparison ({pipe_label})',
                     fontsize=14, fontweight='bold')
        ax.set_xlim(0.3, 1.02)
        ax.set_ylim(0.3, 1.02)
        ax.set_aspect('equal')

        n_af_wins = (af.loc[common] > boltz.loc[common] + 0.01).sum()
        n_boltz_wins = (boltz.loc[common] > af.loc[common] + 0.01).sum()
        n_ties = len(common) - n_af_wins - n_boltz_wins

        ax.text(0.05, 0.95,
                f'AF2 better: {n_af_wins}\nBoltz better: {n_boltz_wins}\nTied: {n_ties}\n'
                f'AF2 mean: {af.loc[common].mean():.4f}\nBoltz mean: {boltz.loc[common].mean():.4f}',
                transform=ax.transAxes, fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

        plt.tight_layout()
        outf = os.path.join(FIGURES_DIR, f"fig5_af_vs_boltz_{pipeline}.pdf")
        plt.savefig(outf, dpi=300, bbox_inches='tight')
        plt.savefig(outf.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Figure 5 ({pipe_label}) -> {outf}")


def fig6_outlier_analysis(df):
    """
    Figure 6: Outlier targets with low TM-score.

    Bar chart showing the worst-performing targets and their TM-scores
    across different source types.

    Generates separate figures for each pipeline (Blue + Green).
    """
    for pipeline in ['blue', 'green']:
        pipe_label = pipeline.capitalize()
        sub = df[df['pipeline'] == pipeline].copy()

        means = sub.groupby(['target', 'source'])['tmscore'].mean().reset_index()
        pivot = means.pivot(index='target', columns='source', values='tmscore')

        outliers = pivot[pivot.min(axis=1) < 0.8].copy()
        available_srcs = [s for s in SOURCE_ORDER if s in outliers.columns]
        outliers = outliers[available_srcs].sort_values(available_srcs[0] if available_srcs else outliers.columns[0])

        if len(outliers) == 0:
            print(f"No outlier targets with TM < 0.8 ({pipe_label})")
            continue

        fig, ax = plt.subplots(figsize=(max(12, len(outliers) * 0.5), 6))

        x = np.arange(len(outliers))
        width = 0.15

        for i, src in enumerate(SOURCE_ORDER):
            if src in outliers.columns:
                ax.bar(x + i * width, outliers[src], width, label=LABELS[src],
                       color=COLORS[src], edgecolor='black', linewidth=0.5)

        ax.set_xticks(x + width * 2)
        ax.set_xticklabels(outliers.index, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('TM-score to Crystal')
        ax.set_title(f'Outlier Targets ({pipe_label}, TM < 0.8 for any source, n={len(outliers)})',
                     fontsize=14, fontweight='bold')
        ax.legend(fontsize=9)
        ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
        ax.axhline(y=0.8, color='red', linestyle='--', alpha=0.3)

        plt.tight_layout()
        outf = os.path.join(FIGURES_DIR, f"fig6_outliers_{pipeline}.pdf")
        plt.savefig(outf, dpi=300, bbox_inches='tight')
        plt.savefig(outf.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Figure 6 ({pipe_label}) -> {outf} ({len(outliers)} outlier targets)")


def write_latex_table(df):
    """
    Generate LaTeX-formatted summary tables for both pipelines.
    """
    for pipeline in ['blue', 'green']:
        pipe_label = pipeline.capitalize()
        sub = df[df['pipeline'] == pipeline].copy()

        results = []
        for src in SOURCE_ORDER:
            s = sub[sub['source'] == src]
            means = s.groupby('target')['tmscore'].mean()
            rmsd_means = s.groupby('target')['rmsd'].mean()
            gdtts_means = s.groupby('target')['gdtts'].mean()

            if len(means) == 0:
                continue

            results.append({
                'Source': LABELS[src],
                'TM-score': f'{means.mean():.3f} ± {means.std():.3f}',
                'RMSD (Å)': f'{rmsd_means.mean():.2f} ± {rmsd_means.std():.2f}',
                'GDT-TS': f'{gdtts_means.mean():.3f} ± {gdtts_means.std():.3f}',
                'N targets': len(means),
            })

        result_df = pd.DataFrame(results)
        latex = result_df.to_latex(index=False, escape=False)

        outf = os.path.join(TABLES_DIR, f"table1_summary_{pipeline}.tex")
        with open(outf, 'w') as f:
            f.write(latex)
        print(f"Table 1 ({pipe_label} LaTeX) -> {outf}")

        result_df.to_csv(os.path.join(TABLES_DIR, f"table1_summary_{pipeline}.tsv"), sep='\t', index=False)


if __name__ == '__main__':
    if not HAS_FULL_STACK:
        print("ERROR: Need pandas/matplotlib/seaborn. Activate red_analysis conda env.")
        sys.exit(1)

    os.makedirs(FIGURES_DIR, exist_ok=True)

    print("Loading data...")
    df = load_data()
    print(f"Loaded {len(df)} rows")

    print("\nGenerating figures...")
    fig1_tmscore_by_source(df)
    fig2_rmsd_by_source(df)
    fig3_blue_vs_green(df)
    fig4_amber_effect(df)
    fig5_af_vs_boltz(df)
    fig6_outlier_analysis(df)

    print("\nGenerating tables...")
    write_latex_table(df)

    print("\n=== All figures generated ===")
