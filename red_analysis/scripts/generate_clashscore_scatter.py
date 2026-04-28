#!/usr/bin/env python3
"""
generate_clashscore_scatter.py - Red Analysis Pipeline

Initial vs Final clashscore scatter: AMBER relaxed vs non-AMBER.

4-panel figure:
  Panel 1: Crystal (reference - no relaxation, plotted on identity)
  Panel 2: AF unrelaxed → AF relaxed (built-in AMBER)
  Panel 3: AF unrelaxed → AMBER(AF) (standalone AMBER)
  Panel 4: Boltz → AMBER(Boltz) (standalone AMBER)

Each point = one target (per-target mean clashscore).
Points below the diagonal = improvement after AMBER.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

METRICS_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"
FIGDIR = "/data/p_csb_meiler/agarwm5/red_analysis/figures"


def save_fig(fig, name):
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(FIGDIR, f"{name}.{fmt}"),
                    dpi=300 if fmt == 'png' else None, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {name}")


def load_data():
    mp = pd.read_csv(os.path.join(METRICS_DIR, "combined_molprobity.tsv"), sep='\t')
    mp['clashscore'] = pd.to_numeric(mp['clashscore'], errors='coerce')
    return mp


def clashscore_initial_vs_final(data):
    """4-panel Initial vs Final clashscore scatter. Both pipelines."""

    panel_configs = [
        ('Crystal\n(no relaxation)', 'crystal', 'crystal', '#e74c3c'),
        ('AF unrelaxed → AF relaxed\n(built-in AMBER)', 'af_unrelaxed', 'af_relaxed', '#1a1a2e'),
        ('AF unrelaxed → AMBER(AF)\n(standalone AMBER)', 'af_unrelaxed', 'amber_af', '#3498db'),
        ('Boltz → AMBER(Boltz)\n(standalone AMBER)', 'boltz', 'amber_boltz', '#2ecc71'),
    ]

    for pipe in ['blue', 'green']:
        pipe_label = pipe.capitalize()
        pipe_data = data[data['pipeline'] == pipe]

        fig, axes = plt.subplots(1, 4, figsize=(16, 4), sharex=True, sharey=True)
        fig.subplots_adjust(wspace=0.08)

        # Determine global max for consistent axes
        all_clash = pipe_data.groupby(['source', 'target'])['clashscore'].mean()
        fixed_max = min(55, np.nanpercentile(all_clash.values, 99) * 1.2)

        for idx, (title, src_before, src_after, color) in enumerate(panel_configs):
            ax = axes[idx]

            before = pipe_data[pipe_data['source'] == src_before].groupby('target')['clashscore'].mean()
            after = pipe_data[pipe_data['source'] == src_after].groupby('target')['clashscore'].mean()
            common = before.index.intersection(after.index)

            if len(common) == 0:
                ax.set_title(title, fontsize=9)
                continue

            x = before.loc[common].values
            y = after.loc[common].values
            mask = np.isfinite(x) & np.isfinite(y)
            x, y = x[mask], y[mask]

            # Color by improvement
            improved = y < x
            worse = y >= x

            if improved.any():
                ax.scatter(x[improved], y[improved], color=color, alpha=0.6,
                          s=25, edgecolors='white', linewidths=0.3,
                          label=f'Improved ({improved.sum()})', zorder=3)
            if worse.any():
                ax.scatter(x[worse], y[worse], color='#cccccc', alpha=0.5,
                          s=25, edgecolors='white', linewidths=0.3,
                          label=f'Worse ({worse.sum()})', zorder=2)

            # Identity line
            ax.plot([0, fixed_max], [0, fixed_max], 'k--', alpha=0.4, lw=1)

            # Stats annotation
            if src_before != src_after:
                delta = np.mean(y - x)
                pct_improved = improved.sum() / len(x) * 100
                ax.text(0.95, 0.05,
                        f'Δ = {delta:+.1f}\n{pct_improved:.0f}% improved\nn = {len(x)}',
                        transform=ax.transAxes, fontsize=8,
                        ha='right', va='bottom',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                  alpha=0.8, edgecolor='gray', linewidth=0.5))
            else:
                ax.text(0.95, 0.05, f'n = {len(x)}\n(reference)',
                        transform=ax.transAxes, fontsize=8,
                        ha='right', va='bottom',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                  alpha=0.8, edgecolor='gray', linewidth=0.5))

            ax.set_title(title, fontsize=9, fontweight='bold')
            ax.set_xlim(0, fixed_max)
            ax.set_ylim(0, fixed_max)
            ax.set_aspect('equal', adjustable='box')
            ax.legend(fontsize=7, loc='upper left', framealpha=0.8)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            if idx == 0:
                ax.set_ylabel('Final Clashscore', fontsize=11)

        # Shared x-label
        fig.text(0.5, -0.02, 'Initial Clashscore', ha='center', fontsize=11)

        fig.suptitle(f'AMBER Relaxation Effect on Clashscore ({pipe_label} pipeline)\n'
                     f'Points below diagonal = improvement',
                     fontsize=13, fontweight='bold', y=1.05)

        save_fig(fig, f'fig21_clashscore_initial_vs_final_{pipe}')


def clashscore_amber_comparison(data):
    """2-panel comparison: AF AMBER vs Boltz AMBER side by side. Both pipelines."""

    for pipe in ['blue', 'green']:
        pipe_label = pipe.capitalize()
        pipe_data = data[data['pipeline'] == pipe]

        fig, axes = plt.subplots(1, 2, figsize=(10, 5))

        pairs = [
            ('AF unrelaxed → AMBER(AF)', 'af_unrelaxed', 'amber_af', '#3498db'),
            ('Boltz → AMBER(Boltz)', 'boltz', 'amber_boltz', '#2ecc71'),
        ]

        fixed_max = 55

        for idx, (title, src_before, src_after, color) in enumerate(pairs):
            ax = axes[idx]

            before = pipe_data[pipe_data['source'] == src_before].groupby('target')['clashscore'].mean()
            after = pipe_data[pipe_data['source'] == src_after].groupby('target')['clashscore'].mean()
            common = before.index.intersection(after.index)

            if len(common) == 0:
                continue

            x = before.loc[common].values
            y = after.loc[common].values
            mask = np.isfinite(x) & np.isfinite(y)
            x, y = x[mask], y[mask]

            # Scatter with error bars (std across models within target)
            ax.scatter(x, y, color=color, alpha=0.6, s=30,
                      edgecolors='black', linewidths=0.3, zorder=3)

            # Identity line
            ax.plot([0, fixed_max], [0, fixed_max], 'k--', alpha=0.4, lw=1)

            # Fill region below diagonal (improvement zone)
            ax.fill_between([0, fixed_max], [0, fixed_max], [0, 0],
                           alpha=0.05, color='green', label='Improvement zone')

            # Stats
            delta = np.mean(y - x)
            improved = (y < x).sum()
            total = len(x)
            median_before = np.median(x)
            median_after = np.median(y)

            stats_text = (
                f'Mean Δ = {delta:+.1f}\n'
                f'Median: {median_before:.1f} → {median_after:.1f}\n'
                f'Improved: {improved}/{total} ({improved/total*100:.0f}%)'
            )
            ax.text(0.95, 0.05, stats_text,
                    transform=ax.transAxes, fontsize=9,
                    ha='right', va='bottom',
                    bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                              alpha=0.9, edgecolor='gray'))

            ax.set_xlabel('Initial Clashscore', fontsize=11)
            if idx == 0:
                ax.set_ylabel('Final Clashscore (after AMBER)', fontsize=11)
            ax.set_title(title, fontsize=11, fontweight='bold')
            ax.set_xlim(0, fixed_max)
            ax.set_ylim(0, fixed_max)
            ax.set_aspect('equal', adjustable='box')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        fig.suptitle(f'AMBER Relaxation: Initial vs Final Clashscore ({pipe_label})\n'
                     f'257 BM5.5 targets - per-target averages',
                     fontsize=13, fontweight='bold')
        plt.tight_layout()
        save_fig(fig, f'fig22_amber_clashscore_comparison_{pipe}')


def main():
    os.makedirs(FIGDIR, exist_ok=True)

    print("Loading MolProbity data...")
    data = load_data()
    print(f"  {len(data)} rows, {data['target'].nunique()} targets")

    print("\nGenerating clashscore scatter figures...")
    clashscore_initial_vs_final(data)
    clashscore_amber_comparison(data)

    print("\nDone!")


if __name__ == '__main__':
    main()
