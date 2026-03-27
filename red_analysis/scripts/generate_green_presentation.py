#!/usr/bin/env python3
"""
generate_green_presentation.py — Red Analysis Pipeline

Generates all presentation figures for the Green pipeline.
For every MolProbity metric, produces:
  1. 4-panel Initial vs Final scatter (AMBER effect)
  2. Per-target bar chart (sorted by initial, showing before/after AMBER)

Output: green_data_analysis/ folder
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

METRICS_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"
OUTDIR = "/data/p_csb_meiler/agarwm5/red_analysis/green_data_analysis"

# All MolProbity metrics to generate figures for
METRICS = [
    ('clashscore', 'Clashscore', True),
    ('molprobity_score', 'MolProbity Score', True),
    ('rama_outliers', 'Ramachandran Outliers %', True),
    ('rama_favored', 'Ramachandran Favored %', False),
    ('rota_outliers', 'Rotamer Outliers %', True),
    ('rms_bonds', 'RMS Bonds (Å)', True),
    ('rms_angles', 'RMS Angles (°)', True),
    ('cbeta_outliers', 'C-beta Outliers', True),
]


def save_fig(fig, name):
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(OUTDIR, f"{name}.{fmt}"),
                    dpi=300 if fmt == 'png' else None, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {name}")


def load_data():
    mp = pd.read_csv(os.path.join(METRICS_DIR, "combined_molprobity.tsv"), sep='\t')
    for col in [m[0] for m in METRICS]:
        mp[col] = pd.to_numeric(mp[col], errors='coerce')
    # Green pipeline only
    mp = mp[mp['pipeline'] == 'green']
    return mp


def make_scatter(data, metric, label, lower_better):
    """
    4-panel Initial vs Final scatter for one metric.

    Panels:
      1. Crystal (reference — on identity line)
      2. AF unrelaxed → AF relaxed (built-in AMBER)
      3. AF unrelaxed → AMBER(AF) (standalone AMBER)
      4. Boltz → AMBER(Boltz) (standalone AMBER)
    """
    panel_configs = [
        ('Crystal\n(reference)', 'crystal', 'crystal', '#e74c3c'),
        ('AF unrelaxed → AF relaxed\n(built-in AMBER)', 'af_unrelaxed', 'af_relaxed', '#1a1a2e'),
        ('AF unrelaxed → AMBER(AF)\n(standalone AMBER)', 'af_unrelaxed', 'amber_af', '#3498db'),
        ('Boltz → AMBER(Boltz)\n(standalone AMBER)', 'boltz', 'amber_boltz', '#2ecc71'),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(16, 4), sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0.08)

    # Global max for axis limits
    all_vals = data.groupby(['source', 'target'])[metric].mean()
    valid = all_vals.dropna()
    if len(valid) == 0:
        plt.close(fig)
        return
    p99 = np.nanpercentile(valid.values, 99)
    fixed_max = max(p99 * 1.2, 1e-6)

    for idx, (title, src_before, src_after, color) in enumerate(panel_configs):
        ax = axes[idx]

        before = data[data['source'] == src_before].groupby('target')[metric].mean()
        after = data[data['source'] == src_after].groupby('target')[metric].mean()
        common = before.index.intersection(after.index)

        if len(common) == 0:
            ax.set_title(title, fontsize=9)
            ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center')
            continue

        x = before.loc[common].values
        y = after.loc[common].values
        mask = np.isfinite(x) & np.isfinite(y)
        x, y = x[mask], y[mask]

        if len(x) == 0:
            ax.set_title(title, fontsize=9)
            continue

        # Color by improvement direction
        if lower_better:
            improved = y < x
        else:
            improved = y > x
        worse = ~improved

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

        # Stats
        if src_before != src_after:
            delta = np.mean(y - x)
            pct = improved.sum() / len(x) * 100
            ax.text(0.95, 0.05,
                    f'Δ = {delta:+.2f}\n{pct:.0f}% improved\nn = {len(x)}',
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
            ax.set_ylabel(f'Final {label}', fontsize=11)

    fig.text(0.5, -0.02, f'Initial {label}', ha='center', fontsize=11)
    direction = '(lower = better)' if lower_better else '(higher = better)'
    fig.suptitle(f'AMBER Relaxation Effect on {label} {direction}\n'
                 f'Green pipeline — points below diagonal = improvement',
                 fontsize=13, fontweight='bold', y=1.06)

    save_fig(fig, f'scatter_{metric}')


def make_per_target_bars(data, metric, label, lower_better):
    """
    Per-target bar chart for one metric.

    3 rows:
      1. Crystal (reference, single bar per target)
      2. AF: initial (af_unrelaxed) as black bars, AF relaxed + AMBER(AF) as colored bars
      3. Boltz: initial (boltz) as black bars, AMBER(Boltz) as colored bars

    Targets sorted by initial clashscore/metric value (descending for lower_better).
    """
    rows = [
        {
            'label': 'Crystal',
            'initial_src': 'crystal',
            'after_sources': [],
            'after_colors': [],
            'after_labels': [],
        },
        {
            'label': 'AlphaFold',
            'initial_src': 'af_unrelaxed',
            'after_sources': ['af_relaxed', 'amber_af'],
            'after_colors': ['#0072B2', '#E69F00'],
            'after_labels': ['AF relaxed (built-in)', 'AMBER(AF) (standalone)'],
        },
        {
            'label': 'Boltz-1',
            'initial_src': 'boltz',
            'after_sources': ['amber_boltz'],
            'after_colors': ['#D55E00'],
            'after_labels': ['AMBER(Boltz)'],
        },
    ]

    # Get sorted target order (by AF unrelaxed value, descending for lower_better)
    af_init = data[data['source'] == 'af_unrelaxed'].groupby('target')[metric].mean()
    if lower_better:
        sorted_targets = af_init.sort_values(ascending=False).index.tolist()
    else:
        sorted_targets = af_init.sort_values(ascending=True).index.tolist()
    # Add any targets not in AF data
    all_targets = data['target'].unique()
    for t in all_targets:
        if t not in sorted_targets:
            sorted_targets.append(t)

    n_targets = len(sorted_targets)
    x = np.arange(n_targets)

    fig, axes = plt.subplots(3, 1, figsize=(max(20, n_targets * 0.1), 12), sharex=True)
    fig.subplots_adjust(hspace=0.08)

    for row_idx, row_cfg in enumerate(rows):
        ax = axes[row_idx]

        # Get initial values
        init_data = data[data['source'] == row_cfg['initial_src']].groupby('target')[metric].mean()
        init_vals = [init_data.get(t, np.nan) for t in sorted_targets]

        # Number of after bars
        n_after = len(row_cfg['after_sources'])

        if n_after == 0:
            # Crystal: just show bars
            ax.bar(x, init_vals, width=0.8, color='#e74c3c', alpha=0.7,
                   edgecolor='white', linewidth=0.2, label='Crystal')
        else:
            width = 0.8 / max(n_after, 1)

            # After values
            for i, (after_src, after_color, after_label) in enumerate(
                    zip(row_cfg['after_sources'], row_cfg['after_colors'], row_cfg['after_labels'])):
                after_data = data[data['source'] == after_src].groupby('target')[metric].mean()
                after_vals = [after_data.get(t, np.nan) for t in sorted_targets]

                offset = (i - (n_after - 1) / 2) * width
                ax.bar(x + offset, after_vals, width, color=after_color,
                       edgecolor='white', linewidth=0.2, label=after_label, zorder=3)

            # Initial values as wide semi-transparent bars (overlay)
            bar_group_width = width * n_after * 1.2
            ax.bar(x, init_vals, bar_group_width, alpha=0.3, color='black',
                   label=f'Initial ({row_cfg["initial_src"].replace("_", " ")})', zorder=10)

        ax.set_ylabel(row_cfg['label'], fontsize=10, fontweight='bold')

        # y-axis limit
        all_row_vals = [v for v in init_vals if not np.isnan(v)]
        if all_row_vals:
            ymax = np.percentile(all_row_vals, 98) * 1.3
            ax.set_ylim(0, max(ymax, 0.1))

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        if row_idx == 0:
            ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
        else:
            ax.legend(loc='upper right', fontsize=7, framealpha=0.9, ncol=n_after + 1)

    # X-axis labels — show every Nth target
    step = max(1, n_targets // 40)
    tick_positions = list(range(0, n_targets, step))
    axes[-1].set_xticks(tick_positions)
    axes[-1].set_xticklabels([sorted_targets[i] for i in tick_positions],
                              rotation=45, ha='right', fontsize=7)

    direction = '(lower = better)' if lower_better else '(higher = better)'
    fig.suptitle(f'{label} Per Target {direction}\n'
                 f'Green pipeline — {n_targets} targets sorted by initial value — '
                 f'black bars = before AMBER, colored = after AMBER',
                 fontsize=13, fontweight='bold')

    # Common y-axis label
    fig.text(0.01, 0.5, label, va='center', rotation='vertical', fontsize=12)

    save_fig(fig, f'bars_{metric}')


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("Loading Green pipeline MolProbity data...")
    data = load_data()
    print(f"  {len(data)} rows, {data['target'].nunique()} targets, "
          f"{data['source'].nunique()} sources")

    for metric, label, lower_better in METRICS:
        n_valid = data[metric].notna().sum()
        if n_valid == 0:
            print(f"\n  SKIP {metric}: no valid data")
            continue

        print(f"\n--- {label} ({metric}) ---")
        make_scatter(data, metric, label, lower_better)
        make_per_target_bars(data, metric, label, lower_better)

    print(f"\n{'='*50}")
    print(f"All Green presentation figures saved to: {OUTDIR}")
    print(f"Total figures: {len([f for f in os.listdir(OUTDIR) if f.endswith('.png')])}")


if __name__ == '__main__':
    main()
