#!/usr/bin/env python3
"""
generate_green_presentation.py — Red Analysis Pipeline

Generates all presentation figures for the Green pipeline.
For every MolProbity metric, produces:
  1. 4-panel scatter: initial vs Rosetta-averaged final, 2 colors = AMBER vs non-AMBER
     Error bars = range (min–max) across 6 Rosetta protocols
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
    """Load pre-Rosetta MolProbity and Rosetta MolProbity, Green pipeline only."""
    metric_cols = [m[0] for m in METRICS]

    # Pre-Rosetta MolProbity (initial values)
    mp = pd.read_csv(os.path.join(METRICS_DIR, "combined_molprobity.tsv"), sep='\t')
    for col in metric_cols:
        mp[col] = pd.to_numeric(mp[col], errors='coerce')
    mp = mp[mp['pipeline'] == 'green']

    # Rosetta MolProbity (post-Rosetta values)
    ros_path = os.path.join(METRICS_DIR, "combined_rosetta_molprobity.tsv")
    ros = pd.read_csv(ros_path, sep='\t')
    for col in metric_cols:
        if col in ros.columns:
            ros[col] = pd.to_numeric(ros[col], errors='coerce')
    ros = ros[ros['pipeline'] == 'green']

    return mp, ros


def compute_rosetta_stats(ros, source, metric):
    """
    For each target, compute Rosetta mean, min, max across 6 protocols.
    Returns DataFrame with columns: target, mean, min, max
    """
    src_data = ros[ros['source'] == source]
    if len(src_data) == 0:
        return pd.DataFrame(columns=['target', 'mean', 'min', 'max'])

    # Per target × protocol mean (averaging across reps first)
    per_proto = src_data.groupby(['target', 'protocol'])[metric].mean().reset_index()
    # Then per target: mean, min, max across protocols
    stats = per_proto.groupby('target')[metric].agg(['mean', 'min', 'max']).reset_index()
    return stats


def make_scatter(mp, ros, metric, label, lower_better):
    """
    4-panel scatter: initial (pre-Rosetta) vs final (Rosetta-averaged).

    Panel 1: Crystal — single color reference
    Panel 2: AF — non-AMBER (af_unrelaxed→Rosetta) vs AMBER (amber_af→Rosetta)
    Panel 3: AF built-in — non-AMBER (af_unrelaxed→Rosetta) vs built-in AMBER (af_relaxed→Rosetta)
    Panel 4: Boltz — non-AMBER (boltz→Rosetta) vs AMBER (amber_boltz→Rosetta)

    Error bars = range (min–max across 6 Rosetta protocols).
    """
    panels = [
        {
            'title': 'Crystal\n(single reference)',
            'sets': [
                ('crystal', 'crystal', '#e74c3c', 'Crystal'),
            ],
        },
        {
            'title': 'AlphaFold\n(standalone AMBER)',
            'sets': [
                ('af_unrelaxed', 'af_unrelaxed', '#56B4E9', 'AF unrelaxed → Rosetta'),
                ('amber_af', 'amber_af', '#E69F00', 'AMBER(AF) → Rosetta'),
            ],
        },
        {
            'title': 'AlphaFold\n(built-in AMBER)',
            'sets': [
                ('af_unrelaxed', 'af_unrelaxed', '#56B4E9', 'AF unrelaxed → Rosetta'),
                ('af_relaxed', 'af_relaxed', '#0072B2', 'AF relaxed → Rosetta'),
            ],
        },
        {
            'title': 'Boltz-1\n(standalone AMBER)',
            'sets': [
                ('boltz', 'boltz', '#009E73', 'Boltz → Rosetta'),
                ('amber_boltz', 'amber_boltz', '#D55E00', 'AMBER(Boltz) → Rosetta'),
            ],
        },
    ]

    fig, axes = plt.subplots(1, 4, figsize=(18, 5))
    fig.subplots_adjust(wspace=0.3)

    # Collect all x and y values to set sensible per-axis limits
    all_x_vals = []
    all_y_vals = []

    for idx, panel in enumerate(panels):
        ax = axes[idx]
        panel_x = []
        panel_y = []

        for initial_src, rosetta_src, color, set_label in panel['sets']:
            # Initial values (pre-Rosetta) — from the source itself
            initial = mp[mp['source'] == initial_src].groupby('target')[metric].mean()

            # Rosetta stats (mean, min, max across protocols)
            ros_stats = compute_rosetta_stats(ros, rosetta_src, metric)
            if len(ros_stats) == 0:
                continue
            ros_stats = ros_stats.set_index('target')

            common = initial.index.intersection(ros_stats.index)
            if len(common) < 3:
                continue

            x = initial.loc[common].values
            y_mean = ros_stats.loc[common, 'mean'].values
            y_min = ros_stats.loc[common, 'min'].values
            y_max = ros_stats.loc[common, 'max'].values

            mask = np.isfinite(x) & np.isfinite(y_mean)
            x = x[mask]
            y_mean = y_mean[mask]
            y_min = y_min[mask]
            y_max = y_max[mask]

            # Error bars: range (min to max), clamp negatives from float precision
            yerr_lo = np.maximum(y_mean - y_min, 0)
            yerr_hi = np.maximum(y_max - y_mean, 0)

            ax.errorbar(x, y_mean, yerr=[yerr_lo, yerr_hi],
                        fmt='o', color=color, alpha=0.6, capsize=2,
                        markersize=5, elinewidth=0.8, markeredgecolor='white',
                        markeredgewidth=0.3, label=set_label, zorder=3)

            panel_x.extend(x)
            panel_y.extend(y_max)  # use max for y-limit
            all_x_vals.extend(x)
            all_y_vals.extend(y_max)

        # X fixed at 90, Y fixed at 10
        ax.set_xlim(0, 90)
        ax.set_ylim(0, 10)
        ax.plot([0, 90], [0, 90], 'k--', alpha=0.4, lw=1)
        ax.set_title(panel['title'], fontsize=9, fontweight='bold')
        ax.legend(fontsize=6.5, loc='upper left', framealpha=0.9)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel(f'Initial {label}', fontsize=9)

        if idx == 0:
            ax.set_ylabel(f'Final {label}\n(Rosetta avg across protocols)', fontsize=10)

    fig.text(0.5, -0.02, f'Initial {label} (pre-Rosetta)', ha='center', fontsize=11)
    direction = '(lower = better)' if lower_better else '(higher = better)'
    fig.suptitle(f'{label} {direction}: Initial vs Rosetta-Averaged Final\n'
                 f'Green pipeline — error bars = range across 6 Rosetta protocols',
                 fontsize=13, fontweight='bold', y=1.06)

    save_fig(fig, f'scatter_{metric}')


def make_per_target_bars(mp, ros, metric, label, lower_better):
    """
    Per-target bar chart for one metric.

    3 rows:
      1. Crystal (reference, single bar per target)
      2. AF: initial (af_unrelaxed) as black bars, AF relaxed + AMBER(AF) as colored bars
      3. Boltz: initial (boltz) as black bars, AMBER(Boltz) as colored bars

    Targets sorted by initial value.
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

    # Sorted target order
    af_init = mp[mp['source'] == 'af_unrelaxed'].groupby('target')[metric].mean()
    if lower_better:
        sorted_targets = af_init.sort_values(ascending=False).index.tolist()
    else:
        sorted_targets = af_init.sort_values(ascending=True).index.tolist()
    all_targets = mp['target'].unique()
    for t in all_targets:
        if t not in sorted_targets:
            sorted_targets.append(t)

    n_targets = len(sorted_targets)
    x = np.arange(n_targets)

    fig, axes = plt.subplots(3, 1, figsize=(max(20, n_targets * 0.1), 12), sharex=True)
    fig.subplots_adjust(hspace=0.08)

    for row_idx, row_cfg in enumerate(rows):
        ax = axes[row_idx]

        init_data = mp[mp['source'] == row_cfg['initial_src']].groupby('target')[metric].mean()
        init_vals = [init_data.get(t, np.nan) for t in sorted_targets]

        n_after = len(row_cfg['after_sources'])

        if n_after == 0:
            ax.bar(x, init_vals, width=0.8, color='#e74c3c', alpha=0.7,
                   edgecolor='white', linewidth=0.2, label='Crystal')
        else:
            width = 0.8 / max(n_after, 1)

            for i, (after_src, after_color, after_label) in enumerate(
                    zip(row_cfg['after_sources'], row_cfg['after_colors'], row_cfg['after_labels'])):
                after_data = mp[mp['source'] == after_src].groupby('target')[metric].mean()
                after_vals = [after_data.get(t, np.nan) for t in sorted_targets]

                offset = (i - (n_after - 1) / 2) * width
                ax.bar(x + offset, after_vals, width, color=after_color,
                       edgecolor='white', linewidth=0.2, label=after_label, zorder=3)

            bar_group_width = width * n_after * 1.2
            ax.bar(x, init_vals, bar_group_width, alpha=0.3, color='black',
                   label=f'Initial ({row_cfg["initial_src"].replace("_", " ")})', zorder=10)

        ax.set_ylabel(row_cfg['label'], fontsize=10, fontweight='bold')

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

    fig.text(0.01, 0.5, label, va='center', rotation='vertical', fontsize=12)

    save_fig(fig, f'bars_{metric}')


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("Loading Green pipeline data...")
    mp, ros = load_data()
    print(f"  MolProbity: {len(mp)} rows, {mp['target'].nunique()} targets")
    print(f"  Rosetta MP: {len(ros)} rows, {ros['target'].nunique()} targets")
    print(f"  Rosetta protocols: {sorted(ros['protocol'].unique())}")

    for metric, label, lower_better in METRICS:
        n_valid_mp = mp[metric].notna().sum()
        n_valid_ros = ros[metric].notna().sum() if metric in ros.columns else 0
        if n_valid_mp == 0:
            print(f"\n  SKIP {metric}: no valid data")
            continue

        print(f"\n--- {label} ({metric}) ---")
        make_scatter(mp, ros, metric, label, lower_better)
        make_per_target_bars(mp, ros, metric, label, lower_better)

    print(f"\n{'='*50}")
    n_figs = len([f for f in os.listdir(OUTDIR) if f.endswith('.png')])
    print(f"All Green presentation figures saved to: {OUTDIR}")
    print(f"Total figures: {n_figs}")


if __name__ == '__main__':
    main()
