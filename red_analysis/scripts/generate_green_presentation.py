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


PROTOCOLS = ['cartesian_beta', 'cartesian_ref15', 'dualspace_beta',
             'dualspace_ref15', 'normal_beta', 'normal_ref15']

PROTO_COLORS = {
    'cartesian_beta':    '#e74c3c',   # Red
    'cartesian_ref15':   '#f39c12',   # Orange
    'dualspace_beta':    '#2ecc71',   # Green
    'dualspace_ref15':   '#1abc9c',   # Teal
    'normal_beta':       '#3498db',   # Blue
    'normal_ref15':      '#9b59b6',   # Purple
}

PROTO_LABELS = {
    'cartesian_beta':    'Cart β',
    'cartesian_ref15':   'Cart REF15',
    'dualspace_beta':    'Dual β',
    'dualspace_ref15':   'Dual REF15',
    'normal_beta':       'Norm β',
    'normal_ref15':      'Norm REF15',
}

# Row configs: one row per source type
BAR_ROWS = [
    {
        'label': 'Crystal',
        'rosetta_src': 'crystal',       # Rosetta input source
        'initial_src': 'crystal',       # Pre-Rosetta initial (grey bar)
        'ref_src': None,                # No AMBER reference
        'amber_color': None,
    },
    {
        'label': 'AF unrelaxed',
        'rosetta_src': 'af_unrelaxed',
        'initial_src': 'af_unrelaxed',
        'ref_src': None,
        'amber_color': None,
    },
    {
        'label': 'AF relaxed\n(built-in AMBER)',
        'rosetta_src': 'af_relaxed',
        'initial_src': 'af_relaxed',
        'ref_src': 'af_unrelaxed',      # Grey reference = before built-in AMBER
        'amber_color': '#0072B2',       # Blue for built-in AMBER
    },
    {
        'label': 'AMBER(AF)\n(standalone)',
        'rosetta_src': 'amber_af',
        'initial_src': 'amber_af',
        'ref_src': 'af_unrelaxed',
        'amber_color': '#E69F00',       # Amber/gold
    },
    {
        'label': 'Boltz',
        'rosetta_src': 'boltz',
        'initial_src': 'boltz',
        'ref_src': None,
        'amber_color': None,
    },
    {
        'label': 'AMBER(Boltz)\n(standalone)',
        'rosetta_src': 'amber_boltz',
        'initial_src': 'amber_boltz',
        'ref_src': 'boltz',
        'amber_color': '#D55E00',       # Red-orange
    },
]


def compute_per_protocol_stats(ros, source, metric):
    """
    For each target × protocol, compute mean/min/max across reps.
    Returns DataFrame with columns: target, protocol, mean, min, max.
    """
    src_data = ros[ros['source'] == source]
    if len(src_data) == 0:
        return pd.DataFrame(columns=['target', 'protocol', 'mean', 'min', 'max'])
    stats = src_data.groupby(['target', 'protocol'])[metric].agg(
        ['mean', 'min', 'max']).reset_index()
    return stats


def make_per_target_bars(mp, ros, metric, label, lower_better):
    """
    Per-target bar chart for one metric.

    6 rows (one per source type — AMBER and non-AMBER kept separate):
      1. Crystal
      2. AF unrelaxed
      3. AF relaxed (built-in AMBER)
      4. AMBER(AF) (standalone)
      5. Boltz
      6. AMBER(Boltz) (standalone)

    Per protein in each row:
      - Grey bar: initial value (for AMBER rows: the pre-AMBER reference)
      - Colored bar (AMBER rows only): post-AMBER initial value
      - 6 colored bars: one per Rosetta protocol (mean across reps)
      - Min–max error bars on each Rosetta bar
      - 1 bar gap between proteins

    Targets sorted by initial value (af_unrelaxed).
    """
    n_rows = len(BAR_ROWS)

    # Sorted target order (by af_unrelaxed initial)
    af_init = mp[mp['source'] == 'af_unrelaxed'].groupby('target')[metric].mean()
    if lower_better:
        sorted_targets = af_init.sort_values(ascending=False).index.tolist()
    else:
        sorted_targets = af_init.sort_values(ascending=True).index.tolist()
    for t in mp['target'].unique():
        if t not in sorted_targets:
            sorted_targets.append(t)

    n_targets = len(sorted_targets)

    # Group layout: bars per protein group
    # Non-AMBER: 1 (grey initial) + 6 (protocols) = 7 bars + 1 gap = 8 positions
    # AMBER:     1 (grey ref) + 1 (AMBER initial) + 6 (protocols) = 8 bars + 1 gap = 9 positions
    # Use 9 for all to keep protein positions aligned across rows
    GROUP_W = 9

    fig_width = max(30, n_targets * GROUP_W * 0.07)
    fig, axes = plt.subplots(n_rows, 1, figsize=(fig_width, n_rows * 2.8), sharex=True)
    fig.subplots_adjust(hspace=0.12)

    for row_idx, row_cfg in enumerate(BAR_ROWS):
        ax = axes[row_idx]
        has_amber = row_cfg['ref_src'] is not None

        # Pre-compute Rosetta stats for this source
        proto_stats = compute_per_protocol_stats(ros, row_cfg['rosetta_src'], metric)

        # Pre-compute initial values
        init_data = mp[mp['source'] == row_cfg['initial_src']].groupby('target')[metric].mean()
        ref_data = None
        if has_amber:
            ref_data = mp[mp['source'] == row_cfg['ref_src']].groupby('target')[metric].mean()

        for t_idx, target in enumerate(sorted_targets):
            base_x = t_idx * GROUP_W
            bar_pos = 0

            # Bar 0: grey reference (pre-AMBER ref for AMBER rows, or initial for non-AMBER)
            if has_amber and ref_data is not None:
                val = ref_data.get(target, np.nan)
                if not np.isnan(val):
                    ax.bar(base_x + bar_pos, val, width=0.8, color='#888888',
                           edgecolor='white', linewidth=0.2, zorder=2)
                bar_pos += 1

                # Bar 1: colored AMBER initial
                val = init_data.get(target, np.nan)
                if not np.isnan(val):
                    ax.bar(base_x + bar_pos, val, width=0.8,
                           color=row_cfg['amber_color'],
                           edgecolor='white', linewidth=0.2, zorder=2)
                bar_pos += 1
            else:
                # Non-AMBER: grey initial bar
                val = init_data.get(target, np.nan)
                if not np.isnan(val):
                    ax.bar(base_x + bar_pos, val, width=0.8, color='#888888',
                           edgecolor='white', linewidth=0.2, zorder=2)
                bar_pos += 1
                bar_pos += 1  # Skip 1 position to align protocols at positions 2-7

            # Bars 2-7: 6 Rosetta protocols
            target_proto = proto_stats[proto_stats['target'] == target]
            for proto in PROTOCOLS:
                p_row = target_proto[target_proto['protocol'] == proto]
                if len(p_row) == 1:
                    y_mean = p_row['mean'].values[0]
                    y_min = p_row['min'].values[0]
                    y_max = p_row['max'].values[0]
                    yerr_lo = max(y_mean - y_min, 0)
                    yerr_hi = max(y_max - y_mean, 0)

                    ax.bar(base_x + bar_pos, y_mean, width=0.8,
                           color=PROTO_COLORS[proto],
                           edgecolor='white', linewidth=0.2, zorder=2)
                    ax.errorbar(base_x + bar_pos, y_mean,
                                yerr=[[yerr_lo], [yerr_hi]],
                                fmt='none', ecolor='black', capsize=1.5,
                                elinewidth=0.5, zorder=3)
                bar_pos += 1

            # Position 8 = gap (no bar drawn)

        # Y-axis limits
        all_vals = []
        for t in sorted_targets:
            v = init_data.get(t, np.nan)
            if not np.isnan(v):
                all_vals.append(v)
            if has_amber and ref_data is not None:
                v = ref_data.get(t, np.nan)
                if not np.isnan(v):
                    all_vals.append(v)
        if all_vals:
            ymax = np.percentile(all_vals, 98) * 1.4
            ax.set_ylim(0, max(ymax, 0.1))

        ax.set_ylabel(row_cfg['label'], fontsize=9, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Legend (only on first row to avoid clutter)
        if row_idx == 0:
            from matplotlib.patches import Patch
            handles = [Patch(facecolor='#888888', label='Initial (pre-Rosetta)')]
            for proto in PROTOCOLS:
                handles.append(Patch(facecolor=PROTO_COLORS[proto],
                                     label=PROTO_LABELS[proto]))
            ax.legend(handles=handles, loc='upper right', fontsize=6,
                      framealpha=0.9, ncol=4)
        elif row_idx == 2:
            # Add AMBER legend on first AMBER row
            from matplotlib.patches import Patch
            handles = [
                Patch(facecolor='#888888', label='Pre-AMBER'),
                Patch(facecolor=row_cfg['amber_color'], label='Post-AMBER'),
            ]
            ax.legend(handles=handles, loc='upper right', fontsize=6,
                      framealpha=0.9, ncol=2)

    # X-axis labels: protein names at center of each group
    tick_step = max(1, n_targets // 50)
    tick_positions = [i * GROUP_W + GROUP_W // 2 for i in range(0, n_targets, tick_step)]
    tick_labels = [sorted_targets[i] for i in range(0, n_targets, tick_step)]
    axes[-1].set_xticks(tick_positions)
    axes[-1].set_xticklabels(tick_labels, rotation=45, ha='right', fontsize=6)
    axes[-1].set_xlim(-1, n_targets * GROUP_W)

    direction = '(lower = better)' if lower_better else '(higher = better)'
    fig.suptitle(f'{label} Per Target {direction}\n'
                 f'Green pipeline — {n_targets} targets — '
                 f'grey = initial, colors = 6 Rosetta protocols (error bars = rep range)',
                 fontsize=12, fontweight='bold')

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
