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

# Per-metric axis overrides for scatter plots.
# Keys: metric name → dict with optional x_min, x_max, y_min, y_max.
SCATTER_AXIS_OVERRIDES = {
    'rota_outliers': {'x_min': 0.0, 'x_max': 30.0, 'y_min': 0.0, 'y_max': 0.5},
    'rama_outliers': {'x_min': 0.0, 'x_max': 20.0, 'y_min': 0.0, 'y_max': 3.0},
}

# Scatter axis bounds are computed from data in make_scatter().
# x captures full initial range (p1–p99 with 5% pad).
# y captures full Rosetta range including error-bar extremes
# (global min of y_min to p99 of y_max, with 5% pad on each side).
# All 4 panels share the same y-axis bounds.


def save_fig(fig, name):
    # Route to scatter/ or bars/ subdirectory based on prefix
    if name.startswith('scatter_'):
        subdir = os.path.join(OUTDIR, 'scatter')
    elif name.startswith('bars_'):
        subdir = os.path.join(OUTDIR, 'bars')
    else:
        subdir = OUTDIR
    os.makedirs(subdir, exist_ok=True)
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(subdir, f"{name}.{fmt}"),
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

    # First pass: collect data per panel for data-driven axis bounds
    panel_data = []  # list of list of (x, y_mean, y_min, y_max, yerr_lo, yerr_hi, color, label)
    all_x_vals = []
    all_y_min_vals = []
    all_y_max_vals = []

    for idx, panel in enumerate(panels):
        sets_data = []
        for initial_src, rosetta_src, color, set_label in panel['sets']:
            initial = mp[mp['source'] == initial_src].groupby('target')[metric].mean()
            ros_stats = compute_rosetta_stats(ros, rosetta_src, metric)
            if len(ros_stats) == 0:
                sets_data.append(None)
                continue
            ros_stats = ros_stats.set_index('target')
            common = initial.index.intersection(ros_stats.index)
            if len(common) < 3:
                sets_data.append(None)
                continue

            x = initial.loc[common].values
            y_mean = ros_stats.loc[common, 'mean'].values
            y_min_v = ros_stats.loc[common, 'min'].values
            y_max_v = ros_stats.loc[common, 'max'].values

            mask = np.isfinite(x) & np.isfinite(y_mean)
            x, y_mean = x[mask], y_mean[mask]
            y_min_v, y_max_v = y_min_v[mask], y_max_v[mask]

            yerr_lo = np.maximum(y_mean - y_min_v, 0)
            yerr_hi = np.maximum(y_max_v - y_mean, 0)

            sets_data.append((x, y_mean, y_min_v, y_max_v, yerr_lo, yerr_hi, color, set_label))
            all_x_vals.extend(x)
            all_y_min_vals.extend(y_min_v)
            all_y_max_vals.extend(y_max_v)
        panel_data.append(sets_data)

    # Compute data-driven bounds (may be overridden per metric)
    overrides = SCATTER_AXIS_OVERRIDES.get(metric, {})

    if 'x_min' in overrides and 'x_max' in overrides:
        x_min, x_max = overrides['x_min'], overrides['x_max']
    elif all_x_vals:
        x_lo = np.nanpercentile(all_x_vals, 1)
        x_hi = np.nanpercentile(all_x_vals, 99)
        x_pad = max((x_hi - x_lo) * 0.05, abs(x_hi) * 0.01)
        x_min = x_lo - x_pad
        x_max = x_hi + x_pad
    else:
        x_min, x_max = 0, 1

    if 'y_min' in overrides and 'y_max' in overrides:
        y_min, y_max = overrides['y_min'], overrides['y_max']
    elif all_y_min_vals:
        y_lo = np.nanmin(all_y_min_vals)
        y_hi = np.nanpercentile(all_y_max_vals, 99)
        y_range = y_hi - y_lo
        top_pad = max(y_range * 0.05, abs(y_hi) * 0.02)
        bot_pad = y_range * 0.08  # 8% of range below the minimum
        y_min = y_lo - bot_pad
        if y_min >= 0:
            y_min = -bot_pad  # force negative axis
        y_max = y_hi + top_pad
    else:
        y_min, y_max = -0.1, 1

    # Second pass: plot with computed bounds
    for idx, panel in enumerate(panels):
        ax = axes[idx]
        for sd in panel_data[idx]:
            if sd is None:
                continue
            x, y_mean, y_min_v, y_max_v, yerr_lo, yerr_hi, color, set_label = sd
            ax.errorbar(x, y_mean, yerr=[yerr_lo, yerr_hi],
                        fmt='o', color=color, alpha=0.6, capsize=2,
                        markersize=5, elinewidth=0.8, markeredgecolor='white',
                        markeredgewidth=0.3, label=set_label, zorder=3)

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        if idx == 0:
            print(f"    DEBUG {metric}: y_min={y_min:.4f} y_max={y_max:.4f} actual={ax.get_ylim()}")
        diag_lo = min(x_min, y_min)
        diag_hi = max(x_max, y_max)
        ax.plot([diag_lo, diag_hi], [diag_lo, diag_hi], 'k--', alpha=0.4, lw=1)
        ax.set_title(panel['title'], fontsize=9, fontweight='bold')
        ax.legend(fontsize=6.5, loc='upper right', framealpha=0.9)
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
    {'label': 'Crystal',                      'source': 'crystal'},
    {'label': 'AF unrelaxed',                 'source': 'af_unrelaxed'},
    {'label': 'AF relaxed\n(built-in AMBER)', 'source': 'af_relaxed'},
    {'label': 'AMBER(AF)\n(standalone)',       'source': 'amber_af'},
    {'label': 'Boltz',                        'source': 'boltz'},
    {'label': 'AMBER(Boltz)\n(standalone)',    'source': 'amber_boltz'},
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

    6 rows (one per source type), shared y-axis across all rows.
    Grey = initial (pre-Rosetta), colored = 6 Rosetta protocols.
    Error bars = min–max across reps within each protocol.
    Targets sorted by initial af_unrelaxed value.

    Fixes applied:
      - Shared y-axis across all 6 rows for cross-source comparison
      - Smart y-limits: show both initial and Rosetta ranges meaningfully
      - rama_favored: floor at 85% since data clusters at 95–100%
      - Every 5th x-label shown (readable at 300 DPI)
      - Legend outside plot area (top of figure)
    """
    from matplotlib.patches import Patch
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
    GROUP_W = 7

    fig_width = max(30, n_targets * GROUP_W * 0.06)
    fig, axes = plt.subplots(n_rows, 1, figsize=(fig_width, n_rows * 2.8),
                             sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.15, top=0.90)

    # Pre-compute global y-range across ALL rows for shared axis
    all_ros_max = []
    all_init_vals = []
    for row_cfg in BAR_ROWS:
        source = row_cfg['source']
        proto_stats = compute_per_protocol_stats(ros, source, metric)
        init_data = mp[mp['source'] == source].groupby('target')[metric].mean()
        if len(proto_stats) > 0:
            all_ros_max.extend(proto_stats['max'].dropna().values)
        all_init_vals.extend(init_data.dropna().values)

    # Compute shared y-limits
    if metric == 'rama_favored':
        # Higher is better — zoom to 85–100%
        y_lo = 85.0
        y_hi = 100.5
    elif all_ros_max and all_init_vals:
        # Show full range: initial + Rosetta
        combined = list(all_ros_max) + list(all_init_vals)
        y_lo = 0
        y_hi = np.nanpercentile(combined, 99) * 1.1
        y_hi = max(y_hi, 0.1)
    elif all_ros_max:
        y_lo = 0
        y_hi = np.nanpercentile(all_ros_max, 99) * 1.3
    else:
        y_lo, y_hi = 0, 1

    for row_idx, row_cfg in enumerate(BAR_ROWS):
        ax = axes[row_idx]
        source = row_cfg['source']

        proto_stats = compute_per_protocol_stats(ros, source, metric)
        init_data = mp[mp['source'] == source].groupby('target')[metric].mean()

        for t_idx, target in enumerate(sorted_targets):
            base_x = t_idx * GROUP_W

            # Grey initial bar
            init_val = init_data.get(target, np.nan)
            if not np.isnan(init_val):
                ax.bar(base_x + 2.5, init_val, width=5.6, color='#888888',
                       alpha=0.25, edgecolor='none', zorder=1)

            # 6 Rosetta protocol bars
            target_proto = proto_stats[proto_stats['target'] == target]
            for p_idx, proto in enumerate(PROTOCOLS):
                p_row = target_proto[target_proto['protocol'] == proto]
                if len(p_row) == 1:
                    y_mean = p_row['mean'].values[0]
                    y_min = p_row['min'].values[0]
                    y_max = p_row['max'].values[0]
                    yerr_lo = max(y_mean - y_min, 0)
                    yerr_hi = max(y_max - y_mean, 0)

                    ax.bar(base_x + p_idx, y_mean, width=0.8,
                           color=PROTO_COLORS[proto],
                           edgecolor='white', linewidth=0.2, zorder=2)
                    ax.errorbar(base_x + p_idx, y_mean,
                                yerr=[[yerr_lo], [yerr_hi]],
                                fmt='none', ecolor='black', capsize=1.5,
                                elinewidth=0.5, zorder=3)

        ax.set_ylim(y_lo, y_hi)
        ax.set_ylabel(row_cfg['label'], fontsize=9, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Legend above the figure
    handles = [Patch(facecolor='#888888', alpha=0.25, label='Initial (pre-Rosetta)')]
    for proto in PROTOCOLS:
        handles.append(Patch(facecolor=PROTO_COLORS[proto],
                             label=PROTO_LABELS[proto]))
    fig.legend(handles=handles, loc='upper center', fontsize=7,
               framealpha=0.9, ncol=7, bbox_to_anchor=(0.5, 0.98))

    # X-axis labels: every 5th target for readability
    tick_positions = [i * GROUP_W + 2.5 for i in range(n_targets)]
    label_every = 5
    shown_labels = [sorted_targets[i] if i % label_every == 0 else ''
                    for i in range(n_targets)]
    axes[-1].set_xticks(tick_positions)
    axes[-1].set_xticklabels(shown_labels, rotation=90, ha='center', fontsize=5)
    axes[-1].set_xlim(-1, n_targets * GROUP_W)

    direction = '(lower = better)' if lower_better else '(higher = better)'
    fig.suptitle(f'{label} Per Target {direction}\n'
                 f'Green pipeline — {n_targets} targets — '
                 f'grey = initial, colors = 6 Rosetta protocols (error bars = rep range)',
                 fontsize=12, fontweight='bold', y=1.0)

    save_fig(fig, f'bars_{metric}')


def make_energy_bars(energy_df):
    """
    Per-target bar chart for Rosetta per-residue energy (REU/residue).

    Same 6-row layout as MolProbity bars, shared y-axis.
    Energy is negative (more negative = better), so bars hang downward.
    Zero line drawn for reference.
    """
    from matplotlib.patches import Patch
    metric = 'per_residue_energy'
    n_rows = len(BAR_ROWS)

    # Sort targets by af_relaxed mean energy (most negative = best)
    af_energy = energy_df[energy_df['source'] == 'af_relaxed'].groupby('target')[metric].mean()
    sorted_targets = af_energy.sort_values(ascending=True).index.tolist()
    for t in energy_df['target'].unique():
        if t not in sorted_targets:
            sorted_targets.append(t)

    n_targets = len(sorted_targets)
    GROUP_W = 7

    # Pre-compute global y-range for shared axis
    all_vals = energy_df.groupby(['source', 'target'])[metric].mean().values
    y_lo = np.nanpercentile(all_vals, 1) * 1.1  # more negative (expand)
    y_hi = np.nanpercentile(all_vals, 99) * 0.9  # less negative (expand toward 0)
    y_hi = max(y_hi, 0.5)  # always show zero line

    fig_width = max(30, n_targets * GROUP_W * 0.06)
    fig, axes = plt.subplots(n_rows, 1, figsize=(fig_width, n_rows * 2.8),
                             sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.15, top=0.90)

    for row_idx, row_cfg in enumerate(BAR_ROWS):
        ax = axes[row_idx]
        source = row_cfg['source']
        src_data = energy_df[energy_df['source'] == source]

        proto_stats = src_data.groupby(['target', 'protocol'])[metric].agg(
            ['mean', 'min', 'max']).reset_index()

        for t_idx, target in enumerate(sorted_targets):
            base_x = t_idx * GROUP_W
            target_proto = proto_stats[proto_stats['target'] == target]

            for p_idx, proto in enumerate(PROTOCOLS):
                p_row = target_proto[target_proto['protocol'] == proto]
                if len(p_row) == 1:
                    y_mean = p_row['mean'].values[0]
                    y_min = p_row['min'].values[0]
                    y_max = p_row['max'].values[0]
                    yerr_lo = max(y_mean - y_min, 0)
                    yerr_hi = max(y_max - y_mean, 0)

                    ax.bar(base_x + p_idx, y_mean, width=0.8,
                           color=PROTO_COLORS[proto],
                           edgecolor='white', linewidth=0.2, zorder=2)
                    ax.errorbar(base_x + p_idx, y_mean,
                                yerr=[[yerr_lo], [yerr_hi]],
                                fmt='none', ecolor='black', capsize=1.5,
                                elinewidth=0.5, zorder=3)

        ax.set_ylim(y_lo, y_hi)
        ax.axhline(y=0, color='black', linewidth=0.5, linestyle='-', alpha=0.3)
        ax.set_ylabel(row_cfg['label'], fontsize=9, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Legend above figure
    handles = [Patch(facecolor=PROTO_COLORS[p], label=PROTO_LABELS[p])
               for p in PROTOCOLS]
    fig.legend(handles=handles, loc='upper center', fontsize=7,
               framealpha=0.9, ncol=6, bbox_to_anchor=(0.5, 0.98))

    # X-axis labels: every 5th
    tick_positions = [i * GROUP_W + 2.5 for i in range(n_targets)]
    label_every = 5
    shown_labels = [sorted_targets[i] if i % label_every == 0 else ''
                    for i in range(n_targets)]
    axes[-1].set_xticks(tick_positions)
    axes[-1].set_xticklabels(shown_labels, rotation=90, ha='center', fontsize=5)
    axes[-1].set_xlim(-1, n_targets * GROUP_W)

    fig.suptitle('Per-Residue Rosetta Energy (REU/residue) Per Target\n'
                 f'Green pipeline — {n_targets} targets — '
                 f'6 Rosetta protocols (error bars = rep range)',
                 fontsize=12, fontweight='bold', y=1.0)

    save_fig(fig, 'bars_energy')


def make_energy_scatter(energy_df, mp):
    """
    Energy scatter: per-residue energy vs clashscore/MolProbity score.

    Shows whether lower Rosetta energy correlates with better structure quality.
    4-panel: one per source type (crystal, af_relaxed, boltz, amber_boltz).
    """
    # Merge energy and MolProbity (both post-Rosetta, per target × source)
    energy_per_target = energy_df.groupby(['source', 'target'])['per_residue_energy'].mean().reset_index()
    mp_per_target = mp.groupby(['source', 'target'])['clashscore'].mean().reset_index()

    sources = [
        ('af_relaxed', 'AF relaxed → Rosetta', '#0072B2'),
        ('af_unrelaxed', 'AF unrelaxed → Rosetta', '#56B4E9'),
        ('boltz', 'Boltz → Rosetta', '#009E73'),
        ('amber_boltz', 'AMBER(Boltz) → Rosetta', '#D55E00'),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(18, 5))
    fig.subplots_adjust(wspace=0.3)

    for idx, (src, label, color) in enumerate(sources):
        ax = axes[idx]

        e_data = energy_per_target[energy_per_target['source'] == src].set_index('target')
        q_data = mp_per_target[mp_per_target['source'] == src].set_index('target')
        common = e_data.index.intersection(q_data.index)

        if len(common) < 5:
            ax.set_title(label, fontsize=9, fontweight='bold')
            ax.text(0.5, 0.5, 'Insufficient data', transform=ax.transAxes,
                    ha='center', fontsize=10)
            continue

        x = e_data.loc[common, 'per_residue_energy'].values
        y = q_data.loc[common, 'clashscore'].values
        mask = np.isfinite(x) & np.isfinite(y)
        x, y = x[mask], y[mask]

        ax.scatter(x, y, color=color, alpha=0.6, s=25, edgecolors='white',
                   linewidths=0.3, zorder=3)

        # Correlation annotation
        if len(x) > 5:
            from scipy import stats as sp_stats
            r, p = sp_stats.pearsonr(x, y)
            ax.text(0.05, 0.95, f'r = {r:.3f}\nn = {len(x)}',
                    transform=ax.transAxes, fontsize=8, va='top',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                              alpha=0.8, edgecolor='gray', linewidth=0.5))

        ax.set_title(label, fontsize=9, fontweight='bold')
        ax.set_xlabel('Per-Residue Energy (REU/res)', fontsize=9)
        if idx == 0:
            ax.set_ylabel('Clashscore (post-Rosetta)', fontsize=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle('Rosetta Energy vs Structure Quality (Clashscore)\n'
                 'Green pipeline — per-target averages',
                 fontsize=13, fontweight='bold', y=1.06)

    save_fig(fig, 'scatter_energy_vs_clashscore')


def load_initial_energy():
    """Load pre-Rosetta energy scores (Rosetta score-only on input structures)."""
    import glob as gl
    files = sorted(gl.glob(os.path.join(METRICS_DIR, "initial_energy", "initial_energy_*.tsv")))
    if not files:
        return None
    dfs = [pd.read_csv(f, sep='\t') for f in files]
    init = pd.concat(dfs, ignore_index=True)
    for col in ['total_score', 'per_residue_energy']:
        init[col] = pd.to_numeric(init[col], errors='coerce')
    return init


def make_energy_before_after_scatter(init_energy, ros_energy):
    """
    4-panel scatter: initial (pre-Rosetta) vs final (Rosetta-averaged) per-residue energy.

    Same layout as MolProbity scatters:
      Panel 1: Crystal
      Panel 2: AF (standalone AMBER) — af_unrelaxed & amber_af
      Panel 3: AF (built-in AMBER) — af_unrelaxed & af_relaxed
      Panel 4: Boltz (standalone AMBER) — boltz & amber_boltz

    Initial energy only available for crystal and af_unrelaxed.
    For Rosetta sources whose input isn't scored (af_relaxed, amber_af, etc.),
    we use the corresponding un-relaxed initial energy where logical:
      - amber_af, af_relaxed: initial af_unrelaxed energy (same base structure)
      - amber_boltz, boltz: no initial energy available
    """
    metric = 'per_residue_energy'

    # Map: rosetta source → initial source for pairing
    initial_source_map = {
        'crystal': 'crystal',
        'af_unrelaxed': 'af_unrelaxed',
        'af_relaxed': 'af_unrelaxed',
        'amber_af': 'af_unrelaxed',
    }

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
                ('af_unrelaxed', 'amber_af', '#E69F00', 'AMBER(AF) → Rosetta'),
            ],
        },
        {
            'title': 'AlphaFold\n(built-in AMBER)',
            'sets': [
                ('af_unrelaxed', 'af_unrelaxed', '#56B4E9', 'AF unrelaxed → Rosetta'),
                ('af_unrelaxed', 'af_relaxed', '#0072B2', 'AF relaxed → Rosetta'),
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

    for idx, panel in enumerate(panels):
        ax = axes[idx]
        has_data = False

        for init_src, ros_src, color, set_label in panel['sets']:
            if init_src is None:
                continue

            # Initial: per-target mean energy
            init_data = init_energy[init_energy['source'] == init_src]
            init_by_t = init_data.groupby('target')[metric].mean()

            # Rosetta: per-target mean/min/max across protocols
            ros_stats = compute_rosetta_energy_stats(ros_energy, ros_src, metric)
            if len(ros_stats) == 0:
                continue
            ros_stats = ros_stats.set_index('target')

            common = init_by_t.index.intersection(ros_stats.index)
            if len(common) < 3:
                continue

            x = init_by_t.loc[common].values
            y_mean = ros_stats.loc[common, 'mean'].values
            y_min_v = ros_stats.loc[common, 'min'].values
            y_max_v = ros_stats.loc[common, 'max'].values

            mask = np.isfinite(x) & np.isfinite(y_mean)
            x, y_mean = x[mask], y_mean[mask]
            y_min_v, y_max_v = y_min_v[mask], y_max_v[mask]

            yerr_lo = np.maximum(y_mean - y_min_v, 0)
            yerr_hi = np.maximum(y_max_v - y_mean, 0)

            ax.errorbar(x, y_mean, yerr=[yerr_lo, yerr_hi],
                        fmt='o', color=color, alpha=0.6, capsize=2,
                        markersize=5, elinewidth=0.8, markeredgecolor='white',
                        markeredgewidth=0.3, label=set_label, zorder=3)
            has_data = True

        if not has_data:
            ax.text(0.5, 0.5, 'No initial energy data',
                    transform=ax.transAxes, ha='center', va='center',
                    fontsize=10, color='gray')

        # Diagonal (no-change line)
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        diag_lo = min(xlim[0], ylim[0])
        diag_hi = max(xlim[1], ylim[1])
        ax.plot([diag_lo, diag_hi], [diag_lo, diag_hi], 'k--', alpha=0.4, lw=1)

        ax.set_title(panel['title'], fontsize=9, fontweight='bold')
        ax.legend(fontsize=6.5, loc='upper right', framealpha=0.9)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel('Initial Energy (REU/res)', fontsize=9)

        if idx == 0:
            ax.set_ylabel('Final Energy (REU/res)\n(Rosetta avg across protocols)', fontsize=10)

    fig.text(0.5, -0.02, 'Initial Per-Residue Energy (pre-Rosetta)', ha='center', fontsize=11)
    fig.suptitle('Per-Residue Energy (lower = better): Initial vs Rosetta-Averaged Final\n'
                 'Green pipeline — error bars = range across 6 Rosetta protocols',
                 fontsize=13, fontweight='bold', y=1.06)

    save_fig(fig, 'scatter_energy_before_after')


def compute_rosetta_energy_stats(ros_energy, source, metric):
    """Per-target mean/min/max of energy across protocols (averaging reps first)."""
    src_data = ros_energy[ros_energy['source'] == source]
    if len(src_data) == 0:
        return pd.DataFrame(columns=['target', 'mean', 'min', 'max'])
    per_proto = src_data.groupby(['target', 'protocol'])[metric].mean().reset_index()
    stats = per_proto.groupby('target')[metric].agg(['mean', 'min', 'max']).reset_index()
    return stats


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("Loading Green pipeline data...")
    mp, ros = load_data()
    print(f"  MolProbity: {len(mp)} rows, {mp['target'].nunique()} targets")
    print(f"  Rosetta MP: {len(ros)} rows, {ros['target'].nunique()} targets")
    print(f"  Rosetta protocols: {sorted(ros['protocol'].unique())}")

    # Load energy data
    energy_path = os.path.join(METRICS_DIR, "combined_rosetta_energy.tsv")
    energy = None
    if os.path.exists(energy_path):
        energy = pd.read_csv(energy_path, sep='\t')
        energy['per_residue_energy'] = pd.to_numeric(energy['per_residue_energy'], errors='coerce')
        energy['total_score'] = pd.to_numeric(energy['total_score'], errors='coerce')
        energy = energy[energy['pipeline'] == 'green']
        print(f"  Energy: {len(energy)} rows, {energy['target'].nunique()} targets")

    # Load initial energy
    init_energy = load_initial_energy()
    if init_energy is not None:
        print(f"  Initial energy: {len(init_energy)} rows, {init_energy['target'].nunique()} targets")

    for metric, label, lower_better in METRICS:
        n_valid_mp = mp[metric].notna().sum()
        n_valid_ros = ros[metric].notna().sum() if metric in ros.columns else 0
        if n_valid_mp == 0:
            print(f"\n  SKIP {metric}: no valid data")
            continue

        print(f"\n--- {label} ({metric}) ---")
        make_scatter(mp, ros, metric, label, lower_better)
        make_per_target_bars(mp, ros, metric, label, lower_better)

    # Energy figures
    if energy is not None and len(energy) > 0:
        print(f"\n--- Energy Analysis ---")
        make_energy_bars(energy)
        make_energy_scatter(energy, ros)

    # Energy before/after scatter
    if energy is not None and init_energy is not None:
        print(f"\n--- Energy Before/After Scatter ---")
        make_energy_before_after_scatter(init_energy, energy)

    print(f"\n{'='*50}")
    n_figs = len([f for f in os.listdir(OUTDIR) if f.endswith('.png')])
    print(f"All Green presentation figures saved to: {OUTDIR}")
    print(f"Total figures: {n_figs}")


if __name__ == '__main__':
    main()
