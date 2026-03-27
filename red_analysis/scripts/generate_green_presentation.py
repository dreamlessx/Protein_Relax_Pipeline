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

# Per-metric scatter axis bounds: (x_min, x_max, y_min, y_max)
# Tight bounds — x captures initial range, y zoomed to Rosetta range
# INCLUDING error bars (min-max range across 6 protocols).
# Based on data: p99/max analysis of Green pipeline (257 targets).
SCATTER_BOUNDS = {
    'clashscore':       (0,     60,     0,     8),     # x p99=50; errbar top max=7.6
    'molprobity_score': (0,     3.5,    0,     1.2),   # x p99=3.15; errbar top max=1.17
    'rama_outliers':    (-0.3,  6,      -0.3,  3.5),   # x p99=5.0; errbar top max=2.96; pad below 0
    'rama_favored':     (82,    100,    90,    100),   # zoom into top range; errbar bot=90.6
    'rota_outliers':    (-0.5,  18,     -0.05,  0.5),   # x p99=16.7; errbar top max=0.47; pad below 0
    'rms_bonds':        (0.005, 0.04,   0.008, 0.04),  # x p99=0.031; errbar top max=0.074 (clipped)
    'rms_angles':       (0.8,   3.0,    1.0,   3.0),  # x p99=2.6; errbar top p99=2.51
    'cbeta_outliers':   (0,     12,     0,     9),     # x p99=10.6; errbar top max=8.2
}


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

        x_min, x_max, y_min, y_max = SCATTER_BOUNDS.get(metric, (0, 90, 0, 10))
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        diag_min = min(x_min, y_min)
        diag_max = max(x_max, y_max)
        ax.plot([diag_min, diag_max], [diag_min, diag_max], 'k--', alpha=0.4, lw=1)
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

    6 rows (one per source type — AMBER and non-AMBER separate):
      1. Crystal, 2. AF unrelaxed, 3. AF relaxed (built-in AMBER),
      4. AMBER(AF) (standalone), 5. Boltz, 6. AMBER(Boltz) (standalone)

    Per protein in each row:
      - Wide transparent grey bar spanning the 6 protocol positions = initial value
      - 6 colored bars on top: one per Rosetta protocol (mean across reps)
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

    # Group layout: 6 protocol bars + 1 gap = 7 positions per protein
    GROUP_W = 7

    fig_width = max(30, n_targets * GROUP_W * 0.06)
    fig, axes = plt.subplots(n_rows, 1, figsize=(fig_width, n_rows * 2.5), sharex=True)
    fig.subplots_adjust(hspace=0.12)

    for row_idx, row_cfg in enumerate(BAR_ROWS):
        ax = axes[row_idx]
        source = row_cfg['source']

        # Pre-compute Rosetta stats for this source
        proto_stats = compute_per_protocol_stats(ros, source, metric)

        # Pre-compute initial values (pre-Rosetta)
        init_data = mp[mp['source'] == source].groupby('target')[metric].mean()

        for t_idx, target in enumerate(sorted_targets):
            base_x = t_idx * GROUP_W

            # Wide transparent grey bar = initial (pre-Rosetta) value
            # Spans the 6 protocol positions (center at base_x + 2.5, width = 6)
            init_val = init_data.get(target, np.nan)
            if not np.isnan(init_val):
                ax.bar(base_x + 2.5, init_val, width=5.6, color='#888888',
                       alpha=0.25, edgecolor='none', zorder=1)

            # 6 Rosetta protocol bars (positions 0–5)
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

            # Position 6 = gap (no bar)

        # Y-axis limits: based on initial values (which are typically larger)
        all_vals = [init_data.get(t, np.nan) for t in sorted_targets]
        all_vals = [v for v in all_vals if not np.isnan(v)]
        if all_vals:
            ymax = np.percentile(all_vals, 98) * 1.4
            ax.set_ylim(0, max(ymax, 0.1))

        ax.set_ylabel(row_cfg['label'], fontsize=9, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Legend on first row
        if row_idx == 0:
            from matplotlib.patches import Patch
            handles = [Patch(facecolor='#888888', alpha=0.25, label='Initial (pre-Rosetta)')]
            for proto in PROTOCOLS:
                handles.append(Patch(facecolor=PROTO_COLORS[proto],
                                     label=PROTO_LABELS[proto]))
            ax.legend(handles=handles, loc='upper right', fontsize=6,
                      framealpha=0.9, ncol=4)

    # X-axis labels: every protein labeled
    tick_positions = [i * GROUP_W + 2.5 for i in range(n_targets)]
    axes[-1].set_xticks(tick_positions)
    axes[-1].set_xticklabels(sorted_targets, rotation=90, ha='center', fontsize=4)
    axes[-1].set_xlim(-1, n_targets * GROUP_W)

    direction = '(lower = better)' if lower_better else '(higher = better)'
    fig.suptitle(f'{label} Per Target {direction}\n'
                 f'Green pipeline — {n_targets} targets — '
                 f'grey = initial, colors = 6 Rosetta protocols (error bars = rep range)',
                 fontsize=12, fontweight='bold')

    save_fig(fig, f'bars_{metric}')


def make_energy_bars(energy_df):
    """
    Per-target bar chart for Rosetta per-residue energy (REU/residue).

    Same 6-row layout as MolProbity bars:
      1 row per source, 6 protocol bars per target, sorted by mean energy.
    No grey initial bar (energy is only post-Rosetta).
    """
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

    fig_width = max(30, n_targets * GROUP_W * 0.06)
    fig, axes = plt.subplots(n_rows, 1, figsize=(fig_width, n_rows * 2.5), sharex=True)
    fig.subplots_adjust(hspace=0.12)

    for row_idx, row_cfg in enumerate(BAR_ROWS):
        ax = axes[row_idx]
        source = row_cfg['source']
        src_data = energy_df[energy_df['source'] == source]

        # Per target × protocol stats
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

        # Y limits
        all_vals = src_data.groupby('target')[metric].mean().values
        if len(all_vals) > 0:
            ymin = np.nanpercentile(all_vals, 2) * 1.15
            ymax = np.nanpercentile(all_vals, 98) * 0.85
            ax.set_ylim(ymin, ymax)

        ax.set_ylabel(row_cfg['label'], fontsize=9, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        if row_idx == 0:
            from matplotlib.patches import Patch
            handles = []
            for proto in PROTOCOLS:
                handles.append(Patch(facecolor=PROTO_COLORS[proto],
                                     label=PROTO_LABELS[proto]))
            ax.legend(handles=handles, loc='lower right', fontsize=6,
                      framealpha=0.9, ncol=3)

    tick_positions = [i * GROUP_W + 2.5 for i in range(n_targets)]
    axes[-1].set_xticks(tick_positions)
    axes[-1].set_xticklabels(sorted_targets, rotation=90, ha='center', fontsize=4)
    axes[-1].set_xlim(-1, n_targets * GROUP_W)

    fig.suptitle('Per-Residue Rosetta Energy (REU/residue) Per Target\n'
                 f'Green pipeline — {n_targets} targets — '
                 f'6 Rosetta protocols (error bars = rep range)',
                 fontsize=12, fontweight='bold')

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

    print(f"\n{'='*50}")
    n_figs = len([f for f in os.listdir(OUTDIR) if f.endswith('.png')])
    print(f"All Green presentation figures saved to: {OUTDIR}")
    print(f"Total figures: {n_figs}")


if __name__ == '__main__':
    main()
