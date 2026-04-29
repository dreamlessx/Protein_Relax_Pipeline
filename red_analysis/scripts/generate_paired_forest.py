#!/usr/bin/env python3
"""
generate_paired_forest.py — Red Analysis Pipeline

fig22_amber_paired_forest: Forest plot of the 16-cell main paired MolProbity test
from amber_paired_table2.tsv. Each row is one (pipeline x pair x protocol x metric)
cell with the mean delta and 95% CI. Filled markers = p<0.05 (none expected at the
locked post-Rosetta frame); hollow markers = ns.

Corrected narrative: zero of 16 cells reach p<0.05; the post-Rosetta paired AMBER
effect is in the noise.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

TABLES_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/tables"
FIGDIR = "/data/p_csb_meiler/agarwm5/red_analysis/figures"

PIPELINE_COLORS = {'blue': '#1f6feb', 'green': '#2ea043'}
PROTOCOL_LABELS = {
    'normal_beta':    'Normal beta',
    'dualspace_beta': 'Dualspace beta',
    'cartesian_beta': 'Cartesian beta',
}
PAIR_LABELS = {
    'amber_af_vs_af_relaxed':     'AMBER(AF) vs AF relaxed',
    'amber_boltz_vs_boltz':       'AMBER(Boltz) vs Boltz',
}


def load():
    f = os.path.join(TABLES_DIR, 'amber_paired_table2.tsv')
    df = pd.read_csv(f, sep='\t')
    return df


def fig22_paired_forest(df):
    # Order rows: pipeline (blue first, green second), then pair, then protocol, then metric
    pair_order = ['amber_af_vs_af_relaxed', 'amber_boltz_vs_boltz']
    proto_order = ['normal_beta', 'dualspace_beta']
    metric_order = ['molprobity_score', 'clashscore']

    # 2 pipelines x 2 pairs x 2 protocols x 2 metrics = 16 cells
    df['pipeline'] = pd.Categorical(df['pipeline'], categories=['blue', 'green'])
    df['pair'] = pd.Categorical(df['pair'], categories=pair_order)
    df['protocol'] = pd.Categorical(df['protocol'], categories=proto_order)
    df['metric'] = pd.Categorical(df['metric'], categories=metric_order)
    df = df.sort_values(['pipeline', 'pair', 'protocol', 'metric']).reset_index(drop=True)

    # Make 16 rows for plotting; group by pipeline visually
    n = len(df)
    fig, ax = plt.subplots(figsize=(11, 8.5))

    y_positions = []
    y_labels = []
    cur_y = 0
    last_pipeline = None
    last_pair = None
    for i, row in df.iterrows():
        if last_pipeline is not None and row['pipeline'] != last_pipeline:
            cur_y += 1.0  # gap between pipelines
        if last_pair is not None and (row['pair'] != last_pair or
                                      row['pipeline'] != last_pipeline):
            cur_y += 0.4
        y_positions.append(cur_y)
        # Label: pipeline + pair + protocol + metric
        lbl = (f"{row['pipeline'].capitalize()} | "
               f"{PAIR_LABELS.get(row['pair'], row['pair'])} | "
               f"{PROTOCOL_LABELS.get(row['protocol'], row['protocol'])} | "
               f"{row['metric']}")
        y_labels.append(lbl)
        cur_y += 1.0
        last_pipeline = row['pipeline']
        last_pair = row['pair']

    y_positions = np.array(y_positions)

    for i, row in df.iterrows():
        y = y_positions[i]
        delta = row['mean_diff']
        ci_lo, ci_hi = row['diff_ci_lo'], row['diff_ci_hi']
        p = row['paired_t_p']
        col = PIPELINE_COLORS[row['pipeline']]
        sig = (p is not None) and (not np.isnan(p)) and (p < 0.05)

        ax.errorbar(delta, y, xerr=[[delta - ci_lo], [ci_hi - delta]],
                    fmt='o', color=col, ecolor=col, elinewidth=1.4,
                    capsize=3.5, markersize=8,
                    markerfacecolor=col if sig else 'white',
                    markeredgecolor=col, markeredgewidth=1.4, zorder=5)

    # Zero reference line
    ax.axvline(x=0, color='black', linestyle='--', linewidth=0.8, alpha=0.7, zorder=2)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel('Mean paired delta (AMBER - reference) with 95% CI', fontsize=10)
    ax.set_title('AMBER paired effect, post-Rosetta, locked frame (n=257 per cell)',
                 fontsize=12, fontweight='bold')

    # Annotation
    p_min = df['paired_t_p'].min()
    p_max = df['paired_t_p'].max()
    n_sig = int((df['paired_t_p'] < 0.05).sum())
    ax.text(0.02, 0.99,
            f'{n_sig} of {n} cells reach p<0.05.\n'
            f'p-value range: {p_min:.3f} to {p_max:.3f}.\n'
            'All 95% CIs straddle zero at the locked\n'
            'post-Rosetta frame: paired AMBER effect\n'
            'is in the noise.',
            transform=ax.transAxes, fontsize=8.5, va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#f7f7f7',
                      edgecolor='#aaa', linewidth=0.6))

    # Legend
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=PIPELINE_COLORS['blue'],
               markeredgecolor=PIPELINE_COLORS['blue'], markersize=8, label='Blue, p>=0.05',
               markeredgewidth=1.4),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=PIPELINE_COLORS['green'],
               markeredgecolor=PIPELINE_COLORS['green'], markersize=8, label='Green, p>=0.05',
               markeredgewidth=1.4),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='white',
               markeredgecolor='black', markersize=8, label='hollow = p>=0.05',
               markeredgewidth=1.4),
    ]
    ax.legend(handles=legend_handles, loc='lower right', fontsize=8, frameon=True)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    save_fig(fig, 'fig22_amber_paired_forest')


def save_fig(fig, name):
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(FIGDIR, f"{name}.{fmt}"),
                    dpi=300 if fmt == 'png' else None,
                    bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {name}")


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    df = load()
    print(f"Loaded {len(df)} paired cells.")
    print(df[['pipeline', 'pair', 'protocol', 'metric', 'mean_diff',
              'diff_ci_lo', 'diff_ci_hi', 'paired_t_p']].to_string(index=False))
    fig22_paired_forest(df)
    print("Done.")


if __name__ == '__main__':
    main()
