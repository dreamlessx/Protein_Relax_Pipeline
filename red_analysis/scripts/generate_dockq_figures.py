#!/usr/bin/env python3
"""
generate_dockq_figures.py — Red Analysis Pipeline

fig21_dockq_by_source: DockQ box+strip across the 7 source buckets, hue by pipeline,
paired-line overlays for the AMBER pairs.

Corrected narrative: AMBER preconditioning preserves the binding interface on AF and
Boltz inputs (mean delta < 1e-3) but introduces a small perturbation on crystal
references (delta = -0.022, p < 1e-17).
"""

import os
import sqlite3
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

DB = "/data/p_csb_meiler/agarwm5/red_analysis/db/bm55.sqlite"
FIGDIR = "/data/p_csb_meiler/agarwm5/red_analysis/figures"

PIPELINE_COLORS = {'blue': '#1f6feb', 'green': '#2ea043'}

SOURCE_ORDER = ['crystal', 'af_unrelaxed', 'af_relaxed', 'amber_af',
                'boltz', 'amber_boltz', 'amber_crystal']
SOURCE_LABELS = {
    'crystal':       'Crystal',
    'af_unrelaxed':  'AF2\n(unrelaxed)',
    'af_relaxed':    'AF2\n(relaxed)',
    'amber_af':      'AMBER\n(AF)',
    'boltz':         'Boltz-1',
    'amber_boltz':   'AMBER\n(Boltz)',
    'amber_crystal': 'AMBER\n(Crystal)',
}

BUCKET_SQL = """
SELECT
  target_id,
  pipeline_id AS pipeline,
  CASE
    WHEN src_type = 'amber_crystal' THEN 'amber_crystal'
    WHEN src_type LIKE 'crystal%' OR src_type = 'crystal' THEN 'crystal'
    WHEN src_type LIKE 'amber_af%' THEN 'amber_af'
    WHEN src_type LIKE 'amber_boltz%' THEN 'amber_boltz'
    WHEN src_type LIKE 'af_relaxed%' THEN 'af_relaxed'
    WHEN src_type LIKE 'af_unrelaxed%' THEN 'af_unrelaxed'
    WHEN src_type LIKE 'boltz%' THEN 'boltz'
    ELSE NULL
  END AS source,
  dockq
FROM dockq_metrics
WHERE dockq IS NOT NULL
"""


def load_data():
    con = sqlite3.connect(DB)
    df = pd.read_sql(BUCKET_SQL, con)
    con.close()
    df = df.dropna(subset=['source'])
    return df


def per_target_means(df):
    """Mean DockQ per (target, pipeline, source). Collapses 1 or 5 model-replicates."""
    return df.groupby(['target_id', 'pipeline', 'source'], as_index=False)['dockq'].mean()


def paired_summary(per_t):
    """Return long-format paired delta means for the 5 AMBER pairs.

    Pairs (ref -> amber_x):
      af_unrelaxed -> amber_af
      af_relaxed   -> amber_af
      boltz        -> amber_boltz
      crystal      -> amber_crystal
    Per pipeline.
    """
    pairs = [
        ('af_unrelaxed', 'amber_af'),
        ('af_relaxed',   'amber_af'),
        ('boltz',        'amber_boltz'),
        ('crystal',      'amber_crystal'),
    ]
    rows = []
    for pipe in ['blue', 'green']:
        for ref, amber in pairs:
            sub_ref = per_t[(per_t.pipeline == pipe) & (per_t.source == ref)].set_index('target_id')['dockq']
            sub_amb = per_t[(per_t.pipeline == pipe) & (per_t.source == amber)].set_index('target_id')['dockq']
            common = sub_ref.index.intersection(sub_amb.index)
            if len(common) < 2:
                continue
            d = sub_amb.loc[common].values - sub_ref.loc[common].values
            t_p = stats.ttest_rel(sub_amb.loc[common].values,
                                  sub_ref.loc[common].values).pvalue
            rows.append({
                'pipeline': pipe,
                'ref': ref,
                'amber': amber,
                'n': len(common),
                'mean_delta': float(np.mean(d)),
                'paired_t_p': float(t_p),
            })
    return pd.DataFrame(rows)


def fig21_dockq_by_source(per_t, paired):
    fig, ax = plt.subplots(figsize=(13, 6.5))

    # Box plot with hue=pipeline
    plot_df = per_t[per_t.source.isin(SOURCE_ORDER)].copy()
    plot_df['source'] = pd.Categorical(plot_df['source'],
                                       categories=SOURCE_ORDER, ordered=True)

    sns.boxplot(data=plot_df, x='source', y='dockq', hue='pipeline',
                order=SOURCE_ORDER, hue_order=['blue', 'green'],
                palette=PIPELINE_COLORS, fliersize=0,
                width=0.7, linewidth=0.9, ax=ax)

    # Strip overlay (jitter)
    sns.stripplot(data=plot_df, x='source', y='dockq', hue='pipeline',
                  order=SOURCE_ORDER, hue_order=['blue', 'green'],
                  palette=PIPELINE_COLORS, alpha=0.25, size=2.2,
                  dodge=True, jitter=0.18, ax=ax, legend=False)

    # Paired-line overlays for the 4 AMBER pairs (showing means)
    pair_to_pos = {
        ('af_unrelaxed', 'amber_af'):    (1, 3),
        ('af_relaxed',   'amber_af'):    (2, 3),
        ('boltz',        'amber_boltz'): (4, 5),
        ('crystal',      'amber_crystal'): (0, 6),
    }

    for _, r in paired.iterrows():
        key = (r['ref'], r['amber'])
        if key not in pair_to_pos:
            continue
        x0, x1 = pair_to_pos[key]
        # dodge offset: blue left (-0.18), green right (+0.18)
        dx = -0.18 if r['pipeline'] == 'blue' else 0.18
        col = PIPELINE_COLORS[r['pipeline']]
        # Take the mean DockQ at each end
        ref_mean = plot_df[(plot_df.pipeline == r['pipeline']) &
                           (plot_df.source == r['ref'])]['dockq'].mean()
        amb_mean = plot_df[(plot_df.pipeline == r['pipeline']) &
                           (plot_df.source == r['amber'])]['dockq'].mean()
        ax.plot([x0 + dx, x1 + dx], [ref_mean, amb_mean],
                color=col, linewidth=1.6, alpha=0.85,
                marker='o', markersize=5, markerfacecolor='white',
                markeredgecolor=col, markeredgewidth=1.4, zorder=10)

    # X-axis labels
    ax.set_xticks(range(len(SOURCE_ORDER)))
    ax.set_xticklabels([SOURCE_LABELS[s] for s in SOURCE_ORDER], fontsize=9)
    ax.set_xlabel('Source bucket', fontsize=11)
    ax.set_ylabel('DockQ', fontsize=11)
    ax.set_ylim(-0.02, 1.05)
    ax.axhline(y=0.23, color='gray', linestyle=':', linewidth=0.7, alpha=0.5)
    ax.text(6.25, 0.24, 'acceptable threshold', fontsize=7,
            color='gray', va='bottom', ha='right')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Headline annotation
    crystal_paired = paired[(paired.ref == 'crystal') &
                            (paired.amber == 'amber_crystal')]
    if len(crystal_paired) > 0:
        c = crystal_paired.iloc[0]
        ax.annotate(f"crystal -> amber_crystal:\nDelta={c['mean_delta']:.4f}, p<1e-17",
                    xy=(6, 0.978), xytext=(5.0, 0.55),
                    fontsize=8, color='#222',
                    arrowprops=dict(arrowstyle='->', lw=0.8, color='#444'))

    # Custom legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], ['Blue pipeline', 'Green pipeline'],
              loc='lower left', frameon=True, fontsize=9)

    fig.suptitle('DockQ across input sources, by pipeline',
                 fontsize=13, fontweight='bold')
    fig.text(0.5, 0.01,
             'AMBER preconditioning preserves the interface on AF/Boltz inputs '
             '(mean DeltaDockQ < 0.001) but perturbs the crystal reference '
             '(DeltaDockQ = -0.022, p < 1e-17).',
             fontsize=8.5, ha='center', style='italic', color='#333')
    plt.tight_layout(rect=[0, 0.04, 1, 0.97])
    save_fig(fig, 'fig21_dockq_by_source')


def save_fig(fig, name):
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(FIGDIR, f"{name}.{fmt}"),
                    dpi=300 if fmt == 'png' else None,
                    bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {name}")


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    print("Loading DockQ data from DB...")
    raw = load_data()
    per_t = per_target_means(raw)
    print(f"  per-target rows: {len(per_t)}")

    paired = paired_summary(per_t)
    print("Paired summary:")
    print(paired.to_string(index=False))

    print("Generating fig21...")
    fig21_dockq_by_source(per_t, paired)
    print("Done.")


if __name__ == '__main__':
    main()
