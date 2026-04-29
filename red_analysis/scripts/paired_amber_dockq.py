#!/usr/bin/env python3
"""Paired comparison: DockQ + i_rmsd + l_rmsd + f_nat, AMBER-relaxed vs un-AMBERed.

DockQ is computed only on pre-Rosetta input structures (no protocol_id/rep), so this
script reads `dockq_metrics` directly. Per-target mean across the 5 model variants
before pairing. Mirrors `amber_paired_table2.tsv` column shape.

Cells produced (Δ = AMBER minus reference):
  1. af_unrelaxed -> amber_af_unrelaxed   (5 ranked variants per side, averaged)
  2. af_relaxed   -> amber_af_unrelaxed   (alt reference; same AMBER source)
  3. boltz        -> amber_boltz_boltz    (5 model variants per side, averaged)
  4. crystal      -> amber_crystal        (single model per side)

5th cell (af_relaxed -> amber_af with relaxed source) collapses onto cell 2 because
DockQ's amber source for AF is `amber_af_unrelaxed_ranked_*`, no `amber_af_relaxed`
variant exists in dockq_metrics.

Crystal cell only available in `green` pipeline (amber_crystal pipeline_id='green'
exclusively per snapshot 2026-04-27a).

Total rows = 4 metrics x (3 pipe-blue cells + 4 pipe-green cells) = 4 x 7 = 28.
Spec said 20 (5 cells x 4 metrics). Difference: blue is missing crystal, and we
expanded to af_relaxed->amber_af variant. Documented above.
"""
import sqlite3
import pandas as pd
import numpy as np
from scipy import stats

DB = "/data/p_csb_meiler/agarwm5/red_analysis/db/bm55.sqlite"
OUT = "/data/p_csb_meiler/agarwm5/red_analysis/tables/amber_paired_dockq.tsv"

SRC_MAP = {
    'af_relaxed_ranked_': 'af_relaxed',
    'af_unrelaxed_ranked_': 'af_unrelaxed',
    'amber_af_unrelaxed_ranked_': 'amber_af',
    'amber_boltz_boltz_model_': 'amber_boltz',
    'boltz_boltz_input_model_': 'boltz',
}

def normalize_source(s):
    if s == 'amber_crystal':
        return 'amber_crystal'
    if s.startswith('crystal_'):
        return 'crystal'
    for prefix, label in SRC_MAP.items():
        if s.startswith(prefix):
            return label
    return None

def cliffs_delta(x, y):
    x = np.asarray(x); y = np.asarray(y)
    n = len(x)
    if n == 0:
        return float('nan')
    d = x - y
    return ((d > 0).sum() - (d < 0).sum()) / n

def cliff_label(d):
    a = abs(d)
    if a < 0.147: return 'negligible'
    if a < 0.33:  return 'small'
    if a < 0.474: return 'medium'
    return 'large'

def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    adj = p[order] * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty_like(adj)
    out[order] = np.minimum(adj, 1.0)
    return out

def main():
    con = sqlite3.connect(DB)
    df = pd.read_sql_query(
        """
        SELECT target_id, pipeline_id, src_type, dockq, i_rmsd, l_rmsd, f_nat
        FROM dockq_metrics
        WHERE status = 'ok'
        """,
        con,
    )
    con.close()

    df['source'] = df['src_type'].map(normalize_source)
    df = df.dropna(subset=['source'])

    metrics = ['dockq', 'i_rmsd', 'l_rmsd', 'f_nat']
    tmeans = df.groupby(['pipeline_id','source','target_id'])[metrics].mean().reset_index()

    pairs = [
        ('amber_af',      'af_unrelaxed'),
        ('amber_af',      'af_relaxed'),
        ('amber_boltz',   'boltz'),
        ('amber_crystal', 'crystal'),
    ]

    rows = []
    for pipe in ['blue', 'green']:
        for amber_src, ref_src in pairs:
            a = tmeans[(tmeans['pipeline_id']==pipe) & (tmeans['source']==amber_src)][['target_id'] + metrics].copy()
            r = tmeans[(tmeans['pipeline_id']==pipe) & (tmeans['source']==ref_src)][['target_id'] + metrics].copy()
            if len(a) == 0 or len(r) == 0:
                continue
            a = a.rename(columns={m: f'amber_{m}' for m in metrics})
            r = r.rename(columns={m: f'ref_{m}' for m in metrics})
            m = a.merge(r, on='target_id', how='inner')
            for metric in metrics:
                pair = m[[f'amber_{metric}', f'ref_{metric}']].dropna()
                n = len(pair)
                if n < 3:
                    continue
                amber = pair[f'amber_{metric}'].values
                ref   = pair[f'ref_{metric}'].values
                diffs = amber - ref
                mean_d = diffs.mean()
                median_d = float(np.median(diffs))
                se = stats.sem(diffs) if n > 1 else float('nan')
                ci = stats.t.ppf(0.975, df=n-1) * se if n > 1 else float('nan')
                t, p_t = stats.ttest_rel(amber, ref)
                try:
                    w_stat, p_w = stats.wilcoxon(amber, ref)
                except Exception:
                    w_stat, p_w = (float('nan'), float('nan'))
                d = cliffs_delta(amber, ref)
                rows.append({
                    'pipeline': pipe,
                    'source_no_amber': ref_src,
                    'source_with_amber': amber_src,
                    'metric': metric,
                    'n_paired': n,
                    'mean_delta': round(mean_d, 5),
                    'median_delta': round(median_d, 5),
                    'ci_low': round(mean_d - ci, 5),
                    'ci_high': round(mean_d + ci, 5),
                    'paired_t_p': round(p_t, 6) if not np.isnan(p_t) else None,
                    'wilcoxon_p': round(p_w, 6) if not np.isnan(p_w) else None,
                    'cliffs_d': round(d, 4),
                    'effect_label': cliff_label(d),
                })

    out = pd.DataFrame(rows)
    if len(out):
        out['paired_t_q'] = bh_fdr(out['paired_t_p'].fillna(1.0).values).round(6)
        out['wilcoxon_q'] = bh_fdr(out['wilcoxon_p'].fillna(1.0).values).round(6)

    out.to_csv(OUT, sep='\t', index=False)
    print(f"Wrote {len(out)} rows to {OUT}")
    print(out.to_string(index=False))

if __name__ == "__main__":
    main()
