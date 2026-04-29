#!/usr/bin/env python3
"""Paired comparison: Rosetta total_score (energy), AMBER vs un-AMBERed.

Reads `rosetta_energy` table directly. Per-target mean across 5 reps x 6 protocols
(or as many as exist) before pairing. Mirrors `amber_paired_table2.tsv` column shape
plus Cliff's d + effect_label + FDR.

Cells produced (Δ = AMBER minus reference):
  1. af_unrelaxed -> amber_af   (5 ranked variants per side, averaged across reps and protocols)
  2. af_relaxed   -> amber_af   (alt reference)
  3. boltz        -> amber_boltz
  4. crystal      -> amber_crystal_relaxed

5th cell collapses (no separate amber_af_relaxed src in rosetta_energy).
n_paired = 257 per cell. 4 cells x 2 pipelines x 1 metric = 8 rows total.
Spec said "5 paired cells, 5 rows" - difference is documented above; pipeline axis
adds another factor since pipelines run independently.
"""
import sqlite3
import pandas as pd
import numpy as np
from scipy import stats

DB  = "/data/p_csb_meiler/agarwm5/red_analysis/db/bm55.sqlite"
OUT = "/data/p_csb_meiler/agarwm5/red_analysis/tables/rosetta_energy_paired.tsv"

SRC_MAP = {
    'af_relaxed_ranked_': 'af_relaxed',
    'af_unrelaxed_ranked_': 'af_unrelaxed',
    'amber_af_ranked_': 'amber_af',
    'amber_boltz_model_': 'amber_boltz',
    'boltz_boltz_input_model_': 'boltz',
}

def normalize_source(s):
    if s == 'amber_crystal_relaxed':
        return 'amber_crystal_relaxed'
    if s.startswith('crystal_'):
        return 'crystal'
    for prefix, label in SRC_MAP.items():
        if s.startswith(prefix):
            return label
    return None

def cliffs_delta(x, y):
    x = np.asarray(x); y = np.asarray(y)
    n = len(x)
    if n == 0: return float('nan')
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
        SELECT target_id, pipeline_id, src_type, protocol_id, rep, total_score
        FROM rosetta_energy
        """,
        con,
    )
    con.close()

    df['source'] = df['src_type'].map(normalize_source)
    df = df.dropna(subset=['source','total_score'])

    # per-target mean across 5 reps x ~6 protocols
    tmeans = df.groupby(['pipeline_id','source','target_id'])['total_score'].mean().reset_index()

    pairs = [
        ('amber_af',                'af_unrelaxed'),
        ('amber_af',                'af_relaxed'),
        ('amber_boltz',             'boltz'),
        ('amber_crystal_relaxed',   'crystal'),
    ]

    rows = []
    for pipe in ['blue', 'green']:
        for amber_src, ref_src in pairs:
            a = tmeans[(tmeans['pipeline_id']==pipe) & (tmeans['source']==amber_src)][['target_id','total_score']].rename(columns={'total_score':'amber'})
            r = tmeans[(tmeans['pipeline_id']==pipe) & (tmeans['source']==ref_src)][['target_id','total_score']].rename(columns={'total_score':'ref'})
            m = a.merge(r, on='target_id', how='inner').dropna()
            n = len(m)
            if n < 3:
                continue
            amber = m['amber'].values
            ref   = m['ref'].values
            diffs = amber - ref
            mean_d = diffs.mean()
            median_d = float(np.median(diffs))
            se = stats.sem(diffs)
            ci = stats.t.ppf(0.975, df=n-1) * se
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
                'metric': 'total_score',
                'n_paired': n,
                'amber_mean': round(float(amber.mean()), 4),
                'ref_mean':   round(float(ref.mean()), 4),
                'mean_delta': round(mean_d, 4),
                'median_delta': round(median_d, 4),
                'ci_low':  round(mean_d - ci, 4),
                'ci_high': round(mean_d + ci, 4),
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
