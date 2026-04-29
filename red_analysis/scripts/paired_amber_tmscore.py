#!/usr/bin/env python3
"""Paired comparison: pre-Rosetta TM-score, AMBER-relaxed input vs un-AMBERed input.

Pairs by target. For each (pipeline, paired-cell), per-target mean over the 5 ranked
variants is taken before pairing. Output mirrors `amber_paired_table2.tsv` column shape
plus `cliffs_d` and `effect_label` per task spec.

Choice of pre-Rosetta TM (`is_post_rosetta=0`):
The headline AMBER-on-TM claim is "AMBER doesn't change TM". Pre-Rosetta TM measures
the input structure's TM relative to crystal *before* Rosetta has had a chance to
re-shuffle backbones. This is the cleanest test of AMBER alone. Post-Rosetta TM
confounds AMBER-input-perturbation with downstream Rosetta protocol behaviour.

Cells produced (Δ = AMBER minus reference):
  1. af_unrelaxed -> amber_af   (AMBER applied to unrelaxed AF as input)
  2. af_relaxed   -> amber_af   (AlphaFold-internal AMBER baseline vs explicit AMBER)
  3. boltz        -> amber_boltz
  (crystal pair skipped: TM(crystal, crystal) is 1.0 trivially; amber_crystal not in
   tm_scores at all)

n_paired per cell = 257 (all targets).
Total rows = 3 cells x 2 pipelines x 1 metric = 6.
"""
import sqlite3
import sys
import pandas as pd
import numpy as np
from scipy import stats

DB = "/data/p_csb_meiler/agarwm5/red_analysis/db/bm55.sqlite"
OUT = "/data/p_csb_meiler/agarwm5/red_analysis/tables/amber_paired_tmscore.tsv"

# Map raw src_type prefix -> normalized source label
SRC_MAP = {
    'af_relaxed_ranked_': 'af_relaxed',
    'af_unrelaxed_ranked_': 'af_unrelaxed',
    'amber_af_ranked_': 'amber_af',
    'amber_boltz_model_': 'amber_boltz',
    'boltz_boltz_input_model_': 'boltz',
}

def normalize_source(s):
    for prefix, label in SRC_MAP.items():
        if s.startswith(prefix):
            return label
    return None

def cliffs_delta(x, y):
    """Cliff's delta on paired arrays: counts of x_i > y_i vs x_i < y_i / N."""
    x = np.asarray(x); y = np.asarray(y)
    n = len(x)
    if n == 0:
        return float('nan')
    d = x - y
    gt = (d > 0).sum()
    lt = (d < 0).sum()
    return (gt - lt) / n

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
    ranked = p[order]
    adj = ranked * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest down
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty_like(adj)
    out[order] = np.minimum(adj, 1.0)
    return out

def main():
    con = sqlite3.connect(DB)
    df = pd.read_sql_query(
        """
        SELECT target_id, pipeline_id, src_type, tm_score
        FROM tm_scores
        WHERE is_post_rosetta = 0
        """,
        con,
    )
    con.close()

    df['source'] = df['src_type'].map(normalize_source)
    df = df.dropna(subset=['source'])

    # per-target mean across the 5 ranked variants
    tmeans = df.groupby(['pipeline_id', 'source', 'target_id'])['tm_score'].mean().reset_index()

    pairs = [
        ('amber_af',    'af_unrelaxed'),
        ('amber_af',    'af_relaxed'),
        ('amber_boltz', 'boltz'),
    ]

    rows = []
    for pipe in ['blue', 'green']:
        for amber_src, ref_src in pairs:
            a = tmeans[(tmeans['pipeline_id']==pipe) & (tmeans['source']==amber_src)][['target_id','tm_score']].rename(columns={'tm_score':'amber'})
            r = tmeans[(tmeans['pipeline_id']==pipe) & (tmeans['source']==ref_src)][['target_id','tm_score']].rename(columns={'tm_score':'ref'})
            m = a.merge(r, on='target_id', how='inner').dropna()
            n = len(m)
            if n < 3:
                continue
            diffs = m['amber'].values - m['ref'].values
            mean_d = diffs.mean()
            median_d = float(np.median(diffs))
            se = stats.sem(diffs)
            ci = stats.t.ppf(0.975, df=n-1) * se
            t, p_t = stats.ttest_rel(m['amber'], m['ref'])
            try:
                w_stat, p_w = stats.wilcoxon(m['amber'], m['ref'])
            except Exception:
                w_stat, p_w = (float('nan'), float('nan'))
            d = cliffs_delta(m['amber'].values, m['ref'].values)
            rows.append({
                'pipeline': pipe,
                'source_no_amber': ref_src,
                'source_with_amber': amber_src,
                'metric': 'tm_score_pre_rosetta',
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
    # BH FDR across all rows on paired_t_p
    if len(out):
        out['paired_t_q'] = bh_fdr(out['paired_t_p'].fillna(1.0).values).round(6)
        out['wilcoxon_q'] = bh_fdr(out['wilcoxon_p'].fillna(1.0).values).round(6)

    out.to_csv(OUT, sep='\t', index=False)
    print(f"Wrote {len(out)} rows to {OUT}")
    print(out.to_string(index=False))

if __name__ == "__main__":
    main()
