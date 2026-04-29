#!/usr/bin/env python3
"""Difficulty-stratified version of `paired_amber_rotamer.py` (the 18-cell rotamer extension).

3 difficulty tiers x 18 cells = 54 rows. Same column shape as
amber_paired_table2_by_difficulty.tsv with Cliff's d + effect_label + FDR.
"""
import sqlite3
import pandas as pd
import numpy as np
from scipy import stats

TSV = "/data/p_csb_meiler/agarwm5/red_analysis/metrics/combined_rosetta_molprobity.tsv"
DB  = "/data/p_csb_meiler/agarwm5/red_analysis/db/bm55.sqlite"
OUT = "/data/p_csb_meiler/agarwm5/red_analysis/tables/amber_paired_rotamer_by_difficulty.tsv"

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
    # The combined TSV upstream coerced target '1E96' -> '1e+96' via float inference.
    # Patch on the way in so the difficulty merge catches all 257 targets.
    df = pd.read_csv(TSV, sep='\t', dtype={'target': str})
    df['target'] = df['target'].replace({'1e+96': '1E96'})
    src_col = 'source' if 'source' in df.columns else 'src_type'

    con = sqlite3.connect(DB)
    diff_df = pd.read_sql_query("SELECT target_id AS target, difficulty FROM targets", con)
    con.close()
    df = df.merge(diff_df, on='target', how='left')

    tmeans = df.groupby(['pipeline', src_col, 'protocol', 'target', 'difficulty'])[
        ['rota_outliers']
    ].mean().reset_index()

    pairs = [('amber_boltz','boltz'), ('amber_af','af_relaxed'), ('amber_crystal','crystal')]
    protocols = ['normal_beta', 'normal_ref15', 'dualspace_beta']
    metric = 'rota_outliers'
    difficulties = ['rigid', 'medium', 'difficult']

    rows = []
    for difficulty in difficulties:
        for pipe in ['blue', 'green']:
            for amber_src, std_src in pairs:
                for proto in protocols:
                    sub = tmeans[(tmeans['difficulty']==difficulty) & (tmeans['pipeline']==pipe) & (tmeans['protocol']==proto)]
                    a = sub[sub[src_col]==amber_src][['target', metric]].rename(columns={metric:'amber'})
                    s = sub[sub[src_col]==std_src][['target', metric]].rename(columns={metric:'std'})
                    m = a.merge(s, on='target', how='inner').dropna()
                    n = len(m)
                    if n < 3:
                        continue
                    diffs = m['amber'].values - m['std'].values
                    mean_d = diffs.mean()
                    median_d = float(np.median(diffs))
                    se = stats.sem(diffs) if n > 1 else float('nan')
                    ci = stats.t.ppf(0.975, df=n-1) * se if n > 1 else float('nan')
                    t, p_t = stats.ttest_rel(m['amber'], m['std'])
                    try:
                        w_stat, p_w = stats.wilcoxon(m['amber'], m['std'])
                    except Exception:
                        w_stat, p_w = (float('nan'), float('nan'))
                    d = cliffs_delta(m['amber'].values, m['std'].values)
                    rows.append({
                        'difficulty': difficulty,
                        'pipeline': pipe,
                        'pair': f'{amber_src}_vs_{std_src}',
                        'protocol': proto,
                        'metric': metric,
                        'n_paired': n,
                        'amber_mean': round(m['amber'].mean(), 5),
                        'std_mean':   round(m['std'].mean(), 5),
                        'mean_delta': round(mean_d, 5),
                        'median_delta': round(median_d, 5),
                        'ci_low':  round(mean_d - ci, 5),
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
