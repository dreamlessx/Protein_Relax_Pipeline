#!/usr/bin/env python3
"""Paired comparison: AMBER(x)+Rosetta vs standalone x+Rosetta on matching targets.
For Table 2 in manuscript."""
import pandas as pd
import numpy as np
from scipy import stats

TSV = "/data/p_csb_meiler/agarwm5/red_analysis/metrics/combined_rosetta_molprobity.tsv"
OUT = "/data/p_csb_meiler/agarwm5/red_analysis/tables/amber_paired_table2.tsv"

df = pd.read_csv(TSV, sep='\t')
src_col = 'source' if 'source' in df.columns else 'src_type'

# per-target means (average reps)
tmeans = df.groupby(['pipeline', src_col, 'protocol', 'target'])[
    ['clashscore', 'molprobity_score', 'rama_favored', 'rota_outliers']
].mean().reset_index()

# pairs: (amber_source, standalone_source)
pairs = [('amber_boltz', 'boltz'), ('amber_af', 'af_relaxed')]
protocols = ['normal_beta', 'dualspace_beta']
metrics = ['molprobity_score', 'clashscore']

rows = []
for pipe in ['blue', 'green']:
    for amber_src, std_src in pairs:
        for proto in protocols:
            for metric in metrics:
                a = tmeans[(tmeans['pipeline']==pipe) & (tmeans[src_col]==amber_src) & (tmeans['protocol']==proto)][['target', metric]].rename(columns={metric: 'amber'})
                s = tmeans[(tmeans['pipeline']==pipe) & (tmeans[src_col]==std_src) & (tmeans['protocol']==proto)][['target', metric]].rename(columns={metric: 'std'})
                m = a.merge(s, on='target', how='inner').dropna()
                n = len(m)
                if n < 3:
                    continue
                diffs = m['amber'] - m['std']
                mean_diff = diffs.mean()
                se = stats.sem(diffs)
                ci = stats.t.ppf(0.975, df=n-1) * se
                t, p = stats.ttest_rel(m['amber'], m['std'])
                w_stat, w_p = stats.wilcoxon(m['amber'], m['std']) if n >= 6 else (np.nan, np.nan)
                rows.append({
                    'pipeline': pipe, 'pair': f'{amber_src}_vs_{std_src}',
                    'protocol': proto, 'metric': metric,
                    'n_paired': n,
                    'amber_mean': round(m['amber'].mean(), 4),
                    'std_mean': round(m['std'].mean(), 4),
                    'mean_diff': round(mean_diff, 4),
                    'diff_ci_lo': round(mean_diff - ci, 4),
                    'diff_ci_hi': round(mean_diff + ci, 4),
                    'paired_t_p': round(p, 5),
                    'wilcoxon_p': round(w_p, 5) if not np.isnan(w_p) else None,
                })

out = pd.DataFrame(rows)
out.to_csv(OUT, sep='\t', index=False)
print(f"Wrote {len(out)} rows to {OUT}")
print(out.to_string(index=False))
