#!/usr/bin/env python3
"""Paired rotamer outlier comparison: AMBER(x)+Rosetta vs x+Rosetta on matched targets.
Extension of paired_amber.py for Table 2 rotamer subsection (iter 32).
Targets the crystal+normal cells flagged in S5 Obs 4 as largest AMBER-driven rotamer drops."""
import pandas as pd
import numpy as np
from scipy import stats

TSV = "/data/p_csb_meiler/agarwm5/red_analysis/metrics/combined_rosetta_molprobity.tsv"
OUT = "/data/p_csb_meiler/agarwm5/red_analysis/tables/amber_paired_rotamer.tsv"

df = pd.read_csv(TSV, sep='\t')
src_col = 'source' if 'source' in df.columns else 'src_type'

tmeans = df.groupby(['pipeline', src_col, 'protocol', 'target'])[
    ['rota_outliers']
].mean().reset_index()

pairs = [('amber_boltz', 'boltz'), ('amber_af', 'af_relaxed'), ('amber_crystal', 'crystal')]
protocols = ['normal_beta', 'normal_ref15', 'dualspace_beta']
metric = 'rota_outliers'

rows = []
for pipe in ['blue', 'green']:
    for amber_src, std_src in pairs:
        for proto in protocols:
            a = tmeans[(tmeans['pipeline']==pipe) & (tmeans[src_col]==amber_src) & (tmeans['protocol']==proto)][['target', metric]].rename(columns={metric: 'amber'})
            s = tmeans[(tmeans['pipeline']==pipe) & (tmeans[src_col]==std_src) & (tmeans['protocol']==proto)][['target', metric]].rename(columns={metric: 'std'})
            m = a.merge(s, on='target', how='inner').dropna()
            n = len(m)
            if n < 3:
                continue
            diffs = m['amber'] - m['std']
            mean_diff = diffs.mean()
            se = stats.sem(diffs) if n > 1 else float('nan')
            ci = stats.t.ppf(0.975, df=n-1) * se if n > 1 else float('nan')
            t, p = stats.ttest_rel(m['amber'], m['std'])
            try:
                w_stat, w_p = stats.wilcoxon(m['amber'], m['std'])
            except Exception:
                w_stat, w_p = (float('nan'), float('nan'))
            rows.append({
                'pipeline': pipe, 'pair': f'{amber_src}_vs_{std_src}',
                'protocol': proto, 'metric': metric,
                'n_paired': n,
                'amber_mean': round(m['amber'].mean(), 5),
                'std_mean': round(m['std'].mean(), 5),
                'mean_diff': round(mean_diff, 5),
                'diff_ci_lo': round(mean_diff - ci, 5),
                'diff_ci_hi': round(mean_diff + ci, 5),
                'paired_t_p': round(p, 5),
                'wilcoxon_p': round(w_p, 5) if not np.isnan(w_p) else None,
            })

out = pd.DataFrame(rows)
out.to_csv(OUT, sep='\t', index=False)
print(f"Wrote {len(out)} rows to {OUT}")
print(out.to_string(index=False))
