#!/usr/bin/env python3
"""Compute per-target MolProbity means for paper reporting.
Per-target = average 5 reps within each (target, pipeline, source, protocol) first,
then compute mean/CI across targets.
"""
import pandas as pd
import numpy as np
from scipy import stats

TSV = "/data/p_csb_meiler/agarwm5/red_analysis/metrics/combined_rosetta_molprobity.tsv"
OUT = "/data/p_csb_meiler/agarwm5/red_analysis/tables/per_target_mp_means.tsv"

df = pd.read_csv(TSV, sep='\t')
print(f"Loaded {len(df)} rows")
print(f"Columns: {list(df.columns)}")
print(f"Unique pipelines: {df['pipeline'].unique()}")
print(f"Unique sources: {df['source'].unique() if 'source' in df.columns else 'N/A'}")

# Use src_type if no 'source' col
src_col = 'source' if 'source' in df.columns else 'src_type'
print(f"Source column: {src_col}")
print(f"Unique {src_col}: {df[src_col].unique()}")
print(f"Unique protocols: {df['protocol'].unique()}")
print(f"Unique targets: {df['target'].nunique()}")

# Step 1: per-target means (average across reps)
target_means = df.groupby(['pipeline', src_col, 'protocol', 'target'])[
    ['clashscore', 'molprobity_score', 'rama_favored', 'rota_outliers']
].mean().reset_index()

# Step 2: mean/CI across targets
def ci95(x):
    n = len(x)
    if n < 2:
        return (x.mean(), np.nan, np.nan, n)
    se = stats.sem(x)
    ci = stats.t.ppf(0.975, df=n-1) * se
    return (x.mean(), x.mean() - ci, x.mean() + ci, n)

rows = []
for (pipe, src, proto), grp in target_means.groupby(['pipeline', src_col, 'protocol']):
    for metric in ['clashscore', 'molprobity_score', 'rota_outliers', 'rama_favored']:
        vals = grp[metric].dropna()
        if len(vals) < 1:
            continue
        m, lo, hi, n = ci95(vals)
        rows.append({
            'pipeline': pipe,
            'source': src,
            'protocol': proto,
            'metric': metric,
            'mean': round(m, 4),
            'ci_lo': round(lo, 4) if not np.isnan(lo) else None,
            'ci_hi': round(hi, 4) if not np.isnan(hi) else None,
            'n_targets': n
        })

result = pd.DataFrame(rows)
result.to_csv(OUT, sep='\t', index=False)
print(f"\nSaved {len(result)} rows to {OUT}")

# Print key paper metrics
print("\n=== KEY PAPER METRICS (per-target means, MolProbity score) ===")
mp = result[result['metric'] == 'molprobity_score'].sort_values(['pipeline', 'source', 'protocol'])
print(mp[['pipeline','source','protocol','mean','ci_lo','ci_hi','n_targets']].to_string(index=False))

print("\n=== KEY PAPER METRICS (per-target means, clashscore) ===")
cs = result[result['metric'] == 'clashscore'].sort_values(['mean'])
print(cs[['pipeline','source','protocol','mean','ci_lo','ci_hi','n_targets']].head(20).to_string(index=False))

