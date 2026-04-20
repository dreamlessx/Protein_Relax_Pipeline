#!/usr/bin/env python3
"""
reaggregate_combined.py

Re-aggregates all per-target rosetta_molprobity_{pipeline}_{target}.tsv files
into a single combined TSV, EXCLUDING legacy collision source directories.

LEGACY SOURCES EXCLUDED:
- amber_af_relaxed   (pre-fix collision dir; pooled with af_relaxed+AMBER)
- amber_boltz_relaxed (pre-fix collision dir; pooled with boltz+AMBER)
These existed only for 2 Blue targets (1K5D, 4GAM) as residue from the
af_relaxed->boltz naming collision fix (see PROJECT_STATUS.md).

Source bucketing uses EXACT-prefix matching to avoid legacy-dir pooling.
"""
import pandas as pd
import glob

LEGACY_SRC = {'amber_af_relaxed', 'amber_boltz_relaxed'}

files = sorted(glob.glob('/data/p_csb_meiler/agarwm5/red_analysis/metrics/rosetta_molprobity_blue_*.tsv') +
               glob.glob('/data/p_csb_meiler/agarwm5/red_analysis/metrics/rosetta_molprobity_green_*.tsv'))
print(f'Files: {len(files)}')
dfs = []
for f in files:
    try:
        d = pd.read_csv(f, sep='\t')
        if len(d) > 0:
            dfs.append(d)
    except Exception as e:
        print(f'skip {f}: {e}')
out = pd.concat(dfs, ignore_index=True)

# Drop legacy collision rows BEFORE source bucketing
before_n = len(out)
out = out[~out['src_type'].isin(LEGACY_SRC)].copy()
dropped = before_n - len(out)
print(f'Dropped {dropped} legacy-source rows ({sorted(LEGACY_SRC)})')

def base_source(s):
    if s.startswith('af_relaxed'):
        return 'af_relaxed'
    if s.startswith('af_unrelaxed'):
        return 'af_unrelaxed'
    if s.startswith('amber_af_ranked'):
        return 'amber_af'
    if s.startswith('amber_boltz_model'):
        return 'amber_boltz'
    if s.startswith('amber_crystal'):
        return 'amber_crystal'
    if s.startswith('boltz_boltz_input'):
        return 'boltz'
    if s.startswith('crystal_'):
        return 'crystal'
    raise ValueError(f'Unknown src_type: {s}')

out['source'] = out['src_type'].apply(base_source)
OUT = '/data/p_csb_meiler/agarwm5/red_analysis/metrics/combined_rosetta_molprobity.tsv'
out.to_csv(OUT, sep='\t', index=False)
print(f'Wrote {len(out)} rows to {OUT}')
print('Source counts:')
print(out.groupby(['pipeline','source']).size())
