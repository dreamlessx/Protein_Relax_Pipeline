#!/usr/bin/env python3
import pandas as pd
import glob
import re

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

def base_source(s):
    if s.startswith('af_relaxed'):
        return 'af_relaxed'
    if s.startswith('af_unrelaxed'):
        return 'af_unrelaxed'
    if s.startswith('amber_af'):
        return 'amber_af'
    if s.startswith('amber_boltz'):
        return 'amber_boltz'
    if s.startswith('amber_crystal'):
        return 'amber_crystal'
    if s.startswith('boltz'):
        return 'boltz'
    if s.startswith('crystal'):
        return 'crystal'
    return s

out['source'] = out['src_type'].apply(base_source)
OUT = '/data/p_csb_meiler/agarwm5/red_analysis/metrics/combined_rosetta_molprobity.tsv'
out.to_csv(OUT, sep='\t', index=False)
print(f'Wrote {len(out)} rows to {OUT}')
print('Source counts:')
print(out.groupby(['pipeline','source']).size())
