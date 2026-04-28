#!/usr/bin/env python3
"""
target_summary.py - Red Analysis Pipeline

Creates a per-target summary combining TM-score and MolProbity data.
Identifies which prediction source is "best" for each target and metric.
Useful for identifying patterns: which targets benefit from AMBER? Which
ones have poor predictions from all sources?

Output: target_summary.tsv with one row per target, columns for each
source's TM-score and MolProbity score, plus "winner" columns.
"""

import os
import numpy as np
import pandas as pd

METRICS_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"
OUTDIR = "/data/p_csb_meiler/agarwm5/red_analysis/tables"

SOURCE_ORDER = ['af_relaxed', 'af_unrelaxed', 'boltz', 'amber_af', 'amber_boltz']


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    pipe = 'blue'

    # Load TM-score
    tm_path = os.path.join(METRICS_DIR, "combined_tmscore.tsv")
    tm = pd.read_csv(tm_path, sep='\t')
    tm['tmscore'] = pd.to_numeric(tm['tmscore'], errors='coerce')
    tm_blue = tm[tm['pipeline'] == pipe]
    tm_by_target = tm_blue.groupby(['target', 'source'])['tmscore'].mean().unstack()

    # Load MolProbity
    mp_path = os.path.join(METRICS_DIR, "combined_molprobity.tsv")
    mp = pd.read_csv(mp_path, sep='\t')
    for col in ['clashscore', 'molprobity_score', 'rota_outliers']:
        mp[col] = pd.to_numeric(mp[col], errors='coerce')
    mp_blue = mp[mp['pipeline'] == pipe]
    mp_by_target = mp_blue.groupby(['target', 'source'])['molprobity_score'].mean().unstack()
    clash_by_target = mp_blue.groupby(['target', 'source'])['clashscore'].mean().unstack()

    # Merge into summary
    summary = pd.DataFrame(index=sorted(set(tm_by_target.index) | set(mp_by_target.index)))
    summary.index.name = 'target'

    # Add TM-score columns
    for src in SOURCE_ORDER:
        if src in tm_by_target.columns:
            summary[f'tm_{src}'] = tm_by_target[src]

    # Add MolProbity score columns
    for src in SOURCE_ORDER:
        if src in mp_by_target.columns:
            summary[f'mp_{src}'] = mp_by_target[src]

    # Add clashscore columns
    for src in SOURCE_ORDER:
        if src in clash_by_target.columns:
            summary[f'clash_{src}'] = clash_by_target[src]

    # Determine winners
    tm_cols = [f'tm_{s}' for s in SOURCE_ORDER if f'tm_{s}' in summary.columns]
    mp_cols = [f'mp_{s}' for s in SOURCE_ORDER if f'mp_{s}' in summary.columns]

    # Best TM (highest) - skip rows with all NaN
    tm_valid = summary[tm_cols].dropna(how='all')
    summary.loc[tm_valid.index, 'best_tm_source'] = tm_valid.idxmax(axis=1).str.replace('tm_', '')
    summary['best_tm_value'] = summary[tm_cols].max(axis=1)

    # Best MP (lowest) - skip rows with all NaN
    mp_valid = summary[mp_cols].dropna(how='all')
    summary.loc[mp_valid.index, 'best_mp_source'] = mp_valid.idxmin(axis=1).str.replace('mp_', '')
    summary['best_mp_value'] = summary[mp_cols].min(axis=1)

    # AMBER improvement metrics
    if 'tm_af_unrelaxed' in summary.columns and 'tm_amber_af' in summary.columns:
        summary['amber_af_tm_delta'] = summary['tm_amber_af'] - summary['tm_af_unrelaxed']
    if 'mp_af_unrelaxed' in summary.columns and 'mp_amber_af' in summary.columns:
        summary['amber_af_mp_delta'] = summary['mp_amber_af'] - summary['mp_af_unrelaxed']
    if 'tm_boltz' in summary.columns and 'tm_amber_boltz' in summary.columns:
        summary['amber_boltz_tm_delta'] = summary['tm_amber_boltz'] - summary['tm_boltz']
    if 'mp_boltz' in summary.columns and 'mp_amber_boltz' in summary.columns:
        summary['amber_boltz_mp_delta'] = summary['mp_amber_boltz'] - summary['mp_boltz']

    # Save
    outpath = os.path.join(OUTDIR, "target_summary.tsv")
    summary.to_csv(outpath, sep='\t', float_format='%.4f')
    print(f"Saved target summary: {len(summary)} targets")

    # Print interesting stats
    print(f"\n=== Summary Statistics ===")
    print(f"Targets with TM data: {summary[tm_cols].notna().any(axis=1).sum()}")
    print(f"Targets with MP data: {summary[mp_cols].notna().any(axis=1).sum()}")
    print(f"Targets with both: {(summary[tm_cols].notna().any(axis=1) & summary[mp_cols].notna().any(axis=1)).sum()}")

    # Best source distribution
    if 'best_tm_source' in summary.columns:
        print(f"\nBest TM-score source distribution:")
        counts = summary['best_tm_source'].value_counts()
        for src, count in counts.items():
            print(f"  {src}: {count} ({count/len(summary)*100:.1f}%)")

    if 'best_mp_source' in summary.columns:
        best_mp = summary['best_mp_source'].dropna().value_counts()
        print(f"\nBest MolProbity source distribution:")
        for src, count in best_mp.items():
            print(f"  {src}: {count} ({count/best_mp.sum()*100:.1f}%)")

    # AMBER improvement consistency
    if 'amber_af_mp_delta' in summary.columns:
        delta = summary['amber_af_mp_delta'].dropna()
        improved = (delta < 0).sum()
        worsened = (delta > 0).sum()
        print(f"\nAMBER(AF) MolProbity: {improved} improved, {worsened} worsened ({improved/(improved+worsened)*100:.1f}% improvement rate)")

    if 'amber_boltz_mp_delta' in summary.columns:
        delta = summary['amber_boltz_mp_delta'].dropna()
        improved = (delta < 0).sum()
        worsened = (delta > 0).sum()
        print(f"AMBER(Boltz) MolProbity: {improved} improved, {worsened} worsened ({improved/(improved+worsened)*100:.1f}% improvement rate)")

    # Targets where AMBER makes the biggest difference
    if 'amber_af_mp_delta' in summary.columns:
        biggest = summary.nlargest(5, 'amber_af_mp_delta')
        print(f"\nTargets where AMBER(AF) WORSENS MolProbity most:")
        for _, row in biggest.iterrows():
            print(f"  {row.name}: Δ={row['amber_af_mp_delta']:+.2f}")

        smallest = summary.nsmallest(5, 'amber_af_mp_delta')
        print(f"\nTargets where AMBER(AF) IMPROVES MolProbity most:")
        for _, row in smallest.iterrows():
            print(f"  {row.name}: Δ={row['amber_af_mp_delta']:+.2f}")


if __name__ == '__main__':
    main()
