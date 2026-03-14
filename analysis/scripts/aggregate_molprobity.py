#!/usr/bin/env python3
"""
aggregate_molprobity.py — Red Analysis Pipeline

Merges per-target MolProbity TSV files into a combined dataset and
computes summary statistics.

COMMENT: This is the companion to aggregate_tmscore.py. Where TM-score
told us "AMBER doesn't change global fold", MolProbity tells us whether
AMBER fixes local geometry problems. The hypothesis from Phase 1:

  - AMBER should dramatically improve clashscore (steric clashes)
  - AMBER should improve rotamer outliers
  - Ramachandran should be relatively unaffected (backbone geometry)
  - RMS bonds/angles should be near-ideal after AMBER

Initial 1A2K test data already supports this:
  Boltz clashscore 8.37 → amber_boltz 0.70  (6x reduction!)
  af_unrelaxed clashscore 24.35 → amber_af 4.82

This script will tell us if the pattern holds across all 257 targets.
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
from scipy import stats

METRICS_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"
OUTDIR = "/data/p_csb_meiler/agarwm5/red_analysis/tables"

def main():
    # Collect all MolProbity TSVs
    files = sorted(glob.glob(os.path.join(METRICS_DIR, "molprobity_*.tsv")))
    print(f"Found {len(files)} MolProbity files")

    if not files:
        print("ERROR: No MolProbity files found")
        sys.exit(1)

    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f, sep='\t')
            if len(df) > 0:
                dfs.append(df)
        except Exception as e:
            print(f"  WARN: Skipping {f}: {e}", file=sys.stderr)

    data = pd.concat(dfs, ignore_index=True)

    # Convert numeric columns
    for col in ['clashscore', 'rama_outliers', 'rama_favored', 'rota_outliers',
                'molprobity_score', 'cbeta_outliers', 'rms_bonds', 'rms_angles']:
        data[col] = pd.to_numeric(data[col], errors='coerce')

    # Drop duplicate rows (from re-runs)
    data = data.drop_duplicates(subset=['target', 'pipeline', 'source', 'model_idx'])

    # Save combined
    combined_path = os.path.join(METRICS_DIR, "combined_molprobity.tsv")
    data.to_csv(combined_path, sep='\t', index=False)
    print(f"\nCombined: {len(data)} rows -> {combined_path}")

    # Unique targets and pipelines
    print(f"Targets: {data['target'].nunique()}")
    print(f"Pipelines: {data['pipeline'].unique().tolist()}")
    print(f"Sources: {data['source'].unique().tolist()}")

    # === Summary by source ===
    # Average across 5 models per target, then across targets
    pred_data = data[data['source'] != 'crystal'].copy()
    metrics = ['clashscore', 'rama_outliers', 'rama_favored', 'rota_outliers',
               'molprobity_score', 'cbeta_outliers', 'rms_bonds', 'rms_angles']

    # Per-target mean (average over 5 models)
    per_target = pred_data.groupby(['pipeline', 'source', 'target'])[metrics].mean().reset_index()

    # Summary across targets
    summary = per_target.groupby(['pipeline', 'source'])[metrics].agg(['mean', 'median', 'std']).reset_index()
    summary.columns = ['_'.join(col).rstrip('_') for col in summary.columns]

    summary_path = os.path.join(OUTDIR, "molprobity_summary.tsv")
    summary.to_csv(summary_path, sep='\t', index=False)
    print(f"\nSummary -> {summary_path}")

    # === Print key comparisons ===
    print("\n" + "="*80)
    print("KEY MOLPROBITY COMPARISONS (mean across targets)")
    print("="*80)

    for pipe in ['blue', 'green']:
        pipe_data = per_target[per_target['pipeline'] == pipe]
        print(f"\n--- {pipe.upper()} ---")
        print(f"{'Source':<15} {'Clashscore':>10} {'Rama Out%':>10} {'Rama Fav%':>10} "
              f"{'Rota Out%':>10} {'MP Score':>10} {'RMS Bond':>10} {'RMS Angle':>10}")
        print("-"*95)
        for src in ['af_relaxed', 'af_unrelaxed', 'boltz', 'amber_af', 'amber_boltz']:
            s = pipe_data[pipe_data['source'] == src]
            if len(s) == 0:
                continue
            print(f"{src:<15} {s['clashscore'].mean():>10.2f} {s['rama_outliers'].mean():>10.2f} "
                  f"{s['rama_favored'].mean():>10.2f} {s['rota_outliers'].mean():>10.2f} "
                  f"{s['molprobity_score'].mean():>10.2f} {s['rms_bonds'].mean():>10.4f} "
                  f"{s['rms_angles'].mean():>10.4f}")

        # Crystal baseline (separate — only 1 model per target)
        crystal = data[(data['pipeline'] == pipe) & (data['source'] == 'crystal')]
        if len(crystal) > 0:
            print(f"{'crystal':<15} {crystal['clashscore'].mean():>10.2f} {crystal['rama_outliers'].mean():>10.2f} "
                  f"{crystal['rama_favored'].mean():>10.2f} {crystal['rota_outliers'].mean():>10.2f} "
                  f"{crystal['molprobity_score'].mean():>10.2f} {crystal['rms_bonds'].mean():>10.4f} "
                  f"{crystal['rms_angles'].mean():>10.4f}")

    # === Statistical tests: AMBER effect on MolProbity ===
    print("\n" + "="*80)
    print("STATISTICAL TESTS: AMBER EFFECT ON MOLPROBITY")
    print("="*80)

    stat_results = []
    for pipe in ['blue', 'green']:
        pipe_data = per_target[per_target['pipeline'] == pipe]

        # AMBER effect on AF: af_unrelaxed → amber_af
        af_unrel = pipe_data[pipe_data['source'] == 'af_unrelaxed'].set_index('target')
        amber_af = pipe_data[pipe_data['source'] == 'amber_af'].set_index('target')
        common = af_unrel.index.intersection(amber_af.index)

        if len(common) > 10:
            for metric in ['clashscore', 'molprobity_score', 'rota_outliers', 'rms_bonds']:
                x = af_unrel.loc[common, metric].values
                y = amber_af.loc[common, metric].values
                mask = np.isfinite(x) & np.isfinite(y)
                if mask.sum() > 10:
                    stat, pval = stats.wilcoxon(x[mask], y[mask], alternative='two-sided')
                    diff = np.mean(y[mask]) - np.mean(x[mask])
                    print(f"  {pipe} af_unrelaxed→amber_af {metric}: Δ={diff:+.4f}, p={pval:.2e}, n={mask.sum()}")
                    stat_results.append({
                        'pipeline': pipe, 'comparison': 'af_unrelaxed→amber_af',
                        'metric': metric, 'delta': diff, 'p_value': pval, 'n': mask.sum()
                    })

        # AMBER effect on Boltz: boltz → amber_boltz
        boltz = pipe_data[pipe_data['source'] == 'boltz'].set_index('target')
        amber_boltz = pipe_data[pipe_data['source'] == 'amber_boltz'].set_index('target')
        common = boltz.index.intersection(amber_boltz.index)

        if len(common) > 10:
            for metric in ['clashscore', 'molprobity_score', 'rota_outliers', 'rms_bonds']:
                x = boltz.loc[common, metric].values
                y = amber_boltz.loc[common, metric].values
                mask = np.isfinite(x) & np.isfinite(y)
                if mask.sum() > 10:
                    stat, pval = stats.wilcoxon(x[mask], y[mask], alternative='two-sided')
                    diff = np.mean(y[mask]) - np.mean(x[mask])
                    print(f"  {pipe} boltz→amber_boltz {metric}: Δ={diff:+.4f}, p={pval:.2e}, n={mask.sum()}")
                    stat_results.append({
                        'pipeline': pipe, 'comparison': 'boltz→amber_boltz',
                        'metric': metric, 'delta': diff, 'p_value': pval, 'n': mask.sum()
                    })

    if stat_results:
        stat_df = pd.DataFrame(stat_results)
        stat_path = os.path.join(OUTDIR, "molprobity_amber_effect.tsv")
        stat_df.to_csv(stat_path, sep='\t', index=False)
        print(f"\nStatistical tests -> {stat_path}")

    # === Blue-Green agreement ===
    print("\n" + "="*80)
    print("BLUE-GREEN AGREEMENT ON MOLPROBITY")
    print("="*80)

    blue_data = per_target[per_target['pipeline'] == 'blue'].copy()
    green_data = per_target[per_target['pipeline'] == 'green'].copy()

    for src in ['af_relaxed', 'af_unrelaxed', 'boltz', 'amber_af', 'amber_boltz']:
        b = blue_data[blue_data['source'] == src].set_index('target')
        g = green_data[green_data['source'] == src].set_index('target')
        common = b.index.intersection(g.index)
        if len(common) > 10:
            for metric in ['clashscore', 'molprobity_score']:
                bv = b.loc[common, metric].values
                gv = g.loc[common, metric].values
                mask = np.isfinite(bv) & np.isfinite(gv)
                if mask.sum() > 10:
                    r, p = stats.pearsonr(bv[mask], gv[mask])
                    print(f"  {src} {metric}: r={r:.3f}, p={p:.2e}, n={mask.sum()}")

    print("\nDone!")


if __name__ == '__main__':
    main()
