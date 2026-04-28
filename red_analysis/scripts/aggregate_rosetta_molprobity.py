#!/usr/bin/env python3
"""
aggregate_rosetta_molprobity.py - Red Analysis Pipeline

Phase 4 aggregation: merges per-target Rosetta MolProbity TSVs and
compares against pre-Rosetta MolProbity (Phase 3).

THE KEY QUESTION: Does Rosetta improve MolProbity enough to justify
the TM-score degradation documented in Phase 2?

We already know:
  - Rosetta degrades TM-score by ~0.02 (medium effect, all 20 targets, Phase 2)
  - AMBER improves clashscore by 21 points (d=-0.99, 257 targets, Phase 3)
  - AMBER has zero TM-score effect (d=-0.01, Phase 1)

Three possible outcomes for Rosetta MolProbity:
  A) Rosetta improves MolProbity MORE than AMBER → Rosetta adds value
     despite TM cost (useful for downstream applications)
  B) Rosetta improves MolProbity SAME as AMBER → Rosetta is pointless
     (AMBER gives the same benefit with zero TM cost)
  C) Rosetta improves MolProbity LESS than AMBER → Rosetta is harmful
     (worse on BOTH metrics)

The 1GCQ pilot data suggests outcome (A): Rosetta achieves clashscore
~0.1-0.5, compared to AMBER's ~1.4. But we need multi-target data.
"""

import os
import sys
import glob
import math
import numpy as np
import pandas as pd
from scipy import stats

METRICS_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"
OUTDIR = "/data/p_csb_meiler/agarwm5/red_analysis/tables"


def classify_src_type(src_name):
    """Map Rosetta directory names back to canonical source types.

    Blue Rosetta dirs: af_relaxed_ranked_0, af_unrelaxed_ranked_0,
                       boltz_boltz_input_model_0, amber_af_relaxed,
                       amber_boltz_relaxed, crystal_{TARGET}
    Green Rosetta dirs: similar but amber_af_ranked_0, amber_boltz_model_0
    """
    s = src_name.lower()
    if s.startswith('crystal'):
        return 'crystal'
    elif s.startswith('amber_af'):
        return 'amber_af'
    elif s.startswith('amber_boltz'):
        return 'amber_boltz'
    elif s.startswith('af_relaxed'):
        return 'af_relaxed'
    elif s.startswith('af_unrelaxed'):
        return 'af_unrelaxed'
    elif s.startswith('boltz'):
        return 'boltz'
    else:
        return 'unknown'


def cliffs_delta(x, y):
    """Cliff's delta effect size (non-parametric)."""
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return 0.0
    dom = 0
    for xi in x:
        for yi in y:
            if xi > yi:
                dom += 1
            elif xi < yi:
                dom -= 1
    return dom / (n1 * n2)


def interpret_delta(d):
    """Interpret Cliff's delta magnitude."""
    ad = abs(d)
    if ad < 0.147:
        return "negligible"
    elif ad < 0.33:
        return "small"
    elif ad < 0.474:
        return "medium"
    else:
        return "large"


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    # Collect all Rosetta MolProbity TSVs
    files = sorted(glob.glob(os.path.join(METRICS_DIR, "rosetta_molprobity_*.tsv")))
    print(f"Found {len(files)} Rosetta MolProbity files")

    if not files:
        print("ERROR: No Rosetta MolProbity files found")
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
    mp_metrics = ['clashscore', 'molprobity_score', 'rama_outliers', 'rama_favored',
                  'rota_outliers', 'rms_bonds', 'rms_angles', 'cbeta_outliers']
    for col in mp_metrics:
        if col in data.columns:
            data[col] = pd.to_numeric(data[col], errors='coerce')

    # Classify source types
    data['source'] = data['src_type'].apply(classify_src_type)

    # Drop duplicates
    data = data.drop_duplicates(subset=['target', 'pipeline', 'src_type', 'protocol', 'rep'])

    # Save combined
    combined_path = os.path.join(METRICS_DIR, "combined_rosetta_molprobity.tsv")
    data.to_csv(combined_path, sep='\t', index=False)
    print(f"\nCombined: {len(data)} rows -> {combined_path}")
    print(f"Targets: {data['target'].nunique()}")
    print(f"Pipelines: {sorted(data['pipeline'].unique().tolist())}")
    print(f"Protocols: {sorted(data['protocol'].unique().tolist())}")
    print(f"Source types: {sorted(data['source'].unique().tolist())}")

    # === Analysis 1: Rosetta MolProbity by protocol ===
    print("\n" + "="*80)
    print("ROSETTA MolProbity BY PROTOCOL (mean across all sources)")
    print("="*80)

    key_metrics = ['clashscore', 'molprobity_score', 'rota_outliers', 'rama_favored']

    for pipe in ['blue', 'green']:
        pipe_data = data[data['pipeline'] == pipe]
        if len(pipe_data) == 0:
            continue
        # Average across reps and models within source, then across targets
        per_target = pipe_data.groupby(['source', 'target', 'protocol'])[key_metrics].mean().reset_index()
        by_proto = per_target.groupby('protocol')[key_metrics].agg(['mean', 'std'])

        print(f"\n{pipe.upper()}: (n={pipe_data['target'].nunique()} targets)")
        print(f"{'Protocol':<20} {'Clashscore':>12} {'MP Score':>12} {'Rota Out%':>12} {'Rama Fav%':>12}")
        print("-"*70)
        for proto in sorted(by_proto.index):
            row = by_proto.loc[proto]
            c = row[('clashscore', 'mean')]
            m = row[('molprobity_score', 'mean')]
            r = row[('rota_outliers', 'mean')]
            rf = row[('rama_favored', 'mean')]
            print(f"{proto:<20} {c:>12.3f} {m:>12.3f} {r:>12.3f} {rf:>12.2f}")

    # === Analysis 2: Pre vs Post Rosetta MolProbity (THE KEY COMPARISON) ===
    print("\n" + "="*80)
    print("PRE vs POST ROSETTA MolProbity COMPARISON")
    print("="*80)

    pre_path = os.path.join(METRICS_DIR, "combined_molprobity.tsv")
    if os.path.exists(pre_path):
        pre_data = pd.read_csv(pre_path, sep='\t')
        for col in key_metrics:
            pre_data[col] = pd.to_numeric(pre_data[col], errors='coerce')

        for pipe in ['blue', 'green']:
            pipe_post = data[data['pipeline'] == pipe]
            pipe_pre = pre_data[pre_data['pipeline'] == pipe]
            if len(pipe_post) == 0:
                continue

            print(f"\n{pipe.upper()} pipeline:")

            for src in ['af_relaxed', 'af_unrelaxed', 'boltz', 'amber_af', 'amber_boltz', 'crystal']:
                # Pre-Rosetta per-target means
                pre_src = pipe_pre[pipe_pre['source'] == src]
                pre_by_target = pre_src.groupby('target')[key_metrics].mean()

                # Post-Rosetta per-target means (average across all protocols)
                post_src = pipe_post[pipe_post['source'] == src]
                if len(post_src) == 0:
                    continue
                post_by_target = post_src.groupby('target')[key_metrics].mean()

                common = pre_by_target.index.intersection(post_by_target.index)
                if len(common) < 3:
                    continue

                print(f"\n  {src} (n={len(common)}):")
                for metric in ['clashscore', 'molprobity_score', 'rota_outliers', 'rama_favored']:
                    pre_vals = pre_by_target.loc[common, metric].values
                    post_vals = post_by_target.loc[common, metric].values
                    mask = np.isfinite(pre_vals) & np.isfinite(post_vals)

                    if mask.sum() < 3:
                        continue

                    pv, pov = pre_vals[mask], post_vals[mask]
                    delta_mean = np.mean(pov) - np.mean(pv)

                    # Lower is better for clashscore, mp_score, rota_outliers
                    # Higher is better for rama_favored
                    if metric == 'rama_favored':
                        improved = np.sum(pov > pv)
                        degraded = np.sum(pov < pv)
                    else:
                        improved = np.sum(pov < pv)
                        degraded = np.sum(pov > pv)

                    # Wilcoxon signed-rank test
                    try:
                        stat, pval = stats.wilcoxon(pv, pov)
                    except:
                        pval = float('nan')

                    # Effect size
                    d = cliffs_delta(pov.tolist(), pv.tolist())

                    print(f"    {metric:<20} pre={np.mean(pv):>7.3f} -> post={np.mean(pov):>7.3f} "
                          f"(Δ={delta_mean:>+7.3f}, d={d:>+.3f} {interpret_delta(d)}, "
                          f"p={pval:.2e}, imp={improved}, deg={degraded})")

    # === Analysis 3: Best protocol for MolProbity ===
    print("\n" + "="*80)
    print("BEST PROTOCOL FOR MolProbity (per-target best)")
    print("="*80)

    blue_data = data[data['pipeline'] == 'blue']
    if len(blue_data) > 0:
        per_target_proto = blue_data.groupby(['target', 'protocol'])[
            ['clashscore', 'molprobity_score']].mean().reset_index()

        # Drop rows with all-NaN metrics to avoid idxmin errors
        valid_clash = per_target_proto.dropna(subset=['clashscore'])
        valid_mp = per_target_proto.dropna(subset=['molprobity_score'])

        best_clash = valid_clash.loc[
            valid_clash.groupby('target')['clashscore'].idxmin()
        ] if len(valid_clash) > 0 else pd.DataFrame()
        best_mp = valid_mp.loc[
            valid_mp.groupby('target')['molprobity_score'].idxmin()
        ] if len(valid_mp) > 0 else pd.DataFrame()

        print("\nBest clashscore protocol distribution:")
        for proto, count in best_clash['protocol'].value_counts().items():
            print(f"  {proto}: {count} ({100*count/len(best_clash):.0f}%)")

        print("\nBest MolProbity score protocol distribution:")
        for proto, count in best_mp['protocol'].value_counts().items():
            print(f"  {proto}: {count} ({100*count/len(best_mp):.0f}%)")

    # === Analysis 4: Rosetta vs AMBER comparison ===
    print("\n" + "="*80)
    print("ROSETTA vs AMBER COMPARISON (the critical question)")
    print("="*80)

    if os.path.exists(pre_path):
        pre_data = pd.read_csv(pre_path, sep='\t')
        for col in key_metrics:
            pre_data[col] = pd.to_numeric(pre_data[col], errors='coerce')

        pipe = 'blue'
        pipe_pre = pre_data[pre_data['pipeline'] == pipe]
        pipe_post = data[data['pipeline'] == pipe]

        if len(pipe_post) > 0:
            # AMBER(AF) MolProbity
            amber_af_mp = pipe_pre[pipe_pre['source'] == 'amber_af'].groupby('target')[key_metrics].mean()

            # Rosetta(AF relaxed, best protocol) MolProbity
            ros_af = pipe_post[pipe_post['source'] == 'af_relaxed']
            if len(ros_af) > 0:
                # Per-target best protocol
                ros_af_by_target = ros_af.groupby('target')[key_metrics].mean()

                common = amber_af_mp.index.intersection(ros_af_by_target.index)
                if len(common) >= 3:
                    print(f"\n  AMBER(AF) vs Rosetta(AF) (n={len(common)} targets):")
                    for metric in ['clashscore', 'molprobity_score', 'rota_outliers']:
                        amber_vals = amber_af_mp.loc[common, metric].values
                        ros_vals = ros_af_by_target.loc[common, metric].values
                        mask = np.isfinite(amber_vals) & np.isfinite(ros_vals)
                        if mask.sum() < 3:
                            continue
                        print(f"    {metric:<20} AMBER={np.mean(amber_vals[mask]):>7.3f} vs "
                              f"Rosetta={np.mean(ros_vals[mask]):>7.3f} "
                              f"(Δ={np.mean(ros_vals[mask])-np.mean(amber_vals[mask]):>+7.3f})")

            # AMBER(Boltz) vs Rosetta(Boltz)
            amber_boltz_mp = pipe_pre[pipe_pre['source'] == 'amber_boltz'].groupby('target')[key_metrics].mean()
            ros_boltz = pipe_post[pipe_post['source'] == 'boltz']
            if len(ros_boltz) > 0:
                ros_boltz_by_target = ros_boltz.groupby('target')[key_metrics].mean()
                common = amber_boltz_mp.index.intersection(ros_boltz_by_target.index)
                if len(common) >= 3:
                    print(f"\n  AMBER(Boltz) vs Rosetta(Boltz) (n={len(common)} targets):")
                    for metric in ['clashscore', 'molprobity_score', 'rota_outliers']:
                        amber_vals = amber_boltz_mp.loc[common, metric].values
                        ros_vals = ros_boltz_by_target.loc[common, metric].values
                        mask = np.isfinite(amber_vals) & np.isfinite(ros_vals)
                        if mask.sum() < 3:
                            continue
                        print(f"    {metric:<20} AMBER={np.mean(amber_vals[mask]):>7.3f} vs "
                              f"Rosetta={np.mean(ros_vals[mask]):>7.3f} "
                              f"(Δ={np.mean(ros_vals[mask])-np.mean(amber_vals[mask]):>+7.3f})")

    # === Analysis 5: Protocol × Source interaction on MolProbity ===
    print("\n" + "="*80)
    print("PROTOCOL × SOURCE INTERACTION (Blue, clashscore)")
    print("="*80)

    if len(blue_data) > 0:
        per_target_detail = blue_data.groupby(['source', 'target', 'protocol'])[
            'clashscore'].mean().reset_index()
        protocols = sorted(per_target_detail['protocol'].unique())
        sources = ['af_relaxed', 'af_unrelaxed', 'boltz', 'amber_af', 'amber_boltz', 'crystal']

        print(f"\n{'Source':<15}", end='')
        for p in protocols:
            print(f"  {p[:12]:>12}", end='')
        print()
        print("-"*(15 + 14*len(protocols)))

        for src in sources:
            src_data = per_target_detail[per_target_detail['source'] == src]
            if len(src_data) == 0:
                continue
            print(f"{src:<15}", end='')
            for proto in protocols:
                proto_data = src_data[src_data['protocol'] == proto]
                if len(proto_data) > 0:
                    print(f"  {proto_data['clashscore'].mean():>12.3f}", end='')
                else:
                    print(f"  {'NA':>12}", end='')
            print()

    # Save summary
    if len(data) > 0:
        summary = data.groupby(['pipeline', 'source', 'protocol'])[key_metrics].agg(
            ['mean', 'std', 'count']).reset_index()
        summary.to_csv(os.path.join(OUTDIR, "rosetta_molprobity_summary.tsv"), sep='\t', index=False)
        print(f"\n\nSummary saved to {OUTDIR}")

    print("\nDone!")


if __name__ == '__main__':
    main()
