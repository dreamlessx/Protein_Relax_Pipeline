#!/usr/bin/env python3
"""
statistical_analysis.py - Red Analysis Pipeline

Comprehensive statistical analysis combining TM-score and MolProbity data.
This is the main analysis script for the paper's Results section.

ANALYSES:
  1. Pre-Rosetta: Does AMBER improve MolProbity while preserving TM-score?
  2. Rosetta effect: How much does each protocol change TM-score?
  3. Protocol ranking: Which Rosetta protocol is best overall?
  4. Input type interaction: Does the best protocol depend on input type?
  5. Convergence: Are 5 replicates sufficient?
  6. Blue-Green agreement: Full reproducibility analysis
  7. Outlier analysis: Which targets behave unusually?

Statistical methods:
  - Wilcoxon signed-rank (paired, non-parametric)
  - Friedman test (repeated measures across protocols)
  - Cliff's delta (non-parametric effect size)
  - Pearson/Spearman correlation (Blue-Green agreement)
  - Bonferroni correction for multiple comparisons

COMMENT: This is where the paper lives or dies. Good science is about
asking the right questions and answering them honestly. The data tells
the story - our job is to present it fairly.
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
from itertools import combinations

METRICS_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"
OUTDIR = "/data/p_csb_meiler/agarwm5/red_analysis/tables"


def cliffs_delta(x, y):
    """Compute Cliff's delta effect size (non-parametric).
    Returns value in [-1, 1]. Thresholds: |d| < 0.147 negligible,
    < 0.33 small, < 0.474 medium, else large.
    """
    n_x, n_y = len(x), len(y)
    if n_x == 0 or n_y == 0:
        return float('nan')
    more = sum(1 for xi in x for yi in y if xi > yi)
    less = sum(1 for xi in x for yi in y if xi < yi)
    return (more - less) / (n_x * n_y)


def effect_size_label(d):
    """Interpret Cliff's delta magnitude."""
    d = abs(d)
    if d < 0.147:
        return "negligible"
    elif d < 0.33:
        return "small"
    elif d < 0.474:
        return "medium"
    else:
        return "large"


def load_all_data():
    """Load all available metric datasets."""
    data = {}

    # Pre-Rosetta TM-score
    path = os.path.join(METRICS_DIR, "combined_tmscore.tsv")
    if os.path.exists(path):
        df = pd.read_csv(path, sep='\t')
        for col in ['rmsd', 'tmscore', 'gdtts', 'gdtha']:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        data['pre_tm'] = df
        print(f"  Pre-Rosetta TM: {len(df)} rows, {df['target'].nunique()} targets")

    # Rosetta TM-score
    path = os.path.join(METRICS_DIR, "combined_rosetta_tmscore.tsv")
    if os.path.exists(path):
        df = pd.read_csv(path, sep='\t')
        for col in ['rmsd', 'tmscore', 'gdtts', 'gdtha']:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        data['ros_tm'] = df
        print(f"  Rosetta TM: {len(df)} rows, {df['target'].nunique()} targets")

    # MolProbity
    path = os.path.join(METRICS_DIR, "combined_molprobity.tsv")
    if os.path.exists(path):
        df = pd.read_csv(path, sep='\t')
        for col in ['clashscore', 'rama_outliers', 'rama_favored', 'rota_outliers',
                     'molprobity_score', 'cbeta_outliers', 'rms_bonds', 'rms_angles']:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        data['molprobity'] = df
        print(f"  MolProbity: {len(df)} rows, {df['target'].nunique()} targets")

    return data


def analysis_1_amber_dual_effect(data):
    """AMBER effect: TM-score unchanged + MolProbity improved = success story."""
    print("\n" + "="*80)
    print("ANALYSIS 1: AMBER DUAL EFFECT (TM-score + MolProbity)")
    print("="*80)

    results = []

    if 'pre_tm' not in data:
        print("  SKIP: No pre-Rosetta TM data")
        return results

    tm = data['pre_tm']

    # TM-score effect
    print("\n--- TM-score effect of AMBER ---")
    pairs = [
        ('af_unrelaxed', 'amber_af', 'AF → AMBER(AF)'),
        ('boltz', 'amber_boltz', 'Boltz → AMBER(Boltz)'),
        ('af_unrelaxed', 'af_relaxed', 'AF unrelaxed → AF AMBER-relaxed'),
    ]

    for src_before, src_after, label in pairs:
        for pipe in ['blue', 'green']:
            before = tm[(tm['pipeline'] == pipe) & (tm['source'] == src_before)]
            after = tm[(tm['pipeline'] == pipe) & (tm['source'] == src_after)]

            b_mean = before.groupby('target')['tmscore'].mean()
            a_mean = after.groupby('target')['tmscore'].mean()
            common = b_mean.index.intersection(a_mean.index)

            if len(common) < 10:
                continue

            bv = b_mean.loc[common].values
            av = a_mean.loc[common].values
            mask = np.isfinite(bv) & np.isfinite(av)

            if mask.sum() < 10:
                continue

            stat, pval = stats.wilcoxon(bv[mask], av[mask])
            delta = np.mean(av[mask]) - np.mean(bv[mask])
            cd = cliffs_delta(av[mask], bv[mask])

            r = {
                'analysis': 'AMBER_TM', 'pipeline': pipe, 'comparison': label,
                'metric': 'tmscore', 'n': mask.sum(),
                'before_mean': np.mean(bv[mask]), 'after_mean': np.mean(av[mask]),
                'delta': delta, 'p_value': pval, 'cliffs_delta': cd,
                'effect_size': effect_size_label(cd)
            }
            results.append(r)
            print(f"  {pipe} {label}: Δ={delta:+.4f}, p={pval:.2e}, d={cd:+.3f} ({effect_size_label(cd)}), n={mask.sum()}")

    # MolProbity effect
    if 'molprobity' in data:
        mp = data['molprobity']
        print("\n--- MolProbity effect of AMBER ---")

        mp_metrics = ['clashscore', 'molprobity_score', 'rota_outliers', 'rms_bonds']
        mp_pairs = [
            ('af_unrelaxed', 'amber_af', 'AF → AMBER(AF)'),
            ('boltz', 'amber_boltz', 'Boltz → AMBER(Boltz)'),
        ]

        for src_before, src_after, label in mp_pairs:
            for pipe in ['blue', 'green']:
                before = mp[(mp['pipeline'] == pipe) & (mp['source'] == src_before)]
                after = mp[(mp['pipeline'] == pipe) & (mp['source'] == src_after)]

                b_mean = before.groupby('target')[mp_metrics].mean()
                a_mean = after.groupby('target')[mp_metrics].mean()
                common = b_mean.index.intersection(a_mean.index)

                if len(common) < 10:
                    continue

                for metric in mp_metrics:
                    bv = b_mean.loc[common, metric].values
                    av = a_mean.loc[common, metric].values
                    mask = np.isfinite(bv) & np.isfinite(av)

                    if mask.sum() < 10:
                        continue

                    stat, pval = stats.wilcoxon(bv[mask], av[mask])
                    delta = np.mean(av[mask]) - np.mean(bv[mask])
                    cd = cliffs_delta(av[mask], bv[mask])

                    r = {
                        'analysis': 'AMBER_MP', 'pipeline': pipe, 'comparison': label,
                        'metric': metric, 'n': mask.sum(),
                        'before_mean': np.mean(bv[mask]), 'after_mean': np.mean(av[mask]),
                        'delta': delta, 'p_value': pval, 'cliffs_delta': cd,
                        'effect_size': effect_size_label(cd)
                    }
                    results.append(r)
                    print(f"  {pipe} {label} [{metric}]: Δ={delta:+.4f}, p={pval:.2e}, d={cd:+.3f} ({effect_size_label(cd)})")

    return results


def analysis_2_rosetta_effect(data):
    """Does Rosetta improve or degrade structures?"""
    print("\n" + "="*80)
    print("ANALYSIS 2: ROSETTA EFFECT ON TM-SCORE")
    print("="*80)

    results = []

    if 'pre_tm' not in data or 'ros_tm' not in data:
        print("  SKIP: Missing pre or Rosetta TM data")
        return results

    pre = data['pre_tm']
    ros = data['ros_tm']

    # Per-target mean for pre-Rosetta
    pre_by_target = pre.groupby(['pipeline', 'source', 'target'])['tmscore'].mean().reset_index()

    # Per-target mean for Rosetta (average across protocols and reps)
    ros_by_target = ros.groupby(['pipeline', 'source', 'target'])['tmscore'].mean().reset_index()

    for pipe in ['blue', 'green']:
        print(f"\n--- {pipe.upper()} ---")
        for src in ['af_relaxed', 'af_unrelaxed', 'boltz']:
            pre_src = pre_by_target[(pre_by_target['pipeline'] == pipe) &
                                    (pre_by_target['source'] == src)].set_index('target')['tmscore']
            ros_src = ros_by_target[(ros_by_target['pipeline'] == pipe) &
                                    (ros_by_target['source'] == src)].set_index('target')['tmscore']
            common = pre_src.index.intersection(ros_src.index)

            if len(common) < 10:
                continue

            pv = pre_src.loc[common].values
            rv = ros_src.loc[common].values
            mask = np.isfinite(pv) & np.isfinite(rv)

            if mask.sum() < 10:
                continue

            stat, pval = stats.wilcoxon(pv[mask], rv[mask])
            delta = np.mean(rv[mask]) - np.mean(pv[mask])
            cd = cliffs_delta(rv[mask], pv[mask])
            improved = np.sum(rv[mask] > pv[mask])
            degraded = np.sum(rv[mask] < pv[mask])

            r = {
                'analysis': 'ROSETTA_EFFECT', 'pipeline': pipe, 'source': src,
                'metric': 'tmscore', 'n': mask.sum(),
                'pre_mean': np.mean(pv[mask]), 'post_mean': np.mean(rv[mask]),
                'delta': delta, 'p_value': pval, 'cliffs_delta': cd,
                'effect_size': effect_size_label(cd),
                'improved': improved, 'degraded': degraded
            }
            results.append(r)
            print(f"  {src}: pre={np.mean(pv[mask]):.4f} → post={np.mean(rv[mask]):.4f}, "
                  f"Δ={delta:+.4f}, p={pval:.2e}, d={cd:+.3f} ({effect_size_label(cd)}), "
                  f"↑{improved}/↓{degraded}")

    return results


def analysis_3_protocol_ranking(data):
    """Friedman test: are protocol differences significant?"""
    print("\n" + "="*80)
    print("ANALYSIS 3: PROTOCOL RANKING (Friedman test)")
    print("="*80)

    results = []

    if 'ros_tm' not in data:
        print("  SKIP: No Rosetta TM data")
        return results

    ros = data['ros_tm']
    protocols = sorted(ros['protocol'].unique())

    for pipe in ['blue', 'green']:
        pipe_data = ros[ros['pipeline'] == pipe]
        # Per-target, per-protocol mean TM-score
        per_target = pipe_data.groupby(['target', 'protocol'])['tmscore'].mean().reset_index()
        pivot = per_target.pivot(index='target', columns='protocol', values='tmscore').dropna()

        if len(pivot) < 10 or len(pivot.columns) < 3:
            continue

        # Friedman test
        groups = [pivot[p].values for p in protocols if p in pivot.columns]
        if len(groups) >= 3:
            stat, pval = stats.friedmanchisquare(*groups)
            print(f"\n  {pipe.upper()}: Friedman χ² = {stat:.2f}, p = {pval:.2e}, n = {len(pivot)}")
            results.append({
                'analysis': 'FRIEDMAN', 'pipeline': pipe, 'stat': stat,
                'p_value': pval, 'n': len(pivot)
            })

            # Pairwise Wilcoxon with Bonferroni correction
            n_pairs = len(list(combinations(protocols, 2)))
            print(f"  Pairwise Wilcoxon ({n_pairs} comparisons, Bonferroni α = {0.05/n_pairs:.4f}):")

            for p1, p2 in combinations(protocols, 2):
                if p1 not in pivot.columns or p2 not in pivot.columns:
                    continue
                v1 = pivot[p1].values
                v2 = pivot[p2].values
                stat_w, pval_w = stats.wilcoxon(v1, v2)
                pval_adj = min(1.0, pval_w * n_pairs)
                delta = np.mean(v1) - np.mean(v2)
                sig = "*" if pval_adj < 0.05 else ""
                print(f"    {p1} vs {p2}: Δ={delta:+.4f}, p_adj={pval_adj:.3f} {sig}")

    return results


def analysis_4_convergence(data):
    """Are 5 replicates sufficient? Report per-target replicate variance."""
    print("\n" + "="*80)
    print("ANALYSIS 4: REPLICATE CONVERGENCE")
    print("="*80)

    if 'ros_tm' not in data:
        print("  SKIP: No Rosetta TM data")
        return

    ros = data['ros_tm']
    rep_std = ros.groupby(['pipeline', 'source', 'target', 'protocol'])[
        'tmscore'].std().reset_index()
    rep_std.columns = list(rep_std.columns[:-1]) + ['rep_std']
    rep_std = rep_std.dropna(subset=['rep_std'])

    for pipe in ['blue', 'green']:
        pipe_std = rep_std[rep_std['pipeline'] == pipe]
        pct_below_001 = 100 * (pipe_std['rep_std'] < 0.01).mean()
        pct_below_005 = 100 * (pipe_std['rep_std'] < 0.005).mean()

        print(f"\n  {pipe.upper()}:")
        print(f"    Mean σ = {pipe_std['rep_std'].mean():.4f}")
        print(f"    Median σ = {pipe_std['rep_std'].median():.4f}")
        print(f"    Max σ = {pipe_std['rep_std'].max():.4f}")
        print(f"    {pct_below_001:.1f}% converged (σ < 0.01)")
        print(f"    {pct_below_005:.1f}% tightly converged (σ < 0.005)")

        # Flag non-converging targets
        high_var = pipe_std[pipe_std['rep_std'] > 0.02]
        if len(high_var) > 0:
            print(f"    WARNING: {len(high_var)} target/protocol combos with σ > 0.02:")
            for _, row in high_var.nlargest(5, 'rep_std').iterrows():
                print(f"      {row['target']} / {row['source']} / {row['protocol']}: σ = {row['rep_std']:.4f}")


def analysis_5_blue_green_agreement(data):
    """Comprehensive Blue-Green reproducibility."""
    print("\n" + "="*80)
    print("ANALYSIS 5: BLUE-GREEN REPRODUCIBILITY")
    print("="*80)

    for name, key, metric in [
        ('Pre-Rosetta TM', 'pre_tm', 'tmscore'),
        ('Rosetta TM', 'ros_tm', 'tmscore'),
    ]:
        if key not in data:
            continue

        df = data[key]
        blue = df[df['pipeline'] == 'blue'].groupby(['source', 'target'])[metric].mean().reset_index()
        green = df[df['pipeline'] == 'green'].groupby(['source', 'target'])[metric].mean().reset_index()

        merged = blue.merge(green, on=['source', 'target'], suffixes=('_blue', '_green'))
        merged = merged.dropna()

        if len(merged) < 10:
            continue

        r_pearson, p_pearson = stats.pearsonr(merged[f'{metric}_blue'], merged[f'{metric}_green'])
        r_spearman, p_spearman = stats.spearmanr(merged[f'{metric}_blue'], merged[f'{metric}_green'])
        mae = np.mean(np.abs(merged[f'{metric}_blue'] - merged[f'{metric}_green']))

        print(f"\n  {name}:")
        print(f"    Pearson r = {r_pearson:.4f} (p = {p_pearson:.2e})")
        print(f"    Spearman ρ = {r_spearman:.4f} (p = {p_spearman:.2e})")
        print(f"    MAE = {mae:.4f}")
        print(f"    n = {len(merged)} (source × target pairs)")

    # MolProbity agreement
    if 'molprobity' in data:
        mp = data['molprobity']
        for metric in ['clashscore', 'molprobity_score']:
            blue = mp[mp['pipeline'] == 'blue'].groupby(['source', 'target'])[metric].mean().reset_index()
            green = mp[mp['pipeline'] == 'green'].groupby(['source', 'target'])[metric].mean().reset_index()
            merged = blue.merge(green, on=['source', 'target'], suffixes=('_blue', '_green')).dropna()

            if len(merged) < 10:
                continue

            r, p = stats.pearsonr(merged[f'{metric}_blue'], merged[f'{metric}_green'])
            print(f"\n  MolProbity {metric}:")
            print(f"    Pearson r = {r:.4f} (p = {p:.2e}), n = {len(merged)}")


def analysis_6_outliers(data):
    """Identify and characterize outlier targets."""
    print("\n" + "="*80)
    print("ANALYSIS 6: OUTLIER TARGETS")
    print("="*80)

    if 'pre_tm' not in data:
        return

    tm = data['pre_tm']
    # Per-target best TM-score across all models
    best_tm = tm.groupby(['target', 'pipeline'])['tmscore'].max().reset_index()
    worst = best_tm.groupby('target')['tmscore'].min().reset_index()
    worst = worst.sort_values('tmscore')

    print(f"\n  Targets with worst best-TM-score (bottom 20):")
    for _, row in worst.head(20).iterrows():
        print(f"    {row['target']}: best TM = {row['tmscore']:.3f}")

    # Catastrophic failures (TM < 0.5)
    catastrophic = worst[worst['tmscore'] < 0.5]
    print(f"\n  Catastrophic failures (TM < 0.5): {len(catastrophic)} targets")
    for _, row in catastrophic.iterrows():
        print(f"    {row['target']}: {row['tmscore']:.3f}")

    # Targets where Boltz fails but AF succeeds (or vice versa)
    for pipe in ['blue']:
        pipe_data = tm[tm['pipeline'] == pipe]
        af = pipe_data[pipe_data['source'] == 'af_relaxed'].groupby('target')['tmscore'].mean()
        boltz = pipe_data[pipe_data['source'] == 'boltz'].groupby('target')['tmscore'].mean()
        common = af.index.intersection(boltz.index)

        if len(common) < 10:
            continue

        diff = af.loc[common] - boltz.loc[common]
        af_wins = diff[diff > 0.1]
        boltz_wins = diff[diff < -0.1]

        print(f"\n  AF >> Boltz (ΔTM > 0.1): {len(af_wins)} targets")
        for t in af_wins.nlargest(5).index:
            print(f"    {t}: AF={af.loc[t]:.3f}, Boltz={boltz.loc[t]:.3f}")

        print(f"\n  Boltz >> AF (ΔTM > 0.1): {len(boltz_wins)} targets")
        for t in boltz_wins.nsmallest(5).index:
            print(f"    {t}: AF={af.loc[t]:.3f}, Boltz={boltz.loc[t]:.3f}")


def main():
    print("Loading all metric data...")
    data = load_all_data()

    if not data:
        print("ERROR: No data found")
        sys.exit(1)

    all_results = []

    r = analysis_1_amber_dual_effect(data)
    all_results.extend(r)

    r = analysis_2_rosetta_effect(data)
    all_results.extend(r)

    analysis_3_protocol_ranking(data)
    analysis_4_convergence(data)
    analysis_5_blue_green_agreement(data)
    analysis_6_outliers(data)

    # Save all pairwise results
    if all_results:
        df = pd.DataFrame(all_results)
        outpath = os.path.join(OUTDIR, "statistical_tests.tsv")
        df.to_csv(outpath, sep='\t', index=False)
        print(f"\n\nAll statistical tests -> {outpath}")

    print("\n" + "="*80)
    print("STATISTICAL ANALYSIS COMPLETE")
    print("="*80)


if __name__ == '__main__':
    main()
