#!/usr/bin/env python3
"""
aggregate_rosetta_tmscore.py — Red Analysis Pipeline

Merges per-target Rosetta TM-score TSVs and performs comprehensive analysis.

COMMENT: This is the main analysis script for the paper. The questions:
  1. Which Rosetta protocol is best (by TM-score to crystal)?
  2. Does Rosetta improve or degrade structures vs pre-Rosetta?
  3. Which input type benefits most from Rosetta?
  4. How do the 5 replicate runs converge?
  5. Are results consistent between Blue and Green?
  6. Is there an interaction between protocol and input type?

The pre-Rosetta analysis showed TM-scores of ~0.94-0.95. If Rosetta
improves things, post-Rosetta should be higher. If it degrades, lower.
Preliminary data (104 targets) shows degradation of ~0.02 — which is
expected since Rosetta optimizes energy, not similarity to crystal.
"""

import os
import sys
import glob
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


def main():
    # Collect all Rosetta TM-score TSVs
    files = sorted(glob.glob(os.path.join(METRICS_DIR, "rosetta_tmscore_*.tsv")))
    print(f"Found {len(files)} Rosetta TM-score files")

    if not files:
        print("ERROR: No Rosetta TM-score files found")
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
    for col in ['rmsd', 'tmscore', 'gdtts', 'gdtha', 'aligned_len', 'seq_len']:
        data[col] = pd.to_numeric(data[col], errors='coerce')

    # Drop duplicates
    data = data.drop_duplicates(subset=['target', 'pipeline', 'src_type', 'protocol', 'rep'])

    # Classify source types
    data['source'] = data['src_type'].apply(classify_src_type)

    # Save combined
    combined_path = os.path.join(METRICS_DIR, "combined_rosetta_tmscore.tsv")
    data.to_csv(combined_path, sep='\t', index=False)
    print(f"\nCombined: {len(data)} rows -> {combined_path}")
    print(f"Targets: {data['target'].nunique()}")
    print(f"Pipelines: {data['pipeline'].unique().tolist()}")
    print(f"Protocols: {data['protocol'].unique().tolist()}")
    print(f"Source types: {sorted(data['source'].unique().tolist())}")

    # === Analysis 1: Protocol comparison ===
    # Average across reps, then across models within source, then across targets
    per_rep = data.groupby(['pipeline', 'source', 'target', 'protocol', 'rep'])[
        ['tmscore', 'rmsd']].mean().reset_index()
    per_target_proto = per_rep.groupby(['pipeline', 'source', 'target', 'protocol'])[
        ['tmscore', 'rmsd']].mean().reset_index()

    print("\n" + "="*80)
    print("PROTOCOL COMPARISON (mean TM-score across all sources and targets)")
    print("="*80)

    for pipe in ['blue', 'green']:
        pipe_data = per_target_proto[per_target_proto['pipeline'] == pipe]
        by_proto = pipe_data.groupby('protocol')['tmscore'].agg(['mean', 'std', 'count'])
        by_proto = by_proto.sort_values('mean', ascending=False)
        print(f"\n{pipe.upper()}:")
        for proto, row in by_proto.iterrows():
            print(f"  {proto:<20} TM={row['mean']:.4f} ± {row['std']:.4f}  (n={int(row['count'])})")

    # === Analysis 2: Protocol × Source interaction ===
    print("\n" + "="*80)
    print("PROTOCOL × SOURCE INTERACTION (Blue pipeline)")
    print("="*80)

    blue_data = per_target_proto[per_target_proto['pipeline'] == 'blue']
    protocols = sorted(blue_data['protocol'].unique())
    sources = ['af_relaxed', 'af_unrelaxed', 'boltz', 'amber_af', 'amber_boltz', 'crystal']

    print(f"\n{'Source':<15}", end='')
    for p in protocols:
        print(f"  {p[:12]:>12}", end='')
    print()
    print("-"*(15 + 14*len(protocols)))

    for src in sources:
        src_data = blue_data[blue_data['source'] == src]
        if len(src_data) == 0:
            continue
        print(f"{src:<15}", end='')
        for proto in protocols:
            proto_data = src_data[src_data['protocol'] == proto]
            if len(proto_data) > 0:
                print(f"  {proto_data['tmscore'].mean():>12.4f}", end='')
            else:
                print(f"  {'NA':>12}", end='')
        print()

    # === Analysis 3: Rosetta vs Pre-Rosetta ===
    print("\n" + "="*80)
    print("ROSETTA vs PRE-ROSETTA COMPARISON")
    print("="*80)

    pre_rosetta_path = os.path.join(METRICS_DIR, "combined_tmscore.tsv")
    if os.path.exists(pre_rosetta_path):
        pre_data = pd.read_csv(pre_rosetta_path, sep='\t')
        pre_data['tmscore'] = pd.to_numeric(pre_data['tmscore'], errors='coerce')

        for pipe in ['blue', 'green']:
            print(f"\n{pipe.upper()}:")
            for src in ['af_relaxed', 'af_unrelaxed', 'boltz', 'amber_af', 'amber_boltz']:
                # Pre-Rosetta
                pre = pre_data[(pre_data['pipeline'] == pipe) & (pre_data['source'] == src)]
                pre_by_target = pre.groupby('target')['tmscore'].mean()

                # Post-Rosetta (best protocol per target)
                post = per_target_proto[(per_target_proto['pipeline'] == pipe) &
                                        (per_target_proto['source'] == src)]
                if len(post) == 0:
                    continue
                # Best single protocol: use mean across all protocols
                post_by_target = post.groupby('target')['tmscore'].mean()

                common = pre_by_target.index.intersection(post_by_target.index)
                if len(common) < 5:
                    continue

                pre_vals = pre_by_target.loc[common].values
                post_vals = post_by_target.loc[common].values
                mask = np.isfinite(pre_vals) & np.isfinite(post_vals)

                if mask.sum() > 10:
                    delta = np.mean(post_vals[mask]) - np.mean(pre_vals[mask])
                    stat, pval = stats.wilcoxon(pre_vals[mask], post_vals[mask])
                    improved = np.sum(post_vals[mask] > pre_vals[mask])
                    degraded = np.sum(post_vals[mask] < pre_vals[mask])
                    tied = np.sum(post_vals[mask] == pre_vals[mask])
                    print(f"  {src:<15} pre={np.mean(pre_vals[mask]):.4f} -> post={np.mean(post_vals[mask]):.4f} "
                          f"(Δ={delta:+.4f}, p={pval:.2e}, improved={improved}, degraded={degraded}, tied={tied})")

    # === Analysis 4: Replicate convergence ===
    print("\n" + "="*80)
    print("REPLICATE CONVERGENCE (std dev across 5 reps)")
    print("="*80)

    rep_std = per_rep.groupby(['pipeline', 'source', 'target', 'protocol'])[
        'tmscore'].std().reset_index()
    rep_std.columns = ['pipeline', 'source', 'target', 'protocol', 'rep_std']

    for pipe in ['blue', 'green']:
        pipe_std = rep_std[rep_std['pipeline'] == pipe]
        print(f"\n{pipe.upper()}: Mean std dev across reps = {pipe_std['rep_std'].mean():.4f}")
        print(f"  Median = {pipe_std['rep_std'].median():.4f}")
        print(f"  Max = {pipe_std['rep_std'].max():.4f}")
        print(f"  % targets with std < 0.01: {100*(pipe_std['rep_std'] < 0.01).mean():.1f}%")
        print(f"  % targets with std < 0.005: {100*(pipe_std['rep_std'] < 0.005).mean():.1f}%")

    # === Analysis 5: Blue-Green agreement ===
    print("\n" + "="*80)
    print("BLUE-GREEN AGREEMENT ON ROSETTA TM-SCORE")
    print("="*80)

    for proto in protocols:
        blue_proto = per_target_proto[(per_target_proto['pipeline'] == 'blue') &
                                      (per_target_proto['protocol'] == proto)]
        green_proto = per_target_proto[(per_target_proto['pipeline'] == 'green') &
                                       (per_target_proto['protocol'] == proto)]

        b_by_target = blue_proto.groupby('target')['tmscore'].mean()
        g_by_target = green_proto.groupby('target')['tmscore'].mean()
        common = b_by_target.index.intersection(g_by_target.index)

        if len(common) > 10:
            bv = b_by_target.loc[common].values
            gv = g_by_target.loc[common].values
            mask = np.isfinite(bv) & np.isfinite(gv)
            if mask.sum() > 10:
                r, p = stats.pearsonr(bv[mask], gv[mask])
                print(f"  {proto:<20} r={r:.4f}, p={p:.2e}, n={mask.sum()}")

    # Save summary tables
    proto_summary = per_target_proto.groupby(['pipeline', 'protocol'])['tmscore'].agg(
        ['mean', 'median', 'std', 'count']).reset_index()
    proto_summary.to_csv(os.path.join(OUTDIR, "rosetta_protocol_summary.tsv"), sep='\t', index=False)

    src_proto_summary = per_target_proto.groupby(['pipeline', 'source', 'protocol'])['tmscore'].agg(
        ['mean', 'std', 'count']).reset_index()
    src_proto_summary.to_csv(os.path.join(OUTDIR, "rosetta_source_protocol_summary.tsv"), sep='\t', index=False)

    print(f"\n\nSummary tables saved to {OUTDIR}")
    print("Done!")


if __name__ == '__main__':
    main()
