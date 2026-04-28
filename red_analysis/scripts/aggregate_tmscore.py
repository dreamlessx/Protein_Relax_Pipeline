#!/usr/bin/env python3
"""
aggregate_tmscore.py - Red Analysis Pipeline

Merges per-target TM-score TSV files into a single combined dataset,
computes summary statistics, and generates initial comparison tables.

Output files:
  - metrics/combined_tmscore.tsv    - all rows merged
  - tables/summary_by_source.tsv   - mean/median/std per source type
  - tables/summary_by_pipeline.tsv - Blue vs Green comparison
  - tables/blue_green_agreement.tsv - per-target correlation metrics

QUESTION: Should we report TM-score or 1-TM-score? TM-score is bounded [0,1]
and most of our values will be very high (>0.95). Reporting 1-TM might make
differences more visible. For now, reporting both.

QUESTION: For multi-model sources (5 models each), do we use best-of-5,
mean-of-5, or report all? Paper convention varies. Using mean-of-5 for
summary stats but keeping all rows in combined file.
"""

import os
import glob
import sys
import numpy as np
from collections import defaultdict

METRICS_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"
TABLES_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/tables"

def load_all_tsvs():
    """Load all per-target TSV files into a list of dicts."""
    rows = []
    files = sorted(glob.glob(os.path.join(METRICS_DIR, "tmscore_*.tsv")))

    if not files:
        print("ERROR: No TSV files found in", METRICS_DIR)
        sys.exit(1)

    print(f"Loading {len(files)} TSV files...")

    for f in files:
        with open(f) as fh:
            header = fh.readline().strip().split('\t')
            for line in fh:
                vals = line.strip().split('\t')
                if len(vals) != len(header):
                    print(f"WARN: Skipping malformed line in {f}: {line.strip()}")
                    continue
                row = dict(zip(header, vals))
                rows.append(row)

    return rows, header


def write_combined(rows, header):
    """Write combined TSV with all rows."""
    outf = os.path.join(METRICS_DIR, "combined_tmscore.tsv")
    with open(outf, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for row in rows:
            f.write('\t'.join(row.get(h, 'NA') for h in header) + '\n')
    print(f"Combined: {len(rows)} rows -> {outf}")
    return outf


def compute_summary_by_source(rows):
    """
    Compute mean/median/std of RMSD, TM-score, GDT-TS, GDT-HA
    grouped by (pipeline, source).

    Uses mean-of-5-models per target, then aggregates across targets.
    This is the standard approach: average over stochastic replicates
    first, then compare across targets.
    """
    # Group by (pipeline, source, target) -> list of metric values
    groups = defaultdict(lambda: defaultdict(list))

    metrics_to_agg = ['rmsd', 'tmscore', 'gdtts', 'gdtha']

    for row in rows:
        key = (row['pipeline'], row['source'], row['target'])
        for m in metrics_to_agg:
            val = row.get(m, 'NA')
            if val != 'NA':
                try:
                    groups[key][m].append(float(val))
                except ValueError:
                    pass

    # Average across models per target
    # key: (pipeline, source) -> list of per-target means
    agg = defaultdict(lambda: defaultdict(list))
    for (pipeline, source, target), metrics in groups.items():
        for m in metrics_to_agg:
            if m in metrics and len(metrics[m]) > 0:
                agg[(pipeline, source)][m].append(np.mean(metrics[m]))

    # Write summary table
    outf = os.path.join(TABLES_DIR, "summary_by_source.tsv")
    with open(outf, 'w') as f:
        header = ['pipeline', 'source', 'n_targets']
        for m in metrics_to_agg:
            header.extend([f'{m}_mean', f'{m}_median', f'{m}_std', f'{m}_min', f'{m}_max'])
        f.write('\t'.join(header) + '\n')

        for (pipeline, source) in sorted(agg.keys()):
            vals = agg[(pipeline, source)]
            n = len(vals.get('rmsd', []))
            row = [pipeline, source, str(n)]
            for m in metrics_to_agg:
                v = np.array(vals.get(m, []))
                if len(v) > 0:
                    row.extend([
                        f'{np.mean(v):.4f}',
                        f'{np.median(v):.4f}',
                        f'{np.std(v):.4f}',
                        f'{np.min(v):.4f}',
                        f'{np.max(v):.4f}'
                    ])
                else:
                    row.extend(['NA'] * 5)
            f.write('\t'.join(row) + '\n')

    print(f"Summary by source -> {outf}")
    return outf


def compute_blue_green_agreement(rows):
    """
    For each (target, source, model_idx), compare Blue vs Green metrics.

    This answers: do independent runs produce the same structural quality?
    We expect high correlation but NOT identical values (AF/Boltz are stochastic).

    COMMENT: This is the core reproducibility analysis. If Blue and Green
    disagree substantially, it indicates sensitivity to stochastic seeds,
    which would be an important finding in itself.
    """
    # Index by (target, source, model_idx, pipeline)
    index = {}
    for row in rows:
        key = (row['target'], row['source'], row['model_idx'], row['pipeline'])
        index[key] = row

    # Find matched pairs
    outf = os.path.join(TABLES_DIR, "blue_green_agreement.tsv")
    with open(outf, 'w') as f:
        f.write('target\tsource\tmodel_idx\t')
        f.write('blue_rmsd\tgreen_rmsd\trmsd_diff\t')
        f.write('blue_tmscore\tgreen_tmscore\ttmscore_diff\t')
        f.write('blue_gdtts\tgreen_gdtts\tgdtts_diff\n')

        n_matched = 0
        n_missing = 0
        diffs_rmsd = []
        diffs_tm = []

        targets = sorted(set(r['target'] for r in rows))
        sources = sorted(set(r['source'] for r in rows))

        for target in targets:
            for source in sources:
                for idx in ['0', '1', '2', '3', '4']:
                    blue_key = (target, source, idx, 'blue')
                    green_key = (target, source, idx, 'green')

                    if blue_key in index and green_key in index:
                        b = index[blue_key]
                        g = index[green_key]

                        try:
                            b_rmsd = float(b.get('rmsd', 'NA'))
                            g_rmsd = float(g.get('rmsd', 'NA'))
                            b_tm = float(b.get('tmscore', 'NA'))
                            g_tm = float(g.get('tmscore', 'NA'))
                            b_gdt = float(b.get('gdtts', 'NA'))
                            g_gdt = float(g.get('gdtts', 'NA'))

                            d_rmsd = b_rmsd - g_rmsd
                            d_tm = b_tm - g_tm
                            d_gdt = b_gdt - g_gdt

                            diffs_rmsd.append(abs(d_rmsd))
                            diffs_tm.append(abs(d_tm))

                            f.write(f'{target}\t{source}\t{idx}\t')
                            f.write(f'{b_rmsd:.4f}\t{g_rmsd:.4f}\t{d_rmsd:.4f}\t')
                            f.write(f'{b_tm:.4f}\t{g_tm:.4f}\t{d_tm:.4f}\t')
                            f.write(f'{b_gdt:.4f}\t{g_gdt:.4f}\t{d_gdt:.4f}\n')
                            n_matched += 1
                        except (ValueError, TypeError):
                            n_missing += 1
                    else:
                        n_missing += 1

    print(f"\nBlue-Green agreement: {n_matched} matched pairs, {n_missing} missing")
    if diffs_rmsd:
        print(f"  RMSD diff:     mean={np.mean(diffs_rmsd):.4f}, "
              f"median={np.median(diffs_rmsd):.4f}, max={np.max(diffs_rmsd):.4f}")
        print(f"  TM-score diff: mean={np.mean(diffs_tm):.6f}, "
              f"median={np.median(diffs_tm):.6f}, max={np.max(diffs_tm):.6f}")

    return outf


def print_quick_summary(rows):
    """Print a quick human-readable summary to stdout."""
    print("\n" + "="*70)
    print("RED ANALYSIS - Pre-Rosetta TM-score Summary")
    print("="*70)

    # Group by (pipeline, source)
    groups = defaultdict(list)
    for row in rows:
        try:
            tm = float(row.get('tmscore', 'NA'))
            rmsd = float(row.get('rmsd', 'NA'))
            groups[(row['pipeline'], row['source'])].append((tm, rmsd))
        except ValueError:
            pass

    print(f"\n{'Pipeline':<10} {'Source':<15} {'N':>5} {'TM mean':>10} {'TM med':>10} {'RMSD mean':>10} {'RMSD med':>10}")
    print("-" * 70)

    for (pipeline, source) in sorted(groups.keys()):
        vals = groups[(pipeline, source)]
        tms = [v[0] for v in vals]
        rmsds = [v[1] for v in vals]
        print(f"{pipeline:<10} {source:<15} {len(vals):>5} "
              f"{np.mean(tms):>10.4f} {np.median(tms):>10.4f} "
              f"{np.mean(rmsds):>10.4f} {np.median(rmsds):>10.4f}")

    # Key question: does AMBER relaxation help or hurt?
    print("\n" + "="*70)
    print("KEY COMPARISON: Does relaxation change structural similarity to crystal?")
    print("="*70)

    for pipeline in ['blue', 'green']:
        print(f"\n--- {pipeline.upper()} ---")

        # AF unrelaxed vs AF relaxed (built-in AMBER)
        af_unrelaxed = [v[0] for v in groups.get((pipeline, 'af_unrelaxed'), [])]
        af_relaxed = [v[0] for v in groups.get((pipeline, 'af_relaxed'), [])]
        amber_af = [v[0] for v in groups.get((pipeline, 'amber_af'), [])]
        boltz = [v[0] for v in groups.get((pipeline, 'boltz'), [])]
        amber_boltz = [v[0] for v in groups.get((pipeline, 'amber_boltz'), [])]

        if af_unrelaxed and af_relaxed:
            print(f"  AF unrelaxed TM: {np.mean(af_unrelaxed):.4f} ± {np.std(af_unrelaxed):.4f}")
            print(f"  AF relaxed (built-in AMBER) TM: {np.mean(af_relaxed):.4f} ± {np.std(af_relaxed):.4f}")
            print(f"  -> AMBER effect on AF: {np.mean(af_relaxed) - np.mean(af_unrelaxed):+.4f}")

        if af_unrelaxed and amber_af:
            print(f"  Standalone AMBER(AF) TM: {np.mean(amber_af):.4f} ± {np.std(amber_af):.4f}")
            print(f"  -> Standalone AMBER effect on AF: {np.mean(amber_af) - np.mean(af_unrelaxed):+.4f}")

        if boltz and amber_boltz:
            print(f"  Boltz TM: {np.mean(boltz):.4f} ± {np.std(boltz):.4f}")
            print(f"  AMBER(Boltz) TM: {np.mean(amber_boltz):.4f} ± {np.std(amber_boltz):.4f}")
            print(f"  -> AMBER effect on Boltz: {np.mean(amber_boltz) - np.mean(boltz):+.4f}")


if __name__ == '__main__':
    rows, header = load_all_tsvs()
    write_combined(rows, header)
    compute_summary_by_source(rows)
    compute_blue_green_agreement(rows)
    print_quick_summary(rows)
