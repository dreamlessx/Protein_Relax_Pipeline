#!/usr/bin/env python3
"""
extract_rosetta_energy.py — Red Analysis Pipeline

Extracts total_score (REU) from Rosetta relax.fasc files for all targets.
Also computes per-residue energy (total_score / seq_len) for cross-target comparison.

Output: metrics/combined_rosetta_energy.tsv
"""

import os
import sys
import glob
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

BLUE_BASE = "/data/p_csb_meiler/agarwm5/af_work"
GREEN_BASE = "/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/data"
OUTDIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"

PROTOCOLS = ['cartesian_beta', 'cartesian_ref15', 'dualspace_beta',
             'dualspace_ref15', 'normal_beta', 'normal_ref15']

# Source mapping: dir name prefix -> simplified source name
def classify_source(src_name):
    """Map src_type directory name to simplified source category."""
    if src_name.startswith('af_relaxed'):
        return 'af_relaxed'
    elif src_name.startswith('af_unrelaxed'):
        return 'af_unrelaxed'
    elif src_name.startswith('amber_af'):
        return 'amber_af'
    elif src_name.startswith('amber_boltz'):
        return 'amber_boltz'
    elif src_name.startswith('boltz'):
        return 'boltz'
    elif src_name.startswith('crystal'):
        return 'crystal'
    return 'unknown'


def get_seq_len(target, base_dir):
    """Get sequence length from FASTA file for per-residue energy."""
    fasta = os.path.join(base_dir, target, 'sequence.fasta')
    if os.path.exists(fasta):
        with open(fasta) as f:
            lines = [l.strip() for l in f if not l.startswith('>')]
            return len(''.join(lines))
    return None


def extract_target(args):
    """Extract energy for one target, one pipeline."""
    target, pipeline = args
    base_dir = BLUE_BASE if pipeline == 'blue' else GREEN_BASE
    rosetta_dir = os.path.join(base_dir, target, 'rosetta_out')

    if not os.path.isdir(rosetta_dir):
        return []

    seq_len = get_seq_len(target, base_dir)
    rows = []

    for src_dir in sorted(glob.glob(os.path.join(rosetta_dir, '*/'))):
        src_name = os.path.basename(src_dir.rstrip('/'))
        source = classify_source(src_name)

        for protocol in PROTOCOLS:
            fasc_path = os.path.join(src_dir, protocol, 'relax.fasc')
            if not os.path.exists(fasc_path):
                continue

            try:
                with open(fasc_path) as f:
                    for line in f:
                        if line.startswith('SCORE:') and 'total_score' not in line:
                            parts = line.split()
                            if len(parts) >= 3:
                                total_score = float(parts[1])
                                desc = parts[-1]  # e.g., ranked_0_r1
                                # Extract rep number from description
                                rep = desc.split('_r')[-1] if '_r' in desc else '1'
                                per_res = total_score / seq_len if seq_len else ''
                                rows.append(f"{target}\t{pipeline}\t{src_name}\t"
                                           f"{protocol}\t{rep}\t{total_score}\t"
                                           f"{per_res}\t{source}")
            except Exception as e:
                print(f"  WARN: {fasc_path}: {e}", file=sys.stderr)

    return rows


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    outfile = os.path.join(OUTDIR, "combined_rosetta_energy.tsv")

    # Get target list from existing MolProbity data
    mp_path = os.path.join(OUTDIR, "combined_rosetta_molprobity.tsv")
    import pandas as pd
    mp = pd.read_csv(mp_path, sep='\t')
    targets = sorted(mp['target'].unique())

    print(f"Extracting Rosetta energy for {len(targets)} targets...")

    # Build task list: (target, pipeline) pairs
    tasks = []
    for target in targets:
        for pipeline in ['blue', 'green']:
            tasks.append((target, pipeline))

    all_rows = []
    with ProcessPoolExecutor(max_workers=16) as executor:
        futures = {executor.submit(extract_target, t): t for t in tasks}
        done = 0
        for future in as_completed(futures):
            done += 1
            rows = future.result()
            all_rows.extend(rows)
            if done % 100 == 0:
                print(f"  {done}/{len(tasks)} tasks, {len(all_rows)} rows...")

    # Write output
    with open(outfile, 'w') as f:
        f.write("target\tpipeline\tsrc_type\tprotocol\trep\ttotal_score\t"
                "per_residue_energy\tsource\n")
        for row in sorted(all_rows):
            f.write(row + '\n')

    print(f"\nDone: {len(all_rows)} rows written to {outfile}")
    unique_targets = set(r.split('\t')[0] for r in all_rows)
    print(f"Targets: {len(unique_targets)}")


if __name__ == '__main__':
    main()
