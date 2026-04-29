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
import gzip
import re
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
    elif src_name.startswith('amber_crystal'):
        return 'amber_crystal'
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


import re

# Patterns for extracting rep number from Rosetta SCORE description column.
# `relax.fasc` uses `<stem>_r<N>` (canonical). When relax.fasc is incomplete,
# fill-in sidecars exist:
#   score_tmp{N}.sc       -> description like "<target>_tmp<N>"
#   score_TMP{N}.sc       -> description like "relaxed_TMP<N>"
#   score_fill{N}.sc      -> description like "<stem>_fill<N>"
# Each holds a single SCORE line for one rep. Combine all source files per
# (src_type, protocol) cell, dedup by rep, prefer relax.fasc when overlap.
REP_RE = re.compile(r'_(?:r|tmp|TMP|fill)(\d+)$')


def _parse_score_file(path, dest, source_priority):
    """Parse one .sc / .fasc file. Append (rep, total_score) into `dest` dict
    only if rep not already present from a higher-priority source."""
    try:
        with open(path) as f:
            for line in f:
                if not line.startswith('SCORE:'):
                    continue
                if 'total_score' in line:
                    continue  # header line
                parts = line.split()
                if len(parts) < 3:
                    continue
                try:
                    total_score = float(parts[1])
                except ValueError:
                    continue
                desc = parts[-1]
                m = REP_RE.search(desc)
                if not m:
                    continue
                rep = int(m.group(1))
                if rep < 1 or rep > 5:
                    continue
                if rep in dest:
                    # Already have this rep from higher- or equal-priority source
                    if dest[rep][1] <= source_priority:
                        continue
                dest[rep] = (total_score, source_priority)
    except (OSError, IOError) as e:
        print(f"  WARN: {path}: {e}", file=sys.stderr)


PDB_REP_RE = re.compile(r'_r(\d+)\.pdb(?:\.gz)?$')


def _parse_pdb_total_score(path):
    """Parse Rosetta-output PDB; return total_score from POSE_ENERGIES_TABLE
    `pose` line (last column). None on failure."""
    try:
        opener = gzip.open if path.endswith('.gz') else open
        in_table = False
        with opener(path, 'rt') as f:
            for line in f:
                if line.startswith('#BEGIN_POSE_ENERGIES_TABLE'):
                    in_table = True
                    continue
                if line.startswith('#END_POSE_ENERGIES_TABLE'):
                    break
                if in_table and line.startswith('pose '):
                    parts = line.split()
                    return float(parts[-1])
    except (OSError, IOError, ValueError):
        return None
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
            proto_dir = os.path.join(src_dir, protocol)
            if not os.path.isdir(proto_dir):
                continue

            # Collect rep -> (total_score, priority) across sources, lowest priority wins:
            #   priority 0: relax.fasc (canonical)
            #   priority 1: score_{tmp,TMP,fill}*.sc sidecars
            #   priority 2: PDB POSE_ENERGIES_TABLE pose line (only if no fasc/sc covers it)
            scores = {}

            fasc = os.path.join(proto_dir, 'relax.fasc')
            if os.path.exists(fasc):
                _parse_score_file(fasc, scores, source_priority=0)

            for sidecar in sorted(glob.glob(os.path.join(proto_dir, 'score_*.sc'))):
                _parse_score_file(sidecar, scores, source_priority=1)

            # Final fallback: parse PDB POSE_ENERGIES_TABLE for any rep still missing.
            for pdb_path in sorted(glob.glob(os.path.join(proto_dir, '*_r*.pdb.gz'))):
                m = PDB_REP_RE.search(pdb_path)
                if not m:
                    continue
                rep = int(m.group(1))
                if rep < 1 or rep > 5 or rep in scores:
                    continue
                ts = _parse_pdb_total_score(pdb_path)
                if ts is not None:
                    scores[rep] = (ts, 2)

            for rep in sorted(scores.keys()):
                total_score = scores[rep][0]
                per_res = total_score / seq_len if seq_len else ''
                rows.append(f"{target}\t{pipeline}\t{src_name}\t"
                           f"{protocol}\t{rep}\t{total_score}\t"
                           f"{per_res}\t{source}")

    return rows


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    outfile = os.path.join(OUTDIR, "combined_rosetta_energy.tsv")

    # Get target list (use target_list.txt as canonical; fallback to MP TSV)
    target_list_path = "/data/p_csb_meiler/agarwm5/red_analysis/target_list.txt"
    if os.path.exists(target_list_path):
        with open(target_list_path) as f:
            targets = sorted(t.strip() for t in f if t.strip())
    else:
        mp_path = os.path.join(OUTDIR, "combined_rosetta_molprobity.tsv")
        seen = set()
        with open(mp_path) as f:
            next(f)
            for line in f:
                seen.add(line.split('\t', 1)[0])
        targets = sorted(seen)

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
