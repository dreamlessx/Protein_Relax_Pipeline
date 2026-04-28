#!/usr/bin/env python3
"""
score_initial_energy.py - Red Analysis Pipeline

Score input PDBs (pre-Rosetta) with Rosetta's score function to get initial energies.
Uses score_jd2 binary to compute total_score for each input structure.

Run as SLURM array: one task per target.
Output: metrics/initial_energy_{target}.tsv
"""

import os
import sys
import glob
import subprocess
import tempfile
import gzip
import shutil

ROSETTA_SCORE = "/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/source/bin/score_jd2.linuxgccrelease"
ROSETTA_DB = "/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/database"
GREEN_BASE = "/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/data"
CRYSTAL_DIR = "/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/cleaned"
TARGET_LIST = "/data/p_csb_meiler/agarwm5/red_analysis/target_list.txt"
OUTDIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics/initial_energy"

# Source directories and their PDB patterns
SOURCES = {
    'af_relaxed':    ('AF/ranked_{i}.pdb', 5),
    'af_unrelaxed':  ('af_out_unrelaxed/ranked_{i}.pdb', 5),
    'amber_af':      ('amber_out/af_unrelaxed_{i}/relaxed.pdb', 5),
    'amber_boltz':   ('amber_out/boltz_model_{i}/relaxed.pdb', 5),
    'boltz':         ('boltz_out/boltz_input_model_{i}.pdb', 5),
    'crystal':       (None, 1),  # special case
}


def get_seq_len(target):
    """Get sequence length for per-residue energy."""
    fasta = os.path.join(GREEN_BASE, target, 'sequence.fasta')
    if os.path.exists(fasta):
        with open(fasta) as f:
            lines = [l.strip() for l in f if not l.startswith('>')]
            return len(''.join(lines))
    return None


def score_pdb(pdb_path, workdir):
    """Score a single PDB with Rosetta and return total_score."""
    if not os.path.exists(pdb_path):
        return None

    # If gzipped, decompress first
    actual_path = pdb_path
    tmpfile = None
    if pdb_path.endswith('.gz'):
        tmpfile = os.path.join(workdir, 'input.pdb')
        with gzip.open(pdb_path, 'rb') as f_in:
            with open(tmpfile, 'wb') as f_out:
                f_out.write(f_in.read())
        actual_path = tmpfile

    scorefile = os.path.join(workdir, 'score.sc')
    if os.path.exists(scorefile):
        os.remove(scorefile)

    cmd = [
        ROSETTA_SCORE,
        '-database', ROSETTA_DB,
        '-in:file:s', actual_path,
        '-out:file:scorefile', scorefile,
        '-ignore_unrecognized_res',
        '-ignore_zero_occupancy', 'false',
    ]

    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if os.path.exists(scorefile):
            with open(scorefile) as f:
                for line in f:
                    if line.startswith('SCORE:') and 'total_score' not in line:
                        parts = line.split()
                        return float(parts[1])
    except Exception as e:
        print(f"  WARN: {pdb_path}: {e}", file=sys.stderr)

    return None


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    # Get target from SLURM or CLI
    slurm_id = os.environ.get('SLURM_ARRAY_TASK_ID')
    if slurm_id:
        with open(TARGET_LIST) as f:
            targets = f.read().strip().split('\n')
        target = targets[int(slurm_id) - 1]
    else:
        target = sys.argv[1] if len(sys.argv) > 1 else '1A2K'

    outfile = os.path.join(OUTDIR, f"initial_energy_{target}.tsv")

    # Skip if ALL 6 sources already scored (expect: 5+5+5+5+5+1 = 26 rows)
    if os.path.exists(outfile):
        with open(outfile) as f:
            n = sum(1 for _ in f) - 1
        if n >= 20:  # at least 4 sources × 5 models
            print(f"SKIP: {target} - {n} rows exist (complete)")
            return

    # Load existing sources to avoid re-scoring
    existing_sources = set()
    if os.path.exists(outfile):
        import pandas as pd
        existing = pd.read_csv(outfile, sep='\t')
        existing_sources = set(existing['source'].unique())
        print(f"  Existing sources: {existing_sources}")

    seq_len = get_seq_len(target)
    print(f"=== Scoring {target} (seq_len={seq_len}) ===")

    rows = []
    with tempfile.TemporaryDirectory() as workdir:
        for source, (pattern, n_models) in SOURCES.items():
            if source in existing_sources:
                print(f"  SKIP {source}: already scored")
                continue
            for i in range(n_models):
                if source == 'crystal':
                    pdb_path = os.path.join(CRYSTAL_DIR, f"{target}.pdb")
                else:
                    pdb_path = os.path.join(GREEN_BASE, target, pattern.format(i=i))

                if not os.path.exists(pdb_path):
                    print(f"  WARN {source} model_{i}: {pdb_path} not found")
                    continue

                total_score = score_pdb(pdb_path, workdir)
                if total_score is not None:
                    per_res = total_score / seq_len if seq_len else ''
                    rows.append(f"{target}\tgreen\t{source}\t{i}\t{total_score}\t{per_res}")
                    print(f"  {source} model_{i}: {total_score:.1f} REU")

    # Append new rows to existing file, or write fresh
    if existing_sources and os.path.exists(outfile):
        with open(outfile, 'a') as f:
            for row in rows:
                f.write(row + '\n')
    else:
        with open(outfile, 'w') as f:
            f.write("target\tpipeline\tsource\tmodel_idx\ttotal_score\tper_residue_energy\n")
            for row in rows:
                f.write(row + '\n')

    print(f"Done: {len(rows)} new rows -> {outfile}")


if __name__ == '__main__':
    main()
