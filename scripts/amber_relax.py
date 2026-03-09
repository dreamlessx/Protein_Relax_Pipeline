#!/usr/bin/env python3
"""
Standalone AMBER relaxation using AlphaFold 2.3.2's relax module.
Relaxes AF unrelaxed + Boltz predictions with GPU-accelerated OpenMM.

Usage: python amber_relax.py <target_dir>
"""

import sys
import os
import glob
import time
import traceback

# Add AF2 to path
sys.path.insert(0, "/sb/apps/alphafold232/alphafold")

from alphafold.common import protein
from alphafold.relax import relax


def relax_pdb(input_path, output_path, use_gpu=True):
    """Relax a single PDB file using AF2 AMBER relaxation."""
    with open(input_path) as f:
        pdb_string = f.read()

    prot = protein.from_pdb_string(pdb_string)

    amber_relaxer = relax.AmberRelaxation(
        max_iterations=0,
        tolerance=2.39,
        stiffness=10.0,
        exclude_residues=[],
        max_outer_iterations=3,
        use_gpu=use_gpu,
    )

    relaxed_pdb, debug_data, violations = amber_relaxer.process(prot=prot)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        f.write(relaxed_pdb)

    return {
        "initial_energy": debug_data.get("initial_energy", None),
        "final_energy": debug_data.get("final_energy", None),
        "rmsd": debug_data.get("rmsd", None),
        "attempts": debug_data.get("attempts", None),
        "num_violations": sum(violations) if violations is not None else 0,
    }


def main():
    if len(sys.argv) < 2:
        print("Usage: python amber_relax.py <target_dir>")
        sys.exit(1)

    target_dir = sys.argv[1]
    target = os.path.basename(target_dir)
    amber_dir = os.path.join(target_dir, "amber_out")

    print(f"=== AMBER relaxation for {target} ===")

    # Collect all PDBs to relax
    jobs = []

    # AF unrelaxed
    for i in range(5):
        pdb = os.path.join(target_dir, "af_out_unrelaxed", f"ranked_{i}.pdb")
        if os.path.exists(pdb):
            out = os.path.join(amber_dir, f"af_unrelaxed_ranked_{i}", "relaxed.pdb")
            jobs.append(("af_unrelaxed", f"ranked_{i}", pdb, out))

    # Boltz
    boltz_dir = os.path.join(
        target_dir,
        "boltz_out",
        "boltz_results_boltz_input",
        "predictions",
        "boltz_input",
    )
    for i in range(5):
        pdb = os.path.join(boltz_dir, f"boltz_input_model_{i}.pdb")
        if os.path.exists(pdb):
            out = os.path.join(amber_dir, f"boltz_model_{i}", "relaxed.pdb")
            jobs.append(("boltz", f"model_{i}", pdb, out))

    if not jobs:
        print(f"No PDBs found for {target}")
        sys.exit(0)

    print(f"Found {len(jobs)} PDBs to relax")

    succeeded = 0
    failed = 0

    for src_type, model_name, input_path, output_path in jobs:
        # Skip if already done
        if os.path.exists(output_path):
            print(f"  SKIP {src_type}/{model_name} (already relaxed)")
            succeeded += 1
            continue

        print(f"  Relaxing {src_type}/{model_name}...", end=" ", flush=True)
        t0 = time.time()

        try:
            result = relax_pdb(input_path, output_path, use_gpu=True)
            elapsed = time.time() - t0
            print(
                f"OK ({elapsed:.0f}s) "
                f"E: {result['initial_energy']:.0f} -> {result['final_energy']:.0f} "
                f"RMSD: {result['rmsd']:.3f} "
                f"violations: {result['num_violations']}"
            )
            succeeded += 1
        except Exception as e:
            elapsed = time.time() - t0
            print(f"FAILED ({elapsed:.0f}s): {e}")
            traceback.print_exc()
            failed += 1

    print(f"\n=== {target}: {succeeded} succeeded, {failed} failed ===")

    if failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
