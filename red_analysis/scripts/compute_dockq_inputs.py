#!/usr/bin/env python3
"""
compute_dockq_inputs.py - BM5.5 interface-quality scoring on input structures.

For each (target, pipeline, src_type) cell, runs DockQ(prediction, crystal_reference)
and emits DockQ, i-RMSD, l-RMSD, f_nat. Scope is the 27 INPUT structures per
(target, pipeline) cell, NOT the post-Rosetta outputs.

Pipelines:
  blue  -> /data/p_csb_meiler/agarwm5/af_work
  green -> /data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/data

Source-type catalogue per (target, pipeline) -> 27 PDBs:
  af_relaxed_ranked_<0..4>     5  af_out/ranked_<N>.pdb
  af_unrelaxed_ranked_<0..4>   5  af_out_unrelaxed/ranked_<N>.pdb
  boltz_input_model_<0..4>     5  boltz_out/.../boltz_input_model_<N>.pdb
                                  (green: Boltz/boltz_input_model_<N>.pdb)
  amber_af_<0..4>              5  amber_out/af_unrelaxed_<N>/relaxed.pdb
                                  (blue uses af_unrelaxed_ranked_<N>)
  amber_boltz_<0..4>           5  amber_out/boltz_model_<N>/relaxed.pdb
  amber_crystal                1  amber_out/crystal/relaxed.pdb (green only)
  crystal                      1  cleaned/<TARGET>.pdb (DockQ vs self -> sanity 1.0)

Crystal reference for all DockQ calls:
  /data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/cleaned/<TARGET>.pdb

Output TSV columns:
  target  pipeline  source  src_type  dockq  i_rmsd  l_rmsd  f_nat  n_interfaces  status

Status:
  ok           DockQ computed
  miss_file    prediction PDB not found
  miss_native  crystal reference not found
  err:<msg>    DockQ raised
"""

import argparse
import gzip
import os
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# DockQ python API
from DockQ.DockQ import (
    load_PDB,
    run_on_all_native_interfaces,
    get_all_chain_maps,
    group_chains,
    count_chain_combinations,
)

BLUE_BASE = "/data/p_csb_meiler/agarwm5/af_work"
GREEN_BASE = "/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/data"
CRYSTAL_BASE = "/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/cleaned"
TARGET_LIST = "/data/p_csb_meiler/agarwm5/red_analysis/target_list.txt"
OUT_TSV = "/data/p_csb_meiler/agarwm5/red_analysis/metrics/combined_dockq_input.tsv"


def expected_inputs(target, pipeline):
    """Return list of (src_type, abs_path) tuples for the 27 input PDBs."""
    base = BLUE_BASE if pipeline == "blue" else GREEN_BASE
    tdir = os.path.join(base, target)
    out = []

    # af_relaxed -> ranked_<N>.pdb
    # blue:  af_out/ranked_<N>.pdb
    # green: AF/ranked_<N>.pdb (af_out/ has only sequence/, not pdbs)
    for n in range(5):
        if pipeline == "blue":
            p = os.path.join(tdir, "af_out", f"ranked_{n}.pdb")
        else:
            p = os.path.join(tdir, "AF", f"ranked_{n}.pdb")
        out.append((f"af_relaxed_ranked_{n}", p))

    # af_unrelaxed -> ranked_<N>.pdb (in af_out_unrelaxed/)
    for n in range(5):
        out.append((f"af_unrelaxed_ranked_{n}",
                    os.path.join(tdir, "af_out_unrelaxed", f"ranked_{n}.pdb")))

    # boltz -> boltz_input_model_<N>.pdb
    # blue:  boltz_out/boltz_results_boltz_input/predictions/boltz_input/boltz_input_model_<N>.pdb
    # green: Boltz/boltz_input_model_<N>.pdb
    for n in range(5):
        if pipeline == "blue":
            p = os.path.join(tdir, "boltz_out", "boltz_results_boltz_input",
                             "predictions", "boltz_input", f"boltz_input_model_{n}.pdb")
        else:
            p = os.path.join(tdir, "Boltz", f"boltz_input_model_{n}.pdb")
        out.append((f"boltz_boltz_input_model_{n}", p))

    # amber_af -> amber_out/af_unrelaxed_<N>/relaxed.pdb (blue uses af_unrelaxed_ranked_<N>)
    for n in range(5):
        if pipeline == "blue":
            sub = f"af_unrelaxed_ranked_{n}"
        else:
            sub = f"af_unrelaxed_{n}"
        out.append((f"amber_af_unrelaxed_ranked_{n}",
                    os.path.join(tdir, "amber_out", sub, "relaxed.pdb")))

    # amber_boltz -> amber_out/boltz_model_<N>/relaxed.pdb
    for n in range(5):
        out.append((f"amber_boltz_boltz_model_{n}",
                    os.path.join(tdir, "amber_out", f"boltz_model_{n}", "relaxed.pdb")))

    # amber_crystal -> amber_out/crystal/relaxed.pdb (both pipelines)
    out.append(("amber_crystal", os.path.join(tdir, "amber_out", "crystal", "relaxed.pdb")))

    # crystal -> cleaned/<TARGET>.pdb (sanity DockQ=1.0 self-comparison)
    out.append((f"crystal_{target}", os.path.join(CRYSTAL_BASE, f"{target}.pdb")))

    return out


def classify_source(src_type):
    """Map src_type -> simplified source bucket (matches BM5.5 7-bucket scheme)."""
    if src_type.startswith("af_relaxed"):
        return "af_relaxed"
    if src_type.startswith("af_unrelaxed"):
        return "af_unrelaxed"
    if src_type.startswith("amber_af"):
        return "amber_af"
    if src_type.startswith("amber_boltz"):
        return "amber_boltz"
    if src_type == "amber_crystal":
        return "amber_crystal"
    if src_type.startswith("boltz"):
        return "boltz"
    if src_type.startswith("crystal"):
        return "crystal"
    return "unknown"


def best_dockq_over_chain_perms(model_path, native_path):
    """Mirror DockQ CLI: load both, enumerate chain maps, pick best total DockQ.
    Returns: (mean_dockq, mean_irmsd, mean_lrmsd, mean_fnat, n_interfaces).

    For multi-interface natives we report the mean of per-interface metrics
    (the "best" mapping is selected by total DockQ summed across interfaces;
    the reported scalars are then averaged for ingestion convenience).
    """
    model_structure = load_PDB(model_path)
    native_structure = load_PDB(native_path)

    model_chains = [c.id for c in model_structure]
    native_chains = [c.id for c in native_structure]

    if len(model_chains) < 2 or len(native_chains) < 2:
        raise RuntimeError(
            f"need >=2 chains, got model={model_chains} native={native_chains}"
        )

    chain_clusters, reverse_map = group_chains(
        model_structure, native_structure,
        model_chains, native_chains,
        allowed_mismatches=0,
    )
    chain_maps = list(get_all_chain_maps(
        chain_clusters, {}, reverse_map, model_chains, native_chains,
    ))
    if not chain_maps:
        raise RuntimeError("no candidate chain maps")

    best_dockq = -1.0
    best_result = None
    for cm in chain_maps:
        result, total = run_on_all_native_interfaces(
            model_structure, native_structure,
            chain_map=cm, no_align=False,
        )
        if result and total > best_dockq:
            best_dockq = total
            best_result = result

    if not best_result:
        raise RuntimeError("no native interfaces found")

    n = len(best_result)
    dockqs = [info["DockQ"] for info in best_result.values()]
    irmsds = [info.get("iRMSD", float("nan")) for info in best_result.values()]
    lrmsds = [info.get("LRMSD", float("nan")) for info in best_result.values()]
    fnats = [info.get("fnat", float("nan")) for info in best_result.values()]

    return (
        sum(dockqs) / n,
        sum(irmsds) / n,
        sum(lrmsds) / n,
        sum(fnats) / n,
        n,
    )


def worker(task):
    target, pipeline, src_type, model_path, native_path = task
    if not os.path.exists(model_path):
        return (target, pipeline, classify_source(src_type), src_type,
                "", "", "", "", 0, "miss_file")
    if not os.path.exists(native_path):
        return (target, pipeline, classify_source(src_type), src_type,
                "", "", "", "", 0, "miss_native")
    try:
        dq, ir, lr, fn, ni = best_dockq_over_chain_perms(model_path, native_path)
        return (target, pipeline, classify_source(src_type), src_type,
                f"{dq:.4f}", f"{ir:.3f}", f"{lr:.3f}", f"{fn:.4f}", ni, "ok")
    except Exception as e:
        msg = str(e).replace("\t", " ").replace("\n", " ")[:200]
        return (target, pipeline, classify_source(src_type), src_type,
                "", "", "", "", 0, f"err:{msg}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--out", default=OUT_TSV)
    ap.add_argument("--target-list", default=TARGET_LIST)
    ap.add_argument("--targets", nargs="+", default=None,
                    help="Override target subset (testing).")
    ap.add_argument("--pipelines", nargs="+", default=["blue", "green"])
    ap.add_argument("--limit", type=int, default=None,
                    help="Only first N targets (testing).")
    args = ap.parse_args()

    if args.targets:
        targets = args.targets
    else:
        with open(args.target_list) as f:
            targets = [t.strip() for t in f if t.strip()]
    if args.limit:
        targets = targets[: args.limit]

    tasks = []
    for t in targets:
        native = os.path.join(CRYSTAL_BASE, f"{t}.pdb")
        for p in args.pipelines:
            for src_type, model_path in expected_inputs(t, p):
                tasks.append((t, p, src_type, model_path, native))

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    print(f"[dockq] {len(tasks)} tasks, {len(targets)} targets, "
          f"pipelines={args.pipelines}, workers={args.workers}", file=sys.stderr)

    header = ("target", "pipeline", "source", "src_type",
              "dockq", "i_rmsd", "l_rmsd", "f_nat", "n_interfaces", "status")

    n_done = n_ok = n_miss = n_err = 0
    with open(args.out, "w") as fout:
        fout.write("\t".join(header) + "\n")
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            for fut in as_completed([ex.submit(worker, t) for t in tasks]):
                row = fut.result()
                fout.write("\t".join(str(x) for x in row) + "\n")
                n_done += 1
                status = row[-1]
                if status == "ok":
                    n_ok += 1
                elif status.startswith("miss"):
                    n_miss += 1
                else:
                    n_err += 1
                if n_done % 500 == 0:
                    print(f"[dockq] {n_done}/{len(tasks)} "
                          f"ok={n_ok} miss={n_miss} err={n_err}",
                          file=sys.stderr, flush=True)
                    fout.flush()

    print(f"[dockq] DONE {n_done}/{len(tasks)} "
          f"ok={n_ok} miss={n_miss} err={n_err} -> {args.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
