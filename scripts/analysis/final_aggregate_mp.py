#!/usr/bin/env python3
"""
final_aggregate_mp.py: Final aggregation of combined_rosetta_molprobity.tsv.

Run this on ACCRE after red_ros_mp (10175020) completes:
    conda activate red_analysis
    python final_aggregate_mp.py \
        --tsv /data/p_csb_meiler/agarwm5/red_analysis/combined_rosetta_molprobity.tsv \
        --out /data/p_csb_meiler/agarwm5/red_analysis/final_aggregate_2026-04-18.json

Outputs:
    1. Per-protocol summary (Table 1 values)
    2. Per-source-protocol summary (for AMBER+Rosetta finding in Results 3.3)
    3. Confidence intervals (95%, t-distribution)
    4. Count of targets per protocol/source combination
    5. NaN row count and list of affected targets
    6. Protocol ranking confirmation

Unit of analysis: per-target means (each target contributes equally; no row-mean bias).
Follows methodology from interim_rosetta_analysis_2026-04-18.md.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def load_tsv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    print(f"Loaded {len(df):,} rows from {path}", file=sys.stderr)
    return df


def check_nans(df: pd.DataFrame) -> dict:
    nan_rows = df[df["mp_score"].isna()]
    result = {
        "nan_count": int(len(nan_rows)),
        "nan_fraction": float(len(nan_rows) / len(df)) if len(df) > 0 else 0.0,
        "affected_targets": sorted(nan_rows["pdb_id"].unique().tolist())
        if "pdb_id" in nan_rows.columns
        else [],
    }
    return result


def per_target_mean(group: pd.DataFrame, metric: str) -> float:
    """Mean of per-target means for a metric (correct unit of analysis)."""
    target_means = group.groupby("pdb_id")[metric].mean()
    return float(target_means.mean())


def confidence_interval_95(group: pd.DataFrame, metric: str) -> tuple[float, float]:
    """95% CI using t-distribution on per-target means."""
    target_means = group.groupby("pdb_id")[metric].mean()
    n = len(target_means)
    if n < 2:
        return float("nan"), float("nan")
    mean = target_means.mean()
    se = target_means.sem()
    t_crit = stats.t.ppf(0.975, df=n - 1)
    return float(mean - t_crit * se), float(mean + t_crit * se)


def n_targets(group: pd.DataFrame) -> int:
    return int(group["pdb_id"].nunique())


def protocol_summary(df: pd.DataFrame) -> list[dict]:
    """Per-protocol summary (Table 1 values). All sources combined."""
    metrics = ["mp_score", "clashscore", "rama_favored", "rota_outlier_pct"]
    rows = []
    for protocol, grp in df.groupby("protocol"):
        row = {"protocol": protocol, "n_targets": n_targets(grp)}
        for m in metrics:
            if m not in grp.columns:
                continue
            clean = grp[grp[m].notna()]
            row[f"{m}_mean"] = per_target_mean(clean, m)
            lo, hi = confidence_interval_95(clean, m)
            row[f"{m}_ci_lo"] = lo
            row[f"{m}_ci_hi"] = hi
        rows.append(row)
    rows.sort(key=lambda r: r.get("mp_score_mean", 999))
    return rows


def source_protocol_summary(df: pd.DataFrame) -> list[dict]:
    """Per-source-protocol summary (for AMBER+Rosetta finding)."""
    metrics = ["mp_score", "clashscore", "rama_favored", "rota_outlier_pct"]
    rows = []
    for (source, protocol), grp in df.groupby(["source_type", "protocol"]):
        row = {
            "source": source,
            "protocol": protocol,
            "n_targets": n_targets(grp),
        }
        for m in metrics:
            if m not in grp.columns:
                continue
            clean = grp[grp[m].notna()]
            row[f"{m}_mean"] = per_target_mean(clean, m)
            lo, hi = confidence_interval_95(clean, m)
            row[f"{m}_ci_lo"] = lo
            row[f"{m}_ci_hi"] = hi
        rows.append(row)
    rows.sort(key=lambda r: r.get("mp_score_mean", 999))
    return rows


def pipeline_summary(df: pd.DataFrame) -> list[dict]:
    """Per-pipeline (blue/green) per-protocol summary for concordance (s3.6)."""
    metrics = ["mp_score", "clashscore"]
    rows = []
    if "pipeline" not in df.columns:
        return rows
    for (pipeline, protocol), grp in df.groupby(["pipeline", "protocol"]):
        row = {
            "pipeline": pipeline,
            "protocol": protocol,
            "n_targets": n_targets(grp),
        }
        for m in metrics:
            if m not in grp.columns:
                continue
            clean = grp[grp[m].notna()]
            row[f"{m}_mean"] = per_target_mean(clean, m)
        rows.append(row)
    return rows


def coverage_check(df: pd.DataFrame, expected_targets: int = 257) -> dict:
    """Confirm all expected targets have data."""
    targets_present = sorted(df["pdb_id"].unique().tolist()) if "pdb_id" in df.columns else []
    return {
        "n_targets_present": len(targets_present),
        "expected_targets": expected_targets,
        "coverage_complete": len(targets_present) >= expected_targets,
        "targets_present": targets_present,
    }


def main():
    parser = argparse.ArgumentParser(description="Final aggregation of Rosetta MolProbity data.")
    parser.add_argument("--tsv", required=True, help="Path to combined_rosetta_molprobity.tsv")
    parser.add_argument("--out", required=True, help="Output JSON path")
    parser.add_argument(
        "--expected-targets",
        type=int,
        default=257,
        help="Expected number of unique PDB targets",
    )
    args = parser.parse_args()

    if not Path(args.tsv).exists():
        print(f"ERROR: input file not found: {args.tsv}", file=sys.stderr)
        sys.exit(1)

    df = load_tsv(args.tsv)

    # Normalize column names to lowercase with underscores
    df.columns = [c.strip().lower().replace(" ", "_").replace("-", "_") for c in df.columns]
    print(f"Columns: {list(df.columns)}", file=sys.stderr)

    nan_info = check_nans(df)
    print(
        f"NaN rows: {nan_info['nan_count']} ({nan_info['nan_fraction']:.4%}) "
        f"in targets: {nan_info['affected_targets']}",
        file=sys.stderr,
    )

    coverage = coverage_check(df, args.expected_targets)
    print(
        f"Coverage: {coverage['n_targets_present']}/{coverage['expected_targets']} targets",
        file=sys.stderr,
    )

    proto_summary = protocol_summary(df)
    src_proto_summary = source_protocol_summary(df)
    pipe_summary = pipeline_summary(df)

    # Print top-line Table 1 results to stdout
    print("\n=== TABLE 1: Per-protocol summary (per-target means) ===")
    print(f"{'Protocol':<20} {'MP':>7} {'CI_lo':>7} {'CI_hi':>7} {'Clash':>7} {'n':>5}")
    for row in proto_summary:
        mp = row.get("mp_score_mean", float("nan"))
        lo = row.get("mp_score_ci_lo", float("nan"))
        hi = row.get("mp_score_ci_hi", float("nan"))
        clash = row.get("clashscore_mean", float("nan"))
        n = row.get("n_targets", 0)
        print(f"{row['protocol']:<20} {mp:>7.3f} {lo:>7.3f} {hi:>7.3f} {clash:>7.3f} {n:>5}")

    print("\n=== PER-SOURCE AMBER+ROSETTA FINDING (s3.3 verification) ===")
    print(f"{'Source':<20} {'Protocol':<20} {'MP':>7} {'CI_lo':>7} {'CI_hi':>7} {'n':>5}")
    for row in src_proto_summary:
        if row.get("source") not in ("amber_boltz", "amber_af", "af_relaxed"):
            continue
        mp = row.get("mp_score_mean", float("nan"))
        lo = row.get("mp_score_ci_lo", float("nan"))
        hi = row.get("mp_score_ci_hi", float("nan"))
        n = row.get("n_targets", 0)
        print(f"{row['source']:<20} {row['protocol']:<20} {mp:>7.3f} {lo:>7.3f} {hi:>7.3f} {n:>5}")

    output = {
        "generated": "2026-04-18",
        "total_rows": int(len(df)),
        "nan_info": nan_info,
        "coverage": coverage,
        "protocol_summary": proto_summary,
        "source_protocol_summary": src_proto_summary,
        "pipeline_summary": pipe_summary,
    }

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(output, fh, indent=2, allow_nan=True)
    print(f"\nOutput written to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
