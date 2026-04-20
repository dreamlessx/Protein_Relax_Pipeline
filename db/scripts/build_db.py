#!/usr/bin/env python3
"""
build_db.py: reproducibly build bm55.sqlite from raw TSVs.

Usage:
    python build_db.py --raw-root /data/p_csb_meiler/agarwm5/red_analysis/metrics \
                       --out       /path/to/bm55.sqlite \
                       --snapshot  2026-04-20

Phases:
    1. ingest:   stream raw TSVs into staging tables (text cols)
    2. validate: QC buckets, quarantine bad rows, fail build if exact + missing > threshold
    3. publish:  typed insert into canonical tables, build views, write manifest

Raw TSV conventions (ACCRE paths):
    metrics/rosetta_molprobity_{pipeline}_{target}.tsv   -- 810 rows each
    metrics/molprobity_{pipeline}_{target}.tsv           -- pre-Rosetta MP
    metrics/rosetta_tmscore_{pipeline}_{target}.tsv      -- Rosetta TM
    metrics/tmscore_{pipeline}_{target}.tsv              -- pre-Rosetta TM
"""
from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import os
import sqlite3
import sys
import time
from pathlib import Path


SCHEMA_PATH = Path(__file__).parent.parent / "sql" / "schema.sql"

SOURCE_BUCKETS = {
    "af_relaxed":    lambda s: s.startswith("af_relaxed"),
    "af_unrelaxed":  lambda s: s.startswith("af_unrelaxed"),
    "amber_af":      lambda s: s.startswith("amber_af_ranked"),
    "amber_boltz":   lambda s: s.startswith("amber_boltz_model"),
    "amber_crystal": lambda s: s.startswith("amber_crystal"),
    "boltz":         lambda s: s.startswith("boltz_boltz_input"),
    "crystal":       lambda s: s.startswith("crystal_"),
}

LEGACY_SRC = {"amber_af_relaxed", "amber_boltz_relaxed"}

PROTOCOLS = {
    "cartesian_beta":   ("beta_nov16", "cartesian"),
    "cartesian_ref15":  ("ref2015",    "cartesian"),
    "dualspace_beta":   ("beta_nov16", "dualspace"),
    "dualspace_ref15":  ("ref2015",    "dualspace"),
    "normal_beta":      ("beta_nov16", "normal"),
    "normal_ref15":     ("ref2015",    "normal"),
}


def bucket(src_type: str) -> str | None:
    if src_type in LEGACY_SRC:
        return None
    for bkt, pred in SOURCE_BUCKETS.items():
        if pred(src_type):
            return bkt
    return None


def sha256_iter(items) -> str:
    h = hashlib.sha256()
    for x in items:
        h.update(str(x).encode())
        h.update(b"\n")
    return h.hexdigest()


def compute_raw_manifest_hash(raw_root: Path) -> str:
    tsvs = sorted(raw_root.glob("rosetta_molprobity_*.tsv"))
    items = [(str(p.relative_to(raw_root)), p.stat().st_size, int(p.stat().st_mtime))
             for p in tsvs]
    return sha256_iter(items)


def init_db(conn: sqlite3.Connection) -> None:
    with open(SCHEMA_PATH) as f:
        conn.executescript(f.read())
    # Seed dimension tables
    conn.executemany("INSERT OR IGNORE INTO pipelines(pipeline_id) VALUES(?)",
                     [("blue",), ("green",)])
    conn.executemany("INSERT OR IGNORE INTO sources(source_id) VALUES(?)",
                     [(s,) for s in SOURCE_BUCKETS.keys()])
    conn.executemany("INSERT OR IGNORE INTO protocols(protocol_id, score_weights, move_set) VALUES(?,?,?)",
                     [(p, sw, mv) for p, (sw, mv) in PROTOCOLS.items()])
    conn.commit()


def seed_targets(conn: sqlite3.Connection, target_list_path: Path) -> None:
    """Seed targets table. Difficulty left NULL for now (can be loaded from BM5.5 metadata later)."""
    with open(target_list_path) as f:
        targets = [t.strip() for t in f.read().splitlines() if t.strip()]
    non_std = {"BAAD", "BOYV", "BP57", "CP57"}
    rows = [(t, None, None, None, None, 1 if t in non_std else 0) for t in targets]
    conn.executemany(
        "INSERT OR IGNORE INTO targets(target_id, difficulty, category, n_chains, n_residues, non_standard_flag) VALUES(?,?,?,?,?,?)",
        rows,
    )
    conn.commit()


def ingest_rosetta_metrics(conn: sqlite3.Connection, raw_root: Path, snapshot_id: str) -> tuple[int, int, int]:
    """
    Load rosetta_molprobity_{pipeline}_{target}.tsv files.
    Returns (inserted, exact_dups, legacy_skipped).
    """
    cur = conn.cursor()
    files = sorted(glob.glob(str(raw_root / "rosetta_molprobity_*.tsv")))
    total = 0
    exact_dups = 0
    legacy = 0
    seen: set[tuple] = set()
    for fp in files:
        with open(fp) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                src_type = row["src_type"]
                bkt = bucket(src_type)
                if bkt is None:
                    legacy += 1
                    continue
                proto = row["protocol"]
                if proto not in PROTOCOLS:
                    continue
                try:
                    rep = int(row["rep"])
                except (TypeError, ValueError):
                    continue
                key = (snapshot_id, row["target"], row["pipeline"], src_type, proto, rep)
                if key in seen:
                    exact_dups += 1
                    continue
                seen.add(key)
                try:
                    cur.execute("""
                        INSERT INTO rosetta_metrics
                        (snapshot_id, target_id, pipeline_id, source_id, src_type, protocol_id, rep,
                         clashscore, rama_outliers, rama_favored, rota_outliers, molprobity_score,
                         cbeta_outliers, rms_bonds, rms_angles)
                        VALUES(?,?,?,?,?,?,?, ?,?,?,?,?, ?,?,?)
                    """, (
                        snapshot_id, row["target"], row["pipeline"], bkt, src_type, proto, rep,
                        _f(row.get("clashscore")), _f(row.get("rama_outliers")),
                        _f(row.get("rama_favored")), _f(row.get("rota_outliers")),
                        _f(row.get("molprobity_score")), _i(row.get("cbeta_outliers")),
                        _f(row.get("rms_bonds")), _f(row.get("rms_angles")),
                    ))
                    total += 1
                except sqlite3.Error as e:
                    cur.execute("""
                        INSERT INTO qc_quarantine(snapshot_id, bucket, raw_tsv_path, raw_row_json, notes)
                        VALUES(?,?,?,?,?)
                    """, (snapshot_id, "conflicting_duplicate", fp, json.dumps(row), str(e)))
    conn.commit()
    return total, exact_dups, legacy


def _f(v):
    try:
        return float(v) if v not in (None, "", "nan", "NaN") else None
    except (TypeError, ValueError):
        return None


def _i(v):
    try:
        return int(float(v)) if v not in (None, "") else None
    except (TypeError, ValueError):
        return None


def qc_coverage(conn: sqlite3.Connection, snapshot_id: str) -> dict:
    """Compute coverage gaps: cells < 810 expected rows."""
    cur = conn.cursor()
    target_count = cur.execute("SELECT COUNT(*) FROM targets").fetchone()[0]
    expected_total = target_count * 2 * 27 * 6 * 5   # target * pipeline * src * protocol * rep
    actual = cur.execute(
        "SELECT COUNT(*) FROM rosetta_metrics WHERE snapshot_id = ?",
        (snapshot_id,),
    ).fetchone()[0]
    gaps = cur.execute("""
        SELECT target_id, pipeline_id, COUNT(*) AS n
        FROM rosetta_metrics
        WHERE snapshot_id = ?
        GROUP BY target_id, pipeline_id
        HAVING n < 810
        ORDER BY n ASC
    """, (snapshot_id,)).fetchall()
    # Quarantine one row per coverage gap cell
    for target_id, pipeline_id, n in gaps:
        cur.execute("""
            INSERT INTO qc_quarantine(snapshot_id, bucket, raw_tsv_path, raw_row_json, notes)
            VALUES(?,?,?,?,?)
        """, (snapshot_id, "coverage_gap", None,
              json.dumps({"target_id": target_id, "pipeline_id": pipeline_id, "n": n}),
              f"{810 - n} rows short"))
    conn.commit()
    return {
        "target_count": target_count,
        "expected_rows": expected_total,
        "actual_rows": actual,
        "coverage_pct": round(100 * actual / expected_total, 3) if expected_total else 0,
        "gap_cells": len(gaps),
        "total_missing_rows": sum(810 - n for _, _, n in gaps),
    }


def publish_build_run(conn: sqlite3.Connection, snapshot_id: str, raw_root: Path,
                       row_counts: dict, qc_status: str) -> None:
    git_commit = os.popen("git -C "
                          f"{Path(__file__).parent.parent.parent} rev-parse --short HEAD").read().strip() or "unknown"
    schema_hash = sha256_iter([SCHEMA_PATH.read_bytes()])
    conn.execute("""
        INSERT OR REPLACE INTO build_runs
        (snapshot_id, git_commit, build_time_utc, raw_manifest_hash, code_env_hash, row_counts_json, qc_status)
        VALUES(?,?,?,?,?,?,?)
    """, (
        snapshot_id, git_commit, time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        compute_raw_manifest_hash(raw_root), schema_hash, json.dumps(row_counts), qc_status,
    ))
    conn.commit()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--snapshot", required=True)
    ap.add_argument("--target-list", type=Path,
                    default=Path("/data/p_csb_meiler/agarwm5/red_analysis/target_list.txt"))
    args = ap.parse_args()

    if args.out.exists():
        args.out.unlink()

    conn = sqlite3.connect(str(args.out))
    print(f"Building {args.out} from {args.raw_root} snapshot {args.snapshot}")
    init_db(conn)
    if args.target_list.exists():
        seed_targets(conn, args.target_list)
        print(f"  seeded {conn.execute('SELECT COUNT(*) FROM targets').fetchone()[0]} targets")

    # Placeholder build_run so rosetta_metrics FK resolves; updated at end
    conn.execute("""
        INSERT OR REPLACE INTO build_runs
        (snapshot_id, git_commit, build_time_utc, raw_manifest_hash, qc_status)
        VALUES(?,?,?,?,?)
    """, (args.snapshot, "pending", "pending", "pending", "warn"))
    conn.commit()

    t0 = time.time()
    inserted, exact_dups, legacy = ingest_rosetta_metrics(conn, args.raw_root, args.snapshot)
    print(f"  ingested {inserted} rosetta_metrics rows ({exact_dups} exact dupes skipped, "
          f"{legacy} legacy-src skipped) in {time.time()-t0:.1f}s")

    cov = qc_coverage(conn, args.snapshot)
    print(f"  coverage: {cov['actual_rows']:,} / {cov['expected_rows']:,} = {cov['coverage_pct']}%")
    print(f"  gap cells (< 810 rows): {cov['gap_cells']}")
    print(f"  total missing rows:    {cov['total_missing_rows']:,}")

    status = "pass" if cov["coverage_pct"] >= 99.0 else "warn" if cov["coverage_pct"] >= 95.0 else "fail"
    row_counts = {"rosetta_metrics": inserted, "exact_dups_skipped": exact_dups,
                  "legacy_skipped": legacy, **cov}
    publish_build_run(conn, args.snapshot, args.raw_root, row_counts, status)
    print(f"  qc_status = {status}")
    print(f"Done: {args.out} ({args.out.stat().st_size / 1024**2:.1f} MB)")
    conn.close()


if __name__ == "__main__":
    main()
