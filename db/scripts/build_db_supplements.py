#!/usr/bin/env python3
"""
build_db_supplements.py: additively populate prerosetta_metrics, tm_scores,
and targets metadata into an existing locked bm55.sqlite snapshot.

Companion to build_db.py. The locked snapshot 2026-04-27a covers rosetta_metrics
only (416,340 rows). This script fills the schema-reserved tables that build_db.py
does not yet ingest, under the SAME snapshot_id, so the manuscript citation point
("snapshot 2026-04-27a") stays singular and the lock remains scientifically
authoritative.

Usage:
    python build_db_supplements.py \\
        --db          /data/.../bm55.sqlite \\
        --raw-root    /data/.../metrics \\
        --difficulty-csv  /path/to/bm55_difficulty.csv \\
        --chains-csv      /path/to/bm55_chains.csv \\
        --snapshot-id 2026-04-27a

Phases:
    B1. prerosetta_metrics  <- combined_molprobity.tsv          (13,344 rows expected)
    B2a. tm_scores (pre)    <- combined_tmscore.tsv             (12,065 rows expected)
    B2b. tm_scores (post)   <- combined_rosetta_tmscore.tsv     (92,700 rows expected)
    B3. targets metadata    <- difficulty-csv + chains-csv      (UPDATE in place)
    B4. build_runs UPDATE   <- recompute manifest hash + row counts JSON
    B5. QC                  <- FK integrity + coverage checks

Idempotent: re-running is safe. INSERT OR IGNORE on prerosetta_metrics + tm_scores;
UPDATE for targets and build_runs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sqlite3
import sys
from pathlib import Path


SOURCE_TO_SRC_TYPE = {
    "af_relaxed":    lambda target, idx: f"af_relaxed_ranked_{idx}",
    "af_unrelaxed":  lambda target, idx: f"af_unrelaxed_ranked_{idx}",
    "amber_af":      lambda target, idx: f"amber_af_ranked_{idx}",
    "amber_boltz":   lambda target, idx: f"amber_boltz_model_{idx}",
    "boltz":         lambda target, idx: f"boltz_boltz_input_model_{idx}",
    "amber_crystal": lambda target, idx: "amber_crystal_relaxed",
    "crystal":       lambda target, idx: f"crystal_{target}",
}

# Upstream bug: combined_molprobity.tsv and combined_rosetta_tmscore.tsv
# write target_id `1E96` as `1e+96` (scientific-notation coercion at TSV
# generation time, likely pandas without dtype=str). Per-target TSVs
# (rosetta_molprobity_*.tsv) preserve the string correctly. Until upstream
# is fixed, normalize on ingest.
TARGET_ALIASES = {
    "1e+96": "1E96",
}


def _normalize_target(target: str) -> str:
    return TARGET_ALIASES.get(target, target)


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


def sha256_iter(items) -> str:
    h = hashlib.sha256()
    for x in items:
        h.update(str(x).encode())
        h.update(b"\n")
    return h.hexdigest()


def compute_supplemental_manifest_hash(raw_root: Path) -> str:
    """Hash including the original rosetta MP TSVs + the supplemental TSVs."""
    patterns = [
        "rosetta_molprobity_*.tsv",
        "combined_molprobity.tsv",
        "combined_tmscore.tsv",
        "combined_rosetta_tmscore.tsv",
    ]
    files = []
    for pat in patterns:
        files.extend(sorted(raw_root.glob(pat)))
    items = [
        (str(p.relative_to(raw_root)), p.stat().st_size, int(p.stat().st_mtime))
        for p in files
    ]
    return sha256_iter(items)


def ingest_prerosetta(conn: sqlite3.Connection, tsv: Path, snapshot_id: str) -> int:
    """B1. Load combined_molprobity.tsv into prerosetta_metrics."""
    cur = conn.cursor()
    inserted = 0
    skipped_unknown_source = 0
    seen: set[tuple] = set()
    with open(tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            source = row["source"]
            if source not in SOURCE_TO_SRC_TYPE:
                skipped_unknown_source += 1
                continue
            target = _normalize_target(row["target"])
            pipeline = row["pipeline"]
            try:
                model_idx = int(row.get("model_idx", "0") or "0")
            except (TypeError, ValueError):
                model_idx = 0
            src_type = SOURCE_TO_SRC_TYPE[source](target, model_idx)
            key = (snapshot_id, target, pipeline, src_type)
            if key in seen:
                continue
            seen.add(key)
            cur.execute(
                """
                INSERT OR IGNORE INTO prerosetta_metrics
                (snapshot_id, target_id, pipeline_id, source_id, src_type,
                 clashscore, rama_outliers, rama_favored, rota_outliers,
                 molprobity_score, cbeta_outliers, rms_bonds, rms_angles)
                VALUES (?,?,?,?,?, ?,?,?,?,?, ?,?,?)
                """,
                (
                    snapshot_id, target, pipeline, source, src_type,
                    _f(row.get("clashscore")),
                    _f(row.get("rama_outliers")),
                    _f(row.get("rama_favored")),
                    _f(row.get("rota_outliers")),
                    _f(row.get("molprobity_score")),
                    _i(row.get("cbeta_outliers")),
                    _f(row.get("rms_bonds")),
                    _f(row.get("rms_angles")),
                ),
            )
            inserted += cur.rowcount  # 1 if inserted, 0 if ignored
    conn.commit()
    if skipped_unknown_source:
        print(f"  warn: skipped {skipped_unknown_source} rows with unrecognized source", file=sys.stderr)
    return inserted


def ingest_pre_tmscore(conn: sqlite3.Connection, tsv: Path, snapshot_id: str) -> int:
    """B2a. Load combined_tmscore.tsv (pre-Rosetta) into tm_scores.

    Idempotency note: tm_scores PK includes (protocol_id, rep) which are NULL
    for pre-Rosetta rows. SQLite treats NULL as distinct in UNIQUE indexes,
    so INSERT OR IGNORE will NOT dedup on re-run. Wipe-then-load for the
    pre-Rosetta slice of this snapshot to keep re-runs idempotent.
    """
    cur = conn.cursor()
    cur.execute(
        "DELETE FROM tm_scores WHERE snapshot_id=? AND is_post_rosetta=0",
        (snapshot_id,),
    )
    inserted = 0
    skipped_unknown_source = 0
    skipped_missing_tm = 0
    seen: set[tuple] = set()
    with open(tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            source = row["source"]
            if source not in SOURCE_TO_SRC_TYPE:
                skipped_unknown_source += 1
                continue
            target = _normalize_target(row["target"])
            pipeline = row["pipeline"]
            try:
                model_idx = int(row.get("model_idx", "0") or "0")
            except (TypeError, ValueError):
                model_idx = 0
            src_type = SOURCE_TO_SRC_TYPE[source](target, model_idx)
            tm = _f(row.get("tmscore"))
            if tm is None:
                skipped_missing_tm += 1
                continue
            key = (snapshot_id, target, pipeline, src_type, None, None, 0)
            if key in seen:
                continue
            seen.add(key)
            cur.execute(
                """
                INSERT OR IGNORE INTO tm_scores
                (snapshot_id, target_id, pipeline_id, source_id, src_type,
                 protocol_id, rep, tm_score, is_post_rosetta)
                VALUES (?,?,?,?,?, NULL, NULL, ?, 0)
                """,
                (snapshot_id, target, pipeline, source, src_type, tm),
            )
            inserted += cur.rowcount
    conn.commit()
    if skipped_unknown_source:
        print(f"  warn: skipped {skipped_unknown_source} pre-Rosetta TM rows with unrecognized source", file=sys.stderr)
    if skipped_missing_tm:
        print(f"  warn: skipped {skipped_missing_tm} pre-Rosetta TM rows with missing tmscore", file=sys.stderr)
    return inserted


def ingest_post_tmscore(conn: sqlite3.Connection, tsv: Path, snapshot_id: str) -> int:
    """B2b. Load combined_rosetta_tmscore.tsv (post-Rosetta) into tm_scores."""
    cur = conn.cursor()
    inserted = 0
    skipped_missing_tm = 0
    seen: set[tuple] = set()
    with open(tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            target = _normalize_target(row["target"])
            pipeline = row["pipeline"]
            src_type = row["src_type"]
            source = row.get("source") or ""
            protocol_id = row["protocol"]
            try:
                rep = int(row["rep"])
            except (TypeError, ValueError):
                continue
            tm = _f(row.get("tmscore"))
            if tm is None:
                skipped_missing_tm += 1
                continue
            key = (snapshot_id, target, pipeline, src_type, protocol_id, rep, 1)
            if key in seen:
                continue
            seen.add(key)
            cur.execute(
                """
                INSERT OR IGNORE INTO tm_scores
                (snapshot_id, target_id, pipeline_id, source_id, src_type,
                 protocol_id, rep, tm_score, is_post_rosetta)
                VALUES (?,?,?,?,?, ?,?, ?, 1)
                """,
                (snapshot_id, target, pipeline, source, src_type, protocol_id, rep, tm),
            )
            inserted += cur.rowcount
    conn.commit()
    if skipped_missing_tm:
        print(f"  warn: skipped {skipped_missing_tm} post-Rosetta TM rows with missing tmscore", file=sys.stderr)
    return inserted


def stage_c_backfill_blue_crystal(conn: sqlite3.Connection, snapshot_id: str) -> tuple[int, list[str]]:
    """
    Stage C: backfill the 20 Blue-pipeline crystal pre-Rosetta MolProbity rows
    that were missing from combined_molprobity.tsv.

    Provenance: Blue and Green pipelines used byte-identical cleaned crystal
    PDBs (verified empirically: 1A2K, 7CEI, 1ACB show identical clashscore /
    MolProbity values across the two pipelines for the crystal source). The
    Blue per-target MolProbity TSV for 1ACB (molprobity_blue_1ACB.tsv) carries
    a crystal row whose values match Green's exactly. The remaining 19 Blue
    crystal rows were never aggregated into combined_molprobity.tsv.

    Resolution: copy Green crystal rows to Blue for the 20 missing target_ids.
    Justified because MolProbity is deterministic on identical input PDBs,
    and the input PDBs were verified identical across pipelines.
    """
    missing = [
        "1ACB", "1ATN", "1CGI", "1FC2", "1MAH", "1MLC", "1PPE", "1RLB",
        "1UDI", "1VFB", "2BTF", "2CFH", "2G77", "2GAF", "2MTA", "2PCC",
        "2SIC", "2SNI", "3FN1", "3R9A",
    ]
    cur = conn.cursor()
    placeholders = ",".join("?" * len(missing))
    cur.execute(
        f"""
        INSERT OR IGNORE INTO prerosetta_metrics
        (snapshot_id, target_id, pipeline_id, source_id, src_type,
         clashscore, rama_outliers, rama_favored, rota_outliers,
         molprobity_score, cbeta_outliers, rms_bonds, rms_angles)
        SELECT ?, target_id, 'blue', source_id, src_type,
               clashscore, rama_outliers, rama_favored, rota_outliers,
               molprobity_score, cbeta_outliers, rms_bonds, rms_angles
        FROM prerosetta_metrics
        WHERE snapshot_id = ?
          AND pipeline_id = 'green'
          AND source_id = 'crystal'
          AND target_id IN ({placeholders})
        """,
        (snapshot_id, snapshot_id, *missing),
    )
    inserted = cur.rowcount
    conn.commit()

    # Verify what we backfilled
    cur.execute(
        f"""
        SELECT target_id FROM prerosetta_metrics
        WHERE snapshot_id = ? AND pipeline_id = 'blue' AND source_id = 'crystal'
          AND target_id IN ({placeholders})
        """,
        (snapshot_id, *missing),
    )
    backfilled = sorted(r[0] for r in cur.fetchall())
    still_missing = sorted(set(missing) - set(backfilled))
    return inserted, still_missing


def update_targets_metadata(
    conn: sqlite3.Connection,
    difficulty_csv: Path,
    chains_csv: Path,
) -> tuple[int, int]:
    """B3. Populate targets.{difficulty, category, n_chains, n_residues}."""
    cur = conn.cursor()

    # Load difficulty + category
    diff_map: dict[str, tuple[str, str]] = {}
    with open(difficulty_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            diff_map[row["target_id"]] = (row.get("difficulty") or None,
                                          row.get("category") or None)

    # Load n_chains + n_residues
    chains_map: dict[str, tuple[int, int]] = {}
    with open(chains_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            chains_map[row["target_id"]] = (int(row["n_chains"]),
                                            int(row["n_residues"]))

    db_targets = {r[0] for r in cur.execute("SELECT target_id FROM targets").fetchall()}

    diff_updates = 0
    chains_updates = 0
    for target in sorted(db_targets):
        if target in diff_map:
            difficulty, category = diff_map[target]
            cur.execute(
                "UPDATE targets SET difficulty=?, category=? WHERE target_id=?",
                (difficulty, category, target),
            )
            diff_updates += cur.rowcount
        if target in chains_map:
            n_chains, n_residues = chains_map[target]
            cur.execute(
                "UPDATE targets SET n_chains=?, n_residues=? WHERE target_id=?",
                (n_chains, n_residues, target),
            )
            chains_updates += cur.rowcount

    missing_diff = sorted(db_targets - set(diff_map.keys()))
    missing_chains = sorted(db_targets - set(chains_map.keys()))
    if missing_diff:
        print(f"  warn: {len(missing_diff)} targets missing difficulty/category: {missing_diff[:5]}{'...' if len(missing_diff) > 5 else ''}", file=sys.stderr)
    if missing_chains:
        print(f"  warn: {len(missing_chains)} targets missing chains/residues: {missing_chains[:5]}{'...' if len(missing_chains) > 5 else ''}", file=sys.stderr)

    conn.commit()
    return diff_updates, chains_updates


def update_build_run(
    conn: sqlite3.Connection,
    snapshot_id: str,
    raw_root: Path,
    new_counts: dict,
    qc_status: str,
) -> None:
    """B4. UPDATE the existing build_runs row in place. Preserve git_commit, build_time_utc."""
    cur = conn.cursor()
    cur.execute(
        "SELECT row_counts_json FROM build_runs WHERE snapshot_id = ?",
        (snapshot_id,),
    )
    existing_json = cur.fetchone()
    if existing_json is None or existing_json[0] is None:
        merged = new_counts
    else:
        try:
            existing = json.loads(existing_json[0])
        except (TypeError, ValueError):
            existing = {}
        merged = {**existing, **new_counts}

    new_manifest_hash = compute_supplemental_manifest_hash(raw_root)
    cur.execute(
        """
        UPDATE build_runs
        SET row_counts_json   = ?,
            raw_manifest_hash = ?,
            qc_status         = ?
        WHERE snapshot_id = ?
        """,
        (json.dumps(merged), new_manifest_hash, qc_status, snapshot_id),
    )
    conn.commit()


def qc_supplements(conn: sqlite3.Connection, snapshot_id: str) -> dict:
    """B5. Coverage + FK integrity checks for supplementary tables."""
    cur = conn.cursor()

    n_pre  = cur.execute("SELECT COUNT(*) FROM prerosetta_metrics WHERE snapshot_id=?", (snapshot_id,)).fetchone()[0]
    n_pre_tm  = cur.execute("SELECT COUNT(*) FROM tm_scores WHERE snapshot_id=? AND is_post_rosetta=0", (snapshot_id,)).fetchone()[0]
    n_post_tm = cur.execute("SELECT COUNT(*) FROM tm_scores WHERE snapshot_id=? AND is_post_rosetta=1", (snapshot_id,)).fetchone()[0]
    n_targets_with_diff = cur.execute("SELECT COUNT(*) FROM targets WHERE difficulty IS NOT NULL").fetchone()[0]
    n_targets_with_chains = cur.execute("SELECT COUNT(*) FROM targets WHERE n_chains IS NOT NULL AND n_residues IS NOT NULL").fetchone()[0]
    n_targets = cur.execute("SELECT COUNT(*) FROM targets").fetchone()[0]

    # FK integrity: prerosetta_metrics target_ids exist in targets
    orphans_pre = cur.execute("""
        SELECT COUNT(*) FROM prerosetta_metrics pm
        WHERE NOT EXISTS (SELECT 1 FROM targets t WHERE t.target_id = pm.target_id)
    """).fetchone()[0]
    orphans_tm = cur.execute("""
        SELECT COUNT(*) FROM tm_scores tm
        WHERE NOT EXISTS (SELECT 1 FROM targets t WHERE t.target_id = tm.target_id)
    """).fetchone()[0]

    return {
        "prerosetta_metrics_rows": n_pre,
        "tm_scores_pre_rows": n_pre_tm,
        "tm_scores_post_rows": n_post_tm,
        "tm_scores_total": n_pre_tm + n_post_tm,
        "targets_with_difficulty": n_targets_with_diff,
        "targets_with_chains_residues": n_targets_with_chains,
        "targets_total": n_targets,
        "orphan_prerosetta_target_rows": orphans_pre,
        "orphan_tm_score_target_rows": orphans_tm,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True, type=Path)
    ap.add_argument("--raw-root", required=True, type=Path,
                    help="Directory containing combined_molprobity.tsv etc")
    ap.add_argument("--difficulty-csv", required=True, type=Path)
    ap.add_argument("--chains-csv", required=True, type=Path)
    ap.add_argument("--snapshot-id", default="2026-04-27a")
    args = ap.parse_args()

    if not args.db.exists():
        sys.exit(f"DB not found: {args.db}")

    pre_mp_tsv = args.raw_root / "combined_molprobity.tsv"
    pre_tm_tsv = args.raw_root / "combined_tmscore.tsv"
    post_tm_tsv = args.raw_root / "combined_rosetta_tmscore.tsv"
    for p in [pre_mp_tsv, pre_tm_tsv, post_tm_tsv, args.difficulty_csv, args.chains_csv]:
        if not p.exists():
            sys.exit(f"missing input: {p}")

    conn = sqlite3.connect(str(args.db))
    conn.execute("PRAGMA foreign_keys = ON")

    print(f"Supplementing {args.db} under snapshot {args.snapshot_id}")

    print("B1. ingesting prerosetta_metrics ...")
    n1 = ingest_prerosetta(conn, pre_mp_tsv, args.snapshot_id)
    print(f"   inserted {n1} rows")

    print("B2a. ingesting pre-Rosetta tm_scores ...")
    n2a = ingest_pre_tmscore(conn, pre_tm_tsv, args.snapshot_id)
    print(f"   inserted {n2a} rows")

    print("B2b. ingesting post-Rosetta tm_scores ...")
    n2b = ingest_post_tmscore(conn, post_tm_tsv, args.snapshot_id)
    print(f"   inserted {n2b} rows")

    print("Stage C. backfilling 20 Blue crystal MP rows from Green ...")
    c_inserted, c_missing = stage_c_backfill_blue_crystal(conn, args.snapshot_id)
    print(f"   inserted {c_inserted} rows; still_missing={c_missing}")

    print("B3. updating targets metadata ...")
    diff_n, chains_n = update_targets_metadata(conn, args.difficulty_csv, args.chains_csv)
    print(f"   updated {diff_n} difficulty/category, {chains_n} chains/residues")

    print("B5. QC ...")
    qc = qc_supplements(conn, args.snapshot_id)
    for k, v in qc.items():
        print(f"   {k}: {v}")

    qc_pass = (
        qc["orphan_prerosetta_target_rows"] == 0 and
        qc["orphan_tm_score_target_rows"] == 0 and
        qc["targets_with_difficulty"] == qc["targets_total"] and
        qc["targets_with_chains_residues"] == qc["targets_total"] and
        qc["prerosetta_metrics_rows"] >= 13_364 and          # 13,344 from TSV + 20 Stage C backfill
        qc["tm_scores_pre_rows"] >= 12_000 and
        qc["tm_scores_post_rows"] >= 92_000
    )
    qc_status = "pass" if qc_pass else "warn"

    print("B4. updating build_runs row in place ...")
    new_counts = {
        "rosetta_metrics": 416340,
        "prerosetta_metrics": qc["prerosetta_metrics_rows"],
        "tm_scores_pre": qc["tm_scores_pre_rows"],
        "tm_scores_post": qc["tm_scores_post_rows"],
        "tm_scores_total": qc["tm_scores_total"],
        "targets_with_metadata": qc["targets_with_difficulty"],
    }
    update_build_run(conn, args.snapshot_id, args.raw_root, new_counts, qc_status)

    final = conn.execute(
        "SELECT snapshot_id, qc_status, raw_manifest_hash, row_counts_json FROM build_runs WHERE snapshot_id=?",
        (args.snapshot_id,),
    ).fetchone()
    print(f"build_runs[{final[0]}]:")
    print(f"   qc_status         = {final[1]}")
    print(f"   raw_manifest_hash = {final[2]}")
    print(f"   row_counts_json   = {final[3]}")
    print(f"final qc_status = {qc_status}")

    conn.close()


if __name__ == "__main__":
    main()
