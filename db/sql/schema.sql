-- BM5.5 Benchmark Database Schema (v1, 2026-04-20)
-- Canonical source: /data/p_csb_meiler/agarwm5/red_analysis/ on ACCRE
-- Build: python scripts/build_db.py --raw-root $RAW --out bm55.sqlite --snapshot YYYY-MM-DD

PRAGMA foreign_keys = ON;
PRAGMA journal_mode = WAL;

-- =========================================================
-- Dimension tables
-- =========================================================

CREATE TABLE IF NOT EXISTS build_runs (
  snapshot_id      TEXT PRIMARY KEY,             -- e.g., "2026-04-20-abc1234"
  git_commit       TEXT NOT NULL,
  build_time_utc   TEXT NOT NULL,                -- ISO-8601
  raw_manifest_hash TEXT NOT NULL,                -- sha256 of sorted raw TSV paths + mtimes
  code_env_hash    TEXT,                          -- sha256 of requirements + schema.sql
  row_counts_json  TEXT,                          -- JSON of per-table row counts
  qc_status        TEXT NOT NULL CHECK(qc_status IN ('pass','warn','fail'))
);

CREATE TABLE IF NOT EXISTS targets (
  target_id        TEXT PRIMARY KEY,             -- e.g., "1ACB", "BAAD", "4Y7M"
  difficulty       TEXT CHECK(difficulty IN ('rigid','medium','difficult')),
  category         TEXT,                          -- zlab BM5.5: AA/AS/EI/ER/ES/OG/OR/OX
  n_chains         INTEGER,
  n_residues       INTEGER,
  non_standard_flag INTEGER DEFAULT 0,            -- BAAD/BOYV/BP57/CP57 = 1
  parent_pdb_id    TEXT                           -- non-standard parent (e.g. BAAD->3AAD_A:B)
);

CREATE TABLE IF NOT EXISTS pipelines (
  pipeline_id      TEXT PRIMARY KEY CHECK(pipeline_id IN ('blue','green'))
);

CREATE TABLE IF NOT EXISTS sources (
  source_id        TEXT PRIMARY KEY              -- af_relaxed, af_unrelaxed, amber_af, amber_boltz, amber_crystal, boltz, crystal
);

CREATE TABLE IF NOT EXISTS protocols (
  protocol_id      TEXT PRIMARY KEY,             -- cartesian_beta, cartesian_ref15, dualspace_beta, dualspace_ref15, normal_beta, normal_ref15
  score_weights    TEXT NOT NULL CHECK(score_weights IN ('beta_nov16','ref2015')),
  move_set         TEXT NOT NULL CHECK(move_set IN ('cartesian','dualspace','normal'))
);

-- =========================================================
-- Fact tables (wide)
-- =========================================================

-- Rosetta-relaxed structure MolProbity metrics
CREATE TABLE IF NOT EXISTS rosetta_metrics (
  snapshot_id      TEXT NOT NULL,
  target_id        TEXT NOT NULL,
  pipeline_id      TEXT NOT NULL,
  source_id        TEXT NOT NULL,
  src_type         TEXT NOT NULL,                -- raw subtype, e.g., "amber_af_ranked_3"
  protocol_id      TEXT NOT NULL,
  rep              INTEGER NOT NULL CHECK(rep BETWEEN 1 AND 5),
  clashscore       REAL,
  rama_outliers    REAL,
  rama_favored     REAL,
  rota_outliers    REAL,
  molprobity_score REAL,
  cbeta_outliers   INTEGER,
  rms_bonds        REAL,
  rms_angles       REAL,
  PRIMARY KEY (snapshot_id, target_id, pipeline_id, src_type, protocol_id, rep),
  FOREIGN KEY (snapshot_id) REFERENCES build_runs(snapshot_id),
  FOREIGN KEY (target_id)   REFERENCES targets(target_id),
  FOREIGN KEY (pipeline_id) REFERENCES pipelines(pipeline_id),
  FOREIGN KEY (source_id)   REFERENCES sources(source_id),
  FOREIGN KEY (protocol_id) REFERENCES protocols(protocol_id)
);

CREATE INDEX IF NOT EXISTS idx_rm_target_pipeline ON rosetta_metrics(target_id, pipeline_id);
CREATE INDEX IF NOT EXISTS idx_rm_source_protocol ON rosetta_metrics(source_id, protocol_id);
CREATE INDEX IF NOT EXISTS idx_rm_snapshot ON rosetta_metrics(snapshot_id);

-- Pre-Rosetta MolProbity (on input structures: crystal, AF, Boltz, AMBER relaxations)
CREATE TABLE IF NOT EXISTS prerosetta_metrics (
  snapshot_id      TEXT NOT NULL,
  target_id        TEXT NOT NULL,
  pipeline_id      TEXT NOT NULL,
  source_id        TEXT NOT NULL,
  src_type         TEXT NOT NULL,
  clashscore       REAL,
  rama_outliers    REAL,
  rama_favored     REAL,
  rota_outliers    REAL,
  molprobity_score REAL,
  cbeta_outliers   INTEGER,
  rms_bonds        REAL,
  rms_angles       REAL,
  PRIMARY KEY (snapshot_id, target_id, pipeline_id, src_type),
  FOREIGN KEY (snapshot_id) REFERENCES build_runs(snapshot_id),
  FOREIGN KEY (target_id)   REFERENCES targets(target_id),
  FOREIGN KEY (pipeline_id) REFERENCES pipelines(pipeline_id)
);

-- Rosetta-relaxed structure energies (total_score + per-residue energy)
-- NOTE: incomplete vs rosetta_metrics. As of 2026-04-28: 184,352 of 416,340
-- expected rows present (~44% coverage). amber_crystal source has 0 energy
-- rows in both pipelines. Missing cells are logged to qc_quarantine bucket
-- 'coverage_gap'. Re-extraction would require re-running Rosetta scoring on
-- the remaining .pdb.gz outputs in /data/.../{blue,green}_relax/ on ACCRE.
CREATE TABLE IF NOT EXISTS rosetta_energy (
  snapshot_id      TEXT NOT NULL,
  target_id        TEXT NOT NULL,
  pipeline_id      TEXT NOT NULL,
  source_id        TEXT NOT NULL,
  src_type         TEXT NOT NULL,
  protocol_id      TEXT NOT NULL,
  rep              INTEGER NOT NULL CHECK(rep BETWEEN 1 AND 5),
  total_score      REAL,
  per_residue_energy REAL,
  PRIMARY KEY (snapshot_id, target_id, pipeline_id, src_type, protocol_id, rep),
  FOREIGN KEY (snapshot_id) REFERENCES build_runs(snapshot_id),
  FOREIGN KEY (target_id)   REFERENCES targets(target_id),
  FOREIGN KEY (pipeline_id) REFERENCES pipelines(pipeline_id),
  FOREIGN KEY (source_id)   REFERENCES sources(source_id),
  FOREIGN KEY (protocol_id) REFERENCES protocols(protocol_id)
);

CREATE INDEX IF NOT EXISTS idx_re_target_pipeline ON rosetta_energy(target_id, pipeline_id);
CREATE INDEX IF NOT EXISTS idx_re_source_protocol ON rosetta_energy(source_id, protocol_id);
CREATE INDEX IF NOT EXISTS idx_re_snapshot ON rosetta_energy(snapshot_id);

-- TM-score (pre + post Rosetta)
CREATE TABLE IF NOT EXISTS tm_scores (
  snapshot_id      TEXT NOT NULL,
  target_id        TEXT NOT NULL,
  pipeline_id      TEXT NOT NULL,
  source_id        TEXT NOT NULL,
  src_type         TEXT NOT NULL,
  protocol_id      TEXT,                         -- NULL for pre-Rosetta
  rep              INTEGER,                      -- NULL for pre-Rosetta
  tm_score         REAL NOT NULL,
  is_post_rosetta  INTEGER NOT NULL DEFAULT 0 CHECK(is_post_rosetta IN (0,1)),
  PRIMARY KEY (snapshot_id, target_id, pipeline_id, src_type, protocol_id, rep, is_post_rosetta),
  FOREIGN KEY (snapshot_id) REFERENCES build_runs(snapshot_id)
);

CREATE INDEX IF NOT EXISTS idx_tm_target_pipeline ON tm_scores(target_id, pipeline_id);

-- =========================================================
-- Quarantine tables (populated by ingest QC)
-- =========================================================

CREATE TABLE IF NOT EXISTS qc_quarantine (
  snapshot_id      TEXT NOT NULL,
  bucket           TEXT NOT NULL CHECK(bucket IN ('exact_duplicate','conflicting_duplicate','missing_key','orphan_row','invalid_enum','coverage_gap')),
  raw_tsv_path     TEXT,
  raw_row_json     TEXT,                          -- full row as JSON
  notes            TEXT
);

CREATE INDEX IF NOT EXISTS idx_qc_bucket ON qc_quarantine(snapshot_id, bucket);

-- =========================================================
-- Derivative views
-- =========================================================

-- Per-target MolProbity means (equal weighting per target)
CREATE VIEW IF NOT EXISTS v_per_target_mp_means AS
SELECT
  snapshot_id,
  pipeline_id,
  source_id,
  protocol_id,
  target_id,
  AVG(clashscore)       AS clashscore_mean,
  AVG(molprobity_score) AS molprobity_score_mean,
  AVG(rama_favored)     AS rama_favored_mean,
  AVG(rota_outliers)    AS rota_outliers_mean,
  COUNT(*)              AS n_reps
FROM rosetta_metrics
GROUP BY snapshot_id, pipeline_id, source_id, protocol_id, target_id;

-- Cell-level summary (aggregated across targets)
CREATE VIEW IF NOT EXISTS v_cell_summary AS
SELECT
  snapshot_id,
  pipeline_id,
  source_id,
  protocol_id,
  AVG(clashscore_mean)       AS clashscore,
  AVG(molprobity_score_mean) AS molprobity_score,
  AVG(rama_favored_mean)     AS rama_favored,
  AVG(rota_outliers_mean)    AS rota_outliers,
  COUNT(DISTINCT target_id)  AS n_targets
FROM v_per_target_mp_means
GROUP BY snapshot_id, pipeline_id, source_id, protocol_id;
