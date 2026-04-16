# Project Status

**Last updated:** 2026-04-16

## Dataset Summary

| Metric | Value |
|--------|-------|
| BM5.5 targets | 257 |
| Total chains | 605 |
| Total residues | 122,966 |
| Rigid-body | 162 |
| Medium difficulty | 60 |
| Difficult | 35 |

## Input Consistency

All 257 FASTAs are crystal-derived and verified identical across all prediction methods.

| Check | Status |
|-------|--------|
| Crystal == AF FASTA | 257/257 |
| AF FASTA == Boltz FASTA | 257/257 |
| No duplicate chains | 257/257 |
| No His-tags | 257/257 |

## Prediction Status

| Method | Progress | Status |
|--------|----------|--------|
| AlphaFold 2.3.2 | 257/257 | COMPLETE |
| AF built-in AMBER | 257/257 | COMPLETE (all have ranked_0..4.pdb) |
| AF unrelaxed | 257/257 | COMPLETE |
| Boltz-1 v0.4.1 | 257/257 | COMPLETE |
| Standalone AMBER (AF) | 257/257 | COMPLETE |
| Standalone AMBER (Boltz) | 257/257 | COMPLETE |
| Standalone AMBER (crystal) | 255/257 | 1ACB, 1ATN blocked (AMBER convergence failure) |

## Canonical Data Matrix

Per complex (both pipelines identical structure):

| Input type | Models | Rosetta outputs (6 protocols x 5 reps) | Files per input |
|------------|--------|----------------------------------------|-----------------|
| crystal | 1 | 30 | 31 |
| amber-crystal | 1 | 30 | 31 |
| AF unrelaxed | 5 | 150 | 155 |
| AF + built-in AMBER | 5 | 150 | 155 |
| AF + standalone AMBER | 5 | 150 | 155 |
| Boltz | 5 | 150 | 155 |
| Boltz + standalone AMBER | 5 | 150 | 155 |
| **Per complex total** | **27** | **810** | **837** |

**Per-pipeline Rosetta target:** 27 x 30 x 257 - 60 (1ACB+1ATN amber-crystal blocked) = **208,050 files**

**Both pipelines combined Rosetta target:** 416,100 .pdb.gz files

**Base structures:** 27 x 257 - 2 = 6,937 per pipeline (13,874 total)

## Rosetta Relaxation Status (2026-04-16, post-legacy-cleanup)

### File counts

| Pipeline | Rosetta .pdb.gz on disk | vs target (208,050) | Gap |
|----------|-------------------------|---------------------|-----|
| Blue | 207,496 | 99.7% | 554 files in partial dirs |
| Green | in progress | TBD | TBD |

### Blue dir status breakdown

| Status | Count |
|--------|-------|
| Full (30 .pdb.gz) | 6,878 |
| Partial (1-29 files) | 57 |
| Empty (0 files) | 1 |
| Missing | 3 (2 blocked: 1ACB/1ATN amber-crystal; 1 in-flight) |
| **Total** | **6,939** |

### Category completeness (Blue)

| Category | Rosetta complete | Notes |
|----------|------------------|-------|
| crystal | 252/257 | 5 partial filling via Job 10145091 |
| amber-crystal | 197/257 -> 255/257 | Job 10107933 (88 targets); 2 permanent blocked |
| AF unrelaxed | 257/257 | COMPLETE |
| AF + standalone AMBER | 257/257 | COMPLETE |
| AF + built-in AMBER | 257/257 | COMPLETE |
| Boltz | 257/257 | COMPLETE |
| Boltz + AMBER | 254/257 | 3 partial filling via Job 10145097 |

### Category completeness (Green)

| Category | Rosetta complete | Notes |
|----------|------------------|-------|
| crystal | 257/257 | COMPLETE |
| amber-crystal | 183/257 -> 255/257 | Job 10107934 (85 targets); 2 permanent blocked |
| AF unrelaxed | 257/257 | COMPLETE |
| AF + standalone AMBER | 257/257 | COMPLETE |
| AF + built-in AMBER | 257/257 | COMPLETE |
| Boltz | 257/257 | COMPLETE |
| Boltz + AMBER | 257/257 | COMPLETE |

### Active SLURM Jobs (2026-04-16)

| Job ID | Name | Purpose | Targets |
|--------|------|---------|---------|
| 10107933 | blue_axtal_v4 | Blue amber-crystal Rosetta | 88 |
| 10107934 | green_axtal_v4 | Green amber-crystal Rosetta | 85 |
| 10145091 | blue_crystal_fin | Blue crystal partial fills | 5 (3WD5, 4DN4, 4G6M, 4GAM, 4Y7M) |
| 10145097 | blue_boltz_amber_fin | Blue Boltz-AMBER partial fills | 3 (1HE8, 1K5D, 4GAM) |

## Cleanup Applied (2026-04-16)

**Legacy collision dirs removed:** 340 dirs, 7,775 .pdb.gz files (Blue only; Green never had them).

The `amber_af_relaxed/` and `amber_boltz_relaxed/` directories were created before the AMBER naming collision fix. They contained stale outputs where all 5 AF-AMBER models (or 5 Boltz-AMBER models) wrote to the same output path, causing collision (last write wins). The authoritative per-model data lives in `amber_af_ranked_{0..4}/` and `amber_boltz_model_{0..4}/`.

These legacy dirs were previously pooled with proper data in `aggregate_rosetta_molprobity.py` and `aggregate_rosetta_tmscore.py`. Aggregate scripts must be re-run after removal to exclude the collision noise from published metrics.

## Blockers (Permanent)

- **1ACB, 1ATN AMBER-crystal convergence failure**
  - AMBER minimization fails after 100 attempts on the crystal PDB
  - Affects both Blue and Green pipelines
  - Impact: -60 Rosetta files + -2 base files per pipeline
  - Final coverage: 255/257 for amber-crystal category
  - Acceptance: paper uses 255/257; other 6 categories remain 257/257

## Pipeline Design

- **Blue pipeline** = primary (Claude). Job prefix: `blue_`
- **Green pipeline** = independent verification of Blue. Job prefix: `green_`
- **Red analysis** = metrics & figures. Job prefix: `red_`
- Both Blue and Green follow identical Rosetta flags, same inputs, same output structure

## Validation

| Validation | Blue | Green |
|------------|------|-------|
| MolProbity (pre-Rosetta) | 257/257 | 257/257 |
| MolProbity (Rosetta outputs) | 174,721 rows aggregated | (in same table) |
| TM-score (pre-Rosetta) | 257/257 | 257/257 |
| Rosetta energy | ~174K scorefiles | (combined) |
| PoseBusters | Not in pipeline (deprecated per paper scope) | — |

## Key Fixes Applied

| # | Fix | Impact |
|---|-----|--------|
| 1 | Crystal-derived FASTAs | 236 targets changed from UniProt full-length; net -3,956 residues |
| 2 | Crystal chain deduplication | 36 PDBs stripped of homo-multimer duplicate chains |
| 3 | FASTA deduplication | 135 Boltz FASTAs had duplicated homo-multimer chains |
| 4 | His-tag removal | 41 targets had expression artifacts removed |
| 5 | 1KTZ template bug | AF workaround using max_template_date=1900 |
| 6 | Rosetta AMBER naming collision | MODEL_LABELS array fix to avoid overwriting |
| 7 | Legacy collision dir purge (2026-04-16) | 340 dirs / 7,775 stale files removed for clean per-model data |

## Remaining Work

- [ ] In-flight Rosetta fills drain (~5-6 days wall clock)
- [ ] Re-run `aggregate_rosetta_molprobity.py` after legacy removal
- [ ] Re-run `aggregate_rosetta_tmscore.py` after legacy removal
- [ ] PDB integrity pass (file size + energy outlier check)
- [ ] Finalize paper tables with clean per-model aggregates

## Timeline

| Date | Event |
|------|-------|
| 2026-02-07 | Initial AF batch (reduced_dbs) |
| 2026-02-09 | Full AF re-run (full_dbs + fallback) |
| 2026-02-11 | AF complete, 7 AMBER failures identified |
| 2026-02-15 | Boltz batch 1 (233/257) |
| 2026-02-20 | AMBER root cause: non-standard residues |
| 2026-02-22 | Rosetta started (Job 9195328) |
| 2026-03-02 | Boltz complete (257/257) |
| 2026-03-03 | FASTA dedup (135 targets) + consistency fix |
| 2026-03-04 | Crystal-derived FASTAs (all 257) |
| 2026-03-04 | Crystal stripping (36 PDBs) |
| 2026-03-05 | AF + Boltz re-predictions with crystal FASTAs |
| 2026-03-08 | All predictions complete (AF, Boltz, built-in AMBER) |
| 2026-03-14 | Red analysis Phases 1-3, 5-6 complete |
| 2026-03-24 | MolProbity 257/257 both pipelines; PAPER_FINDINGS.md committed |
| 2026-03-29 | Blue 100K / Green 99K .pdb.gz (~50% Rosetta) |
| 2026-04-07 | 1ACB/1ATN AMBER-crystal fill attempts fail (convergence error) |
| 2026-04-16 | Amber-crystal Rosetta fill jobs submitted (88 Blue, 85 Green) |
| 2026-04-16 | Legacy collision dirs purged (340 dirs, 7,775 files) |
| 2026-04-16 | Crystal + Boltz-AMBER partial fills submitted (Blue only; Green clean) |
