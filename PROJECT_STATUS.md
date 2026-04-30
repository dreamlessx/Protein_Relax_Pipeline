# Project Status

**Last verified:** 2026-04-29. Snapshot 2026-04-27a locked. Every fact table at 100% coverage post data-gap fixes (157 green af_unrelaxed TM rows backfilled; 257 Blue amber_crystal DockQ rows added; 1Y64 Boltz outlier identified). qc_quarantine clean.

## Locked DB (single citation point)

| Table | Rows | Coverage |
|---|---|---|
| `rosetta_metrics` | 416,340 | 100.000% (locked) |
| `prerosetta_metrics` | 13,364 | 13,344 base + 20 Stage C Blue crystal backfill |
| `tm_scores` | 105,550 | 12,850 pre-Rosetta (post green af_unrelaxed backfill) + 92,700 post-Rosetta |
| `rosetta_energy` | 416,340 | 100.000% (matches `rosetta_metrics` 1:1 on every join key) |
| `dockq_metrics` | 13,878 | 100.000% across 7 sources × 257 targets × 2 pipelines (post Blue amber_crystal backfill) |
| `targets` | 257 | difficulty + category + n_chains + n_residues + non_standard_flag + parent_pdb_id (4 non-standard) |
| `qc_quarantine` | 0 | clean |
| dimension tables | 2 / 7 / 6 | pipelines / sources / protocols |
| views | 84 / 21,588 | `v_cell_summary` / `v_per_target_mp_means` |

0 orphan rows on every FK relationship. 663 exact-duplicate + 27 legacy-source rows filtered at `rosetta_metrics` ingest. qc_status = pass.

## Dataset coverage

| Quantity | Value |
|---|---|
| BM5.5 targets | 257 |
| Rigid-body / Medium / Difficult | 162 / 60 / 35 |
| Non-standard zlab IDs | 4 (BAAD, BOYV, BP57, CP57; parent_pdb_id populated) |
| Total chains across 257 targets | 605 |
| Total residues across 257 targets | 122,966 |

## Pipeline matrix

| Per target | Per pipeline | Both pipelines |
|---|---|---|
| 27 input structures | 27 × 257 = 6,939 base | 13,878 |
| 810 Rosetta runs (27 × 6 × 5) | 810 × 257 = 208,170 | 416,340 |

| Source bucket | Models | Per-pipeline `rosetta_metrics` rows |
|---|---|---|
| `crystal` | 1 | 7,710 |
| `amber_crystal` | 1 | 7,710 |
| `af_relaxed` | 5 | 38,550 |
| `af_unrelaxed` | 5 | 38,550 |
| `amber_af` | 5 | 38,550 |
| `boltz` | 5 | 38,550 |
| `amber_boltz` | 5 | 38,550 |

Both pipelines (Blue + Green) populated identically across all 7 sources × 6 protocols × 5 reps × 257 targets.

## Pipelines

- **Blue** (primary): scripts under `scripts/`, prefix `blue_`, ACCRE root `/data/p_csb_meiler/agarwm5/protein_pipeline/`
- **Green** (independent verification): in companion repo [`Protein_Ideal`](https://github.com/dreamlessx/Protein_Ideal), prefix `green_`, ACCRE root `/data/p_csb_meiler/agarwm5/protein_ideal_test/`
- **Red analysis** (canonical metrics + figures): under `red_analysis/`, prefix `red_`

Blue and Green use identical Rosetta flags, the same 257 targets, the same 27 input structures per target, the same 6 protocols, and the same 5 replicates. The DB unifies both under snapshot 2026-04-27a.

## Five findings (`red_analysis/PAPER_FINDINGS.md` for full numbers + figure pointers)

1. **AMBER fixes prediction defects without meaningfully changing global fold or interface.** Clashscore Cliff's d = -0.99, mean drop 13-21 clashes per target. TM mean Δ < 0.001 absolute. DockQ mean Δ < 0.0002 (well below DockQ's practitioner threshold of ~0.05). 257/257 AlphaFold + 256/257 Boltz targets improve.
2. **AMBER on crystal references produces a small structural perturbation comparable to crystallographic refinement uncertainty.** Median Δ i-RMSD 0.246 Å (vs Luzzati / Sigma-A 0.2-0.3 Å typical for the BM5.5 ~2.2 Å median resolution); Δ DockQ -0.022 sits below the practitioner decision threshold. The Cliff's d = -1.0 is mathematically tautological (crystal vs crystal = 1.0 by construction). Skip AMBER on crystal sources where the crystal is held as the unmodified reference.
3. **Post-Rosetta paired AMBER effect: small, directionally consistent, below practitioner noise floor.** 0/16 cells reach paired-t p<0.05; 2/16 reach Wilcoxon p<0.05 (consistent direction toward AMBER, magnitude ~1-2% of cell mean). Mean Δ MP sits below the n=257 paired-t MDE (~0.008). AMBER as geometry pre-conditioner (massive Section 1 effect), not as Rosetta-MolProbity scoring enhancer.
4. **`dualspace_beta` and `normal_beta` tie at the top of integrated MolProbity (within 0.01 MP units).** `boltz + normal_beta` MP 0.213 (no pre-AMBER); `amber_boltz + dualspace_beta` MP 0.222 (with pre-AMBER). beta_nov16 dominates ref2015 on 40-42 of 42 (pipeline, source, move-set) triples for MP, clashscore, Rama-favored; ref2015 wins on rotamer outliers (21 of 42 triples). cartesian_beta is the TM-retentive runner-up.
5. **Crystal MolProbity calibration.** Crystal clashscore 13.85 vs predictions 0.68-2.82. Frame as score-function idealization (predictions trained on energy-minimized targets; crystals at cryogenic packing); not a critique of crystallography. 1Y64 named as the single Boltz target where AMBER pre-conditioning makes MolProbity worse (Δ MP +0.46).

## Reproducibility (Blue vs Green agreement)

| Metric | Pearson r | n |
|---|---|---|
| Pre-Rosetta TM | 0.997 | 1,128 |
| Pre-Rosetta RMSD | 0.994 | 1,128 |
| Post-Rosetta TM | 0.999 | 60 |
| Per-source clashscore | 0.867 to 0.991 | 257 |
| Per-source MP score | 0.941 to 0.984 | 257 |

The Green run statistically reproduces Blue.

## Resolved issues at lock

| Issue | Resolution |
|---|---|
| 1ACB / 1ATN AMBER-crystal divergence | `amber_relax_crystal_v5.py` detects peptide bonds > 2.5 Å, splits affected chains. 1ACB E energy went from 8,350 to -5,986 kcal/mol; 1ATN E from 20,790 to -12,079. |
| 20 missing Blue crystal pre-Rosetta MP rows | Stage C backfill from Green crystal MP. Verified byte-identical Blue/Green crystal MP on 1A2K, 7CEI, 1ACB sample; MolProbity is deterministic on identical input PDBs. |
| `1E96` target ID coerced to `1e+96` in two combined TSVs | Pandas-without-dtype-str scientific-notation coercion in upstream aggregation. Loader normalizes via `TARGET_ALIASES` until upstream is fixed. Per-target TSVs preserve `1E96` correctly. |
| `rosetta_energy` extraction at 44% coverage | `extract_rosetta_energy.py` patched: added `amber_crystal` `classify_source` branch (recovered 15,420 cells), added per-rep `score_*.sc` sidecar walk (recovered ~217k cells), added `POSE_ENERGIES_TABLE` PDB fallback (recovered 1 final cell). Now 100%. |
| Legacy collision dirs (`amber_af_relaxed`, `amber_boltz_relaxed`) | 340 dirs and 7,775 stale `.pdb.gz` files removed 2026-04-16. Authoritative per-model data is in `amber_af_ranked_{0..4}/` and `amber_boltz_model_{0..4}/`. Legacy src_types filtered at every ingest path. |

## Validation status

| Validation | Blue | Green |
|---|---|---|
| MolProbity (pre-Rosetta) | 257/257 | 257/257 |
| MolProbity (post-Rosetta) | 208,170 rows | 208,170 rows |
| TM-score (pre-Rosetta) | 257/257 | 257/257 |
| TM-score (post-Rosetta) | full | full |
| Rosetta total energy | 208,170 cells (100%) | 208,170 cells (100%) |

## Schema additions vs original 2026-04-20 schema

- `rosetta_energy` table: total_score + per_residue_energy per (target, pipeline, src_type, protocol, rep) cell.
- `targets.parent_pdb_id` column: maps non-standard zlab IDs to canonical parent PDB + chain selection.

Both applied idempotently by `db/scripts/build_db_supplements.py` (CREATE TABLE IF NOT EXISTS + ALTER TABLE ADD COLUMN with try/except on duplicate).

## Release artifact

`db-2026-04-27a-supp` on this repo. 6 assets: `bm55.sqlite` (~225 MB, VACUUMed) + 5 raw TSVs (`combined_molprobity.tsv`, `combined_tmscore.tsv`, `combined_rosetta_molprobity.tsv`, `combined_rosetta_tmscore.tsv`, `combined_rosetta_energy.tsv`).

## Timeline

| Date | Event |
|---|---|
| 2026-02-07 | Initial AlphaFold batch (reduced_dbs) |
| 2026-02-09 | Full AlphaFold re-run with full_dbs + reduced_dbs fallback |
| 2026-02-15 | Boltz-1 batch (single-sequence, 233/257 first pass) |
| 2026-02-22 | Rosetta started |
| 2026-03-02 | Boltz complete (257/257) |
| 2026-03-04 | Crystal-derived FASTAs replace UniProt; 36 PDBs stripped of homo-multimer duplicates |
| 2026-03-08 | All predictions complete (AF + Boltz + built-in AMBER + standalone AMBER) |
| 2026-03-24 | MolProbity 257/257 both pipelines |
| 2026-04-16 | 1ACB/1ATN AMBER-crystal v5 fix; legacy collision dirs purged |
| 2026-04-20 | Data integrity audit; legacy contamination filtered at ingest |
| 2026-04-27 | Snapshot 2026-04-27a locked at 100.000% (416,340 / 416,340 `rosetta_metrics` rows) |
| 2026-04-28 | Stage B + C: `prerosetta_metrics` + `tm_scores` + `targets` metadata loaded under same snapshot |
| 2026-04-28 | Stage B6 + B7 + B8: `rosetta_energy` recovered to 100%, `targets.parent_pdb_id` populated, `qc_quarantine` cleaned |
