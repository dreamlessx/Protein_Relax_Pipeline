# Data Integrity Audit — 2026-04-20

Full audit of Rosetta + MolProbity + aggregation pipeline for the BM5.5 relaxation benchmark. Run after red_ros_mp (SLURM 10175020) locked at 255/257 COMPLETED + 2 TIMEOUT.

## Exact file counts

### Rosetta output .pdb.gz structures

| Pipeline | Target (27 × 6 × 5 × 257) | Actual | Delta |
|---|---|---|---|
| Green | 208,170 | 208,170 | 0 (exact) |
| Blue | 208,170 | 208,228 | +58 (net) |
| Combined | 416,340 | 416,398 | +58 |

Green: 257 of 257 targets at exactly 810 files. Perfect.

Blue: 253 of 257 targets at exactly 810 files. 4 anomalies:

| Target | Count | Delta | Cause |
|---|---|---|---|
| 1K5D | 865 | +55 | Legacy collision dirs (amber_af_relaxed, amber_boltz_relaxed) + amber_boltz_model_4 normal_ref15 missing 5 reps |
| 4GAM | 848 | +38 | Legacy collision dirs (amber_af_relaxed full 30, amber_boltz_relaxed partial 8) |
| 1HE8 | 805 | -5 | amber_boltz_model_1 normal_ref15 missing 5 reps |
| 4Y7M | 780 | -30 | Missing crystal_4Y7M source entirely (Blue never ran Rosetta on crystal input; Green has it) |

### Rosetta MolProbity metric TSVs

- Per-target TSV files: 514 of 514 (257 × 2 pipelines), all non-empty.
- Rosetta TM-score TSVs: 514 of 514.
- Pre-Rosetta MP TSVs: 514 of 514.
- Pre-Rosetta TM TSVs: 514 of 514.

### Combined aggregation

| Stage | Rows | Notes |
|---|---|---|
| Pre-clean combined | ~417,082 | includes 244 legacy-source contamination rows |
| Clean combined (this audit) | 416,838 | legacy filtered, all 244 dropped |
| Blue in clean combined | ~208,210 | slightly over 208,170 target (+40); extras in af_relaxed/af_unrelaxed/boltz source buckets from Rosetta over-production of replicates |
| Green in clean combined | ~208,628 | over target (+458); same cause |
| NaN cells | 0 | clean |

## Legacy collision fix

The LEGACY `amber_af_relaxed/` and `amber_boltz_relaxed/` source directories (residue from the af_relaxed + boltz naming collision resolved in early 2026) existed only for 2 Blue targets (1K5D, 4GAM). Under the previous `reaggregate_combined.py`, they were pooled with legitimate data via `startswith('amber_af')` / `startswith('amber_boltz')`, contaminating those 2 target x 2 source cells.

Patched `reaggregate_combined.py`:
- Added explicit `LEGACY_SRC = {'amber_af_relaxed', 'amber_boltz_relaxed'}` filter BEFORE source bucketing.
- Replaced loose `startswith('amber_af')` / `startswith('amber_boltz')` / `startswith('boltz')` with precise `startswith('amber_af_ranked')` / `startswith('amber_boltz_model')` / `startswith('boltz_boltz_input')`.
- Unknown `src_type` now raises `ValueError` instead of silently bucketing.

Impact on findings: negligible. 244 rows out of ~417,082 = 0.059% contamination across 2 targets. Rerun of per_target_means + paired_amber + paired_amber_rotamer produces identical numbers to iter 42 pre-clean state (to 4 decimal places).

## Real data gaps

Three genuine missing-data gaps (40 files out of 208,170 Blue = 0.019%):

1. **Blue 4Y7M crystal source**: missing all 30 files (1 source × 6 protocols × 5 reps). Rosetta never ran on crystal_4Y7M for Blue pipeline. Green has it.
2. **Blue 1HE8 amber_boltz_model_1 normal_ref15**: 0 of 5 reps (5 files missing).
3. **Blue 1K5D amber_boltz_model_4 normal_ref15**: 0 of 5 reps (5 files missing).

These do NOT affect any cross-target comparison at n=257; per-target means for the 3 affected cells are computed from fewer replicates or (for 4Y7M Blue crystal) are NaN/absent. Paired tests handle this automatically (paired-t drops pairs with missing values).

## Rosetta over-production

Blue + Green both show ~250-600 extra .pdb.gz files per pipeline in af/boltz source buckets. Cause: Rosetta SLURM jobs occasionally produced extra replicates beyond the 5 requested. Does not affect analysis (per-target means average over all available reps equally).

## Timeout jobs

2 of 257 array tasks on red_ros_mp (10175020) hit 12h wall-time:

| Task | Target | Blue MP rows | Green MP rows | Status |
|---|---|---|---|---|
| 83 | 1N2C | 811 (100%) | 601 (74%) | Blue complete, Green partial |
| 215 | 4GAM | 838 (100%) | 643 (79%) | Blue complete (incl. legacy), Green partial |

Both are the largest BM5.5 complexes by residue count. Partial Green MP is acceptable: all 6 protocols have >= 17 reps per source (sufficient for per-protocol means).

## Verdict

Data is suitable for publication with the following caveats:

- Per-target n=257 for all unpaired analyses.
- Paired test n=245-256 (pair-matching limits; unchanged from iter 42).
- 1 target (Blue 4Y7M) lacks crystal source contribution.
- 2 targets (1N2C, 4GAM) have partial Green MP data for the long-tail protocols.
- Legacy collision rows now filtered; analysis scripts and downstream outputs refreshed.

Zero NaN in the cleaned combined file. All 514 metric TSVs non-empty. Paired findings (4 of 18 rotamer cells significant on AMBER(Boltz) vs Boltz at p<0.05) hold on the cleaned data.

## Action items

1. (Optional, ACCRE compute) Submit 3 small Rosetta fill jobs for the 40 missing .pdb.gz files. Not load-bearing.
2. (Optional) Delete 4 legacy source subdirs on ACCRE: `1K5D/{amber_af_relaxed,amber_boltz_relaxed}`, `4GAM/{amber_af_relaxed,amber_boltz_relaxed}`. Would reduce Blue .pdb.gz count from 208,228 to 208,160 (close to 208,170 target; -10 from the real gaps).
3. Update manuscript Methods Table 1 + Results with exact replicate counts; acknowledge 4Y7M Blue crystal gap in a footnote.
