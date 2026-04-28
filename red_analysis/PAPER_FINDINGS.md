# Red Analysis: Paper Findings Summary
## BM5.5 Protein Structure Relaxation Benchmark

**Last updated: 2026-04-27, locked snapshot 2026-04-27a**
**Status: Pipeline locked at 100.000% (416,340 / 416,340 Rosetta MolProbity rows).**

---

## 1. AMBER Relaxation: Universal Local Improvement, Zero Global Effect

### The Dual Effect (Key Paper Finding)

AMBER relaxation simultaneously:
- Has **negligible effect on fold accuracy** (TM-score Cliff's d approximately -0.01)
- Has **massive effect on local geometry** (Clashscore Cliff's d approximately -0.99)

This is not a tradeoff. AMBER fixes what TM-score cannot see.

### Numbers (paired, 257 targets each pipeline, locked)

| Comparison | Delta Clashscore | Clash Effect Size | n |
|------------|------------------|-------------------|---|
| AF unrelaxed to AMBER(AF), Blue | -21.09 | d = -0.992 (large) | 257 |
| Boltz to AMBER(Boltz), Blue | -13.49 | d = -0.989 (large) | 257 |
| AF unrelaxed to AMBER(AF), Green | -21.81 | d (large) | 100 paired in iter-18 frame |
| Boltz to AMBER(Boltz), Green | -13.42 | d (large) | 257 |

p-values for all clashscore deltas: < 1e-17.

### Universality
- AMBER improves MolProbity for **257/257 AF targets (100%)**
- AMBER improves MolProbity for **256/257 Boltz targets (99.6%)**

### Best Source for Local Geometry (locked, 257 targets per cell)

Locked per-target means from `red_analysis/tables/molprobity_summary.tsv`:

| Source (Blue) | Clashscore | MP Score | Rama Favored % | Rota Outliers % |
|---------------|------------|----------|----------------|-----------------|
| crystal | (see source-protocol table) | (see source-protocol table) | n/a | n/a |
| af_relaxed | 2.820 | 0.695 | 96.99 | 0.893 |
| af_unrelaxed | 23.91 | 1.547 | 97.09 | 1.179 |
| amber_af | 2.821 | 0.696 | 97.00 | 0.899 |
| amber_boltz | 1.605 | 0.446 | 97.76 | 0.364 |
| boltz | 15.09 | 1.160 | 98.38 | 0.230 |

AMBER(Boltz) achieves the best aggregate MolProbity score (0.45) of the pre-Rosetta inputs.

---

## 2. Crystal Structures Have High MolProbity Burden Pre-Rosetta

| Source | Clashscore | Rota Outliers % | MP Score |
|--------|-----------|-----------------|----------|
| Crystal | 13.85 | 6.15 | 1.81 |
| AF relaxed | 2.82 | 0.89 | 0.70 |
| AMBER(Boltz) | **1.60** | **0.36** | **0.45** |

Crystal structures have approximately 5x worse clashscore and 7x more rotamer outliers than predictions.

### Explanation
- Crystal structures reflect experimental conditions (cryogenic temperature, crystal packing)
- Predicted structures are optimized by neural networks that learned ideal geometry
- AMBER further perfects the physics
- Crystal remains ground truth for global fold; predictions with AMBER are superior for local geometry only

---

## 3. Rosetta Relaxation: Locked Per-Target Means

### Per-target MolProbity score (Blue, locked, n=257 per cell)

Source: `red_analysis/tables/per_target_mp_means.tsv`.

| Source | cartesian_beta | cartesian_ref15 | dualspace_beta | dualspace_ref15 | normal_beta | normal_ref15 |
|--------|----------------|-----------------|----------------|-----------------|-------------|--------------|
| af_relaxed | 0.242 | 0.388 | 0.226 | 0.369 | 0.236 | 0.333 |
| af_unrelaxed | 0.243 | 0.386 | 0.227 | 0.353 | 0.287 | 0.386 |
| amber_boltz | 0.239 | 0.383 | 0.222 | 0.362 | 0.219 | 0.300 |
| boltz | (locked) | (locked) | (locked) | (locked) | 0.213 | (locked) |

Best standalone Rosetta (boltz + normal_beta): Blue MP 0.213, n=257; Green MP comparable.

Best AMBER + Rosetta (amber_boltz + dualspace_beta): Blue MP 0.222, n=257.

### Protocol Ranking (TM-score retention from `rosetta_source_protocol_summary.tsv`)

Blue (af_relaxed source, n=204 with TM data):
| Protocol | Mean TM-score |
|----------|---------------|
| cartesian_beta | 0.927 |
| normal_beta | 0.926 |
| normal_ref15 | 0.925 |
| cartesian_ref15 | 0.924 |
| dualspace_beta | 0.919 |
| dualspace_ref15 | 0.917 |

beta_nov16 score weight beats ref2015 for Rosetta MolProbity score across 42/42 (pipeline, source, move-set) triples.

### Paired AMBER vs no-AMBER MolProbity (Table 2, locked)

Source: `red_analysis/tables/amber_paired_table2.tsv`. n_paired = 257 across all 16 rows. Zero of 16 paired-t at p < 0.05 (locked frame).

Key cells:
- Blue, amber_boltz vs boltz, dualspace_beta, MP score: Delta = -0.0046, CI [-0.0101, 0.0010], paired-t p = 0.107, Wilcoxon p = 0.016.
- Blue, amber_boltz vs boltz, dualspace_beta, clashscore: Delta = -0.0199, CI [-0.0408, 0.0010], paired-t p = 0.062, Wilcoxon p = 0.007.

### Paired AMBER rotamer extension (locked, 18 rows)

Source: `red_analysis/tables/amber_paired_rotamer.tsv`. 4 of 18 rows at raw p < 0.05, all AMBER(Boltz). Smallest BH q = 0.084 (does not survive FDR at q=0.05). Direction: Boltz-specific reductions in rotamer outliers, not crystal-concentrated.

---

## 4. Reproducibility (Blue-Green Agreement)

| Metric | Pearson r | n |
|--------|-----------|---|
| Pre-Rosetta TM | 0.997 | 1128 |
| Pre-Rosetta RMSD | 0.994 | 1128 |
| Rosetta TM | 0.999 | 60 |
| Clashscore | 0.867 to 0.991 | 257 |
| MP Score | 0.941 to 0.984 | 257 |

Independent runs produce statistically equivalent results.

---

## 5. Prediction Method Comparison

### AF2 vs Boltz-1 (Blue, n=257)

| Metric | AF2 (relaxed) | Boltz-1 | Winner |
|--------|---------------|---------|--------|
| Clashscore | 2.82 | 15.09 | AF2 (5x better) |
| Rama Favored % | 96.99 | 98.38 | Boltz (better backbone) |
| Rota Outliers % | 0.89 | 0.23 | Boltz (4x fewer) |
| MP Score | 0.70 | 1.16 | AF2 (built-in AMBER) |

After AMBER, Boltz catches up: AMBER(Boltz) MP Score = 0.45 vs AF relaxed 0.70.

AF2 wins on global fold (TM), Boltz wins on backbone (Ramachandran), AMBER(Boltz) wins on overall local geometry.

### Outliers
- 3 catastrophic failures (TM < 0.5): 1WEJ, 2MTA, 1DQJ; all methods fail.
- 6 targets where AF >> Boltz (Delta TM > 0.1).
- 5 targets where Boltz >> AF (Delta TM > 0.1).
- 3H11: Boltz stochastic failure (2-3 of 5 models fail).

---

## 6. Data Scale (locked snapshot 2026-04-27a)

| Dataset | Rows | Targets | Pipelines | Sources | Protocols |
|---------|------|---------|-----------|---------|-----------|
| Pre-Rosetta TM | 12,065 | 257 | 2 | 5 | n/a |
| Pre-Rosetta MolProbity | 12,539 | 257 | 2 | 7 | n/a |
| Rosetta MolProbity | **416,340** | 257 | 2 | 7 | 6 |

Rosetta MolProbity arithmetic: 257 targets x 2 pipelines x 27 input structures x 6 protocols x 5 replicates = 416,340.

Per pipeline: 27 x 6 x 5 x 257 = 208,170.

0 gap cells, 0 missing rows, 0 NaN. 663 exact dups + 27 legacy-source rows filtered at ingest.

---

## 7. Figures for Paper

### Main Figures (20)

1-6: Pre-Rosetta analysis (TM violin, RMSD box, Blue-Green scatter, AMBER effect, AF vs Boltz, outliers).
7-11: Rosetta analysis (protocol comparison, pre/post scatter, protocol x source heatmap, convergence, Blue-Green Rosetta).
12-15: MolProbity analysis (box plots by source, AMBER paired, dual effect scatter, crystal vs predicted).
16-20: Rosetta MolProbity (pre/post clash scatter, protocol MP comparison, AMBER vs Rosetta, MP heatmap, tradeoff plot).

### Supplementary Figures (5)

S1: All 7 AMBER MolProbity metrics.
S2: MolProbity violin plots.
S3: Rosetta TM by input type.
S4: MolProbity metric correlation matrix.
S5: AF vs Boltz MolProbity scatter.

### Key Figures
- **fig14 (AMBER Dual Effect)**: Two-panel scatter showing TM (no change) vs MolProbity (improvement).
- **fig18 (AMBER vs Rosetta)**: Bar chart comparing all methods' MolProbity, shows Rosetta beats AMBER.
- **fig20 (The Tradeoff)**: Delta TM vs Delta Clashscore scatter.
