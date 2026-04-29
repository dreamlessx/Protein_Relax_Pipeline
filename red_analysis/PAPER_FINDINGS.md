# BM5.5 Docking-Relaxation Benchmark: Findings

**Last verified:** 2026-04-28. Snapshot 2026-04-27a locked. Every fact table at 100% coverage. qc_status = pass.

The DB at `db/bm55.sqlite` is the canonical artifact: 416,340 `rosetta_metrics` rows + 13,364 `prerosetta_metrics` rows + 104,765 `tm_scores` rows + 416,340 `rosetta_energy` rows + 257 `targets` with full metadata + 0 `qc_quarantine` rows. All numbers below are computed from this DB unless noted.

---

## 1. AMBER fixes local geometry without touching global fold

AMBER relaxation simultaneously
- has negligible effect on global fold (TM-score Cliff's d ≈ -0.01)
- has near-saturating effect on local geometry (clashscore Cliff's d ≈ -0.99)

This is not a tradeoff. AMBER fixes what TM-score does not see.

### Paired clashscore deltas (matched by target, n = 257 per pipeline)

| Comparison | Δ Clashscore | Cliff's d | n | p (paired t) |
|---|---|---|---|---|
| AF unrelaxed → AMBER(AF), Blue | -21.09 | -0.992 | 257 | < 1e-17 |
| Boltz → AMBER(Boltz), Blue | -13.49 | -0.989 | 257 | < 1e-17 |
| AF unrelaxed → AMBER(AF), Green | -21.81 | (large) | 257 | < 1e-17 |
| Boltz → AMBER(Boltz), Green | -13.42 | (large) | 257 | < 1e-17 |

### Universality

- AMBER improves MolProbity for 257 of 257 AlphaFold targets (100%).
- AMBER improves MolProbity for 256 of 257 Boltz targets (99.6%).

### Per-source pre-Rosetta MolProbity (Blue, locked, n = 257 per cell)

| Source | Clashscore | MP score | Rama favored % | Rota outliers % |
|---|---|---|---|---|
| `crystal` | 13.85 | 1.81 | 94.94 | 6.15 |
| `amber_crystal` | (clean) | (clean) | (high) | (low) |
| `af_relaxed` | 2.82 | 0.70 | 96.99 | 0.89 |
| `af_unrelaxed` | 23.91 | 1.55 | 97.09 | 1.18 |
| `amber_af` | 2.82 | 0.70 | 97.00 | 0.90 |
| `boltz` | 15.09 | 1.16 | 98.38 | 0.23 |
| `amber_boltz` | **1.60** | **0.45** | **97.76** | 0.36 |

`amber_boltz` produces the cleanest aggregate pre-Rosetta MP (0.45). Defensible recommendation: AMBER is a near-complete clash-removal step at near-zero accuracy cost.

Figures: `fig4_amber_effect`, `fig13_amber_clashscore_effect`, `fig14_amber_dual_effect`.

---

## 2. Crystal structures carry the worst pre-Rosetta MolProbity

| Source | Clashscore | Rota outliers % | MP score |
|---|---|---|---|
| `crystal` | 13.85 | 6.15 | 1.81 |
| `af_relaxed` | 2.82 | 0.89 | 0.70 |
| `amber_boltz` | **1.60** | **0.36** | **0.45** |

Crystal structures show roughly 5× worse clashscore and 7× more rotamer outliers than predictions.

The framing is idealization, not failure. Predicted structures are produced by neural networks that learned ideal geometry. AMBER perfects the physics further. Crystal structures are observed at cryogenic temperatures, in packed lattices, with constraints that make local geometry deviate from gas-phase ideality. The Discussion frames this as the predictions being closer to ideal local geometry, while crystal remains the canonical reference for global fold.

Figures: `fig15_crystal_vs_predicted`.

---

## 3. dualspace_beta wins integrated MolProbity at small TM cost

### Per-target post-Rosetta MP score (Blue, locked, n = 257 per cell)

| Source | cartesian_beta | cartesian_ref15 | dualspace_beta | dualspace_ref15 | normal_beta | normal_ref15 |
|---|---|---|---|---|---|---|
| `af_relaxed` | 0.242 | 0.388 | 0.226 | 0.369 | 0.236 | 0.333 |
| `af_unrelaxed` | 0.243 | 0.386 | 0.227 | 0.353 | 0.287 | 0.386 |
| `amber_boltz` | 0.239 | 0.383 | **0.222** | 0.362 | 0.219 | 0.300 |
| `boltz` | (locked) | (locked) | (locked) | (locked) | 0.213 | (locked) |

Best AMBER + Rosetta combination: `amber_boltz` + `dualspace_beta`, MP = 0.222.
Best standalone Rosetta (no pre-AMBER): `boltz` + `normal_beta`, MP = 0.213.

### Protocol TM-score retention (Blue, `af_relaxed` source, n = 204 with TM data)

| Protocol | Mean TM-score |
|---|---|
| `cartesian_beta` | 0.927 |
| `normal_beta` | 0.926 |
| `normal_ref15` | 0.925 |
| `cartesian_ref15` | 0.924 |
| `dualspace_beta` | 0.919 |
| `dualspace_ref15` | 0.917 |

`dualspace_beta` costs 0.008 to 0.010 TM relative to the top-retentive `cartesian_beta`, and 0.019 to 0.025 relative to the most-retentive cells. The win on MP score outweighs the TM cost.

### Score-function comparison (beta_nov16 vs ref2015)

`beta_nov16` beats `ref2015` on:
- MP score: 42 of 42 (pipeline, source, move-set) triples
- Clashscore: 42 of 42
- Rama favored: 40 of 42

`ref2015` beats `beta_nov16` on:
- Rotamer outliers: 21 of 42 (one of the few `ref2015` wins; consider when sidechain-packing fidelity dominates the use case)

Default recommendation: `beta_nov16`. Use `ref2015` only when rotamer fidelity is the binding constraint.

### Paired AMBER vs no-AMBER tests at the locked frame

Main paired-t table (n_paired = 257 across all 16 cells): zero rows reach raw p < 0.05 at the locked frame. Source: `red_analysis/tables/amber_paired_table2.tsv`.

Rotamer extension (18 rows): 4 rows reach raw p < 0.05, all on AMBER(Boltz) inputs. Smallest BH q = 0.084 (does not survive FDR at q = 0.05). Source: `red_analysis/tables/amber_paired_rotamer.tsv`. Direction: AMBER reduces sidechain packing errors that Rosetta alone (without AMBER pre-conditioning) does not fully correct from Boltz's prediction geometry.

Figures: `fig17_rosetta_protocol_mp`, `figS2_molprobity_violins`, `figS3_rosetta_by_source`, `fig18_amber_vs_rosetta`, `fig20_tm_mp_tradeoff`.

---

## 4. Reproducibility (Blue vs Green)

| Metric | Pearson r | n |
|---|---|---|
| Pre-Rosetta TM | 0.997 | 1,128 |
| Pre-Rosetta RMSD | 0.994 | 1,128 |
| Post-Rosetta TM | 0.999 | 60 |
| Per-source clashscore | 0.867 to 0.991 | 257 |
| Per-source MP score | 0.941 to 0.984 | 257 |

Green's matched-parameter run reproduces Blue's three findings statistically and on a per-target basis.

Figures: `fig11_blue_green_rosetta`.

---

## 5. AlphaFold 2.3.2 vs Boltz-1 v0.4.1 (Blue, n = 257)

| Metric | AF (relaxed) | Boltz | Winner |
|---|---|---|---|
| Clashscore | 2.82 | 15.09 | AF (5× better, AF includes built-in AMBER) |
| Rama favored % | 96.99 | 98.38 | Boltz (better backbone) |
| Rota outliers % | 0.89 | 0.23 | Boltz (4× fewer) |
| MP score | 0.70 | 1.16 | AF |

After standalone AMBER, Boltz catches up: AMBER(Boltz) MP = 0.45 vs AF-relaxed MP = 0.70.

Outliers: 3 catastrophic failures (TM < 0.5) where both methods fail (1WEJ, 2MTA, 1DQJ). 6 targets where AF >> Boltz (ΔTM > 0.1). 5 targets where Boltz >> AF (ΔTM > 0.1). 3H11 shows Boltz stochastic failure (2-3 of 5 models fail).

---

## 6. Data scale (locked snapshot 2026-04-27a)

| Dataset | Rows | Targets | Pipelines | Sources | Protocols |
|---|---|---|---|---|---|
| Pre-Rosetta TM | 12,065 | 257 | 2 | 5 | n/a |
| Pre-Rosetta MolProbity | 13,364 | 257 | 2 | 7 | n/a |
| Post-Rosetta TM | 92,700 | 257 | 2 | 7 | 6 |
| Post-Rosetta MolProbity (`rosetta_metrics`) | **416,340** | 257 | 2 | 7 | 6 |
| Post-Rosetta total energy (`rosetta_energy`) | **416,340** | 257 | 2 | 7 | 6 |

Rosetta arithmetic per cell: 257 targets × 2 pipelines × 27 input structures × 6 protocols × 5 replicates = 416,340. Per pipeline: 27 × 6 × 5 × 257 = 208,170.

0 gap cells. 0 missing rows. 0 NaN. 663 exact-duplicate + 27 legacy-source rows filtered at ingest. Stage C resolved 20 missing Blue crystal pre-Rosetta MP rows from Green source (PDBs verified byte-identical, MolProbity deterministic on identical input).

---

## 7. Figures and tables

### Main figures (20)

1-6: Pre-Rosetta analysis (TM violin, RMSD box, Blue-Green scatter, AMBER effect, AF vs Boltz, outliers).
7-11: Rosetta analysis (protocol comparison, pre/post scatter, protocol × source heatmap, replicate convergence, Blue-Green Rosetta).
12-15: MolProbity analysis (box plots by source, AMBER paired, dual effect scatter, crystal vs predicted).
16-20: Rosetta MolProbity (pre/post clash scatter, protocol MP comparison, AMBER vs Rosetta, MP heatmap, TM-MP tradeoff).

### Supplementary figures (5)

S1: All 7 AMBER MolProbity metrics.
S2: MolProbity violin plots.
S3: Rosetta TM by input type.
S4: MolProbity metric correlation matrix.
S5: AF vs Boltz MolProbity scatter.

### Headline figures for the paper

- `fig14_amber_dual_effect`: two-panel scatter showing TM (no change) vs MolProbity (improvement). Captures Finding 1.
- `fig15_crystal_vs_predicted`: crystal MP burden vs predictions and AMBER. Captures Finding 2.
- `fig17_rosetta_protocol_mp`: Rosetta protocol ranking on MolProbity. Captures Finding 3.
- `fig18_amber_vs_rosetta`: bar chart comparing AMBER vs Rosetta protocols across MP metrics. Headline comparison plot.
- `fig20_tm_mp_tradeoff`: ΔTM vs ΔClashscore scatter. Frames the protocol tradeoff at a glance.

### Tables

| Table | File | Rows |
|---|---|---|
| Per-cell summary | `tables/v_cell_summary` view | 84 |
| Per-target MP means | `tables/per_target_mp_means.tsv` | 21,588 |
| AMBER paired (main) | `tables/amber_paired_table2.tsv` | 16 |
| AMBER paired (rotamer extension) | `tables/amber_paired_rotamer.tsv` | 18 |
| Rosetta protocol summary | `tables/rosetta_source_protocol_summary.tsv` | per-(pipeline, source, protocol) |
| Pairwise Wilcoxon | `tables/pairwise_wilcoxon.tsv` | per-comparison |

LaTeX tables for the manuscript live alongside as `table{1..6}_*.tex`.

---

## 8. What's not in the DB

The DB is complete for the 5 fact tables defined in the schema. Two scope choices worth noting in Methods:

- **AlphaFold 3 not used.** AF2 was chosen for benchmark continuity with prior docking-evaluation literature. Mention as a one-sentence justification in the Discussion.
- **Boltz-1 single-sequence mode.** Intentional, to remove the MSA-quality confound across 257 targets with heterogeneous sequence-database coverage. Methods note required.
- **PoseBusters not run.** Out of paper scope; deprecated for this release. Methodology established in the [`Protein_Data_Analysis`](https://github.com/dreamlessx/Protein_Data_Analysis) Phase 1 pilot (20 proteins, 6,820 structures), retained for the historical record.

The 43 AlphaFold targets that fell back to `reduced_dbs` (HHblits hit the 32,763-residue limit on antibody/immunoglobulin sequences) are not separately tracked in the current DB schema. Source-list trace is pending; mention as a non-blocking caveat in the Discussion.
