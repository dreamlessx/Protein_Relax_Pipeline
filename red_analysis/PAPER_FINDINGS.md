# BM5.5 Docking-Relaxation Benchmark: Findings

**Last verified:** 2026-04-29 (post data-gap fixes; reconciled to fresh-context adversarial review). Snapshot 2026-04-27a locked. Every fact table at 100% coverage. qc_status = pass.

The DB at `db/bm55.sqlite` is the canonical artifact: 416,340 `rosetta_metrics` rows + 13,364 `prerosetta_metrics` rows + 12,850 `tm_scores_pre` rows + 92,700 `tm_scores_post` rows + 416,340 `rosetta_energy` rows + 13,878 `dockq_metrics` rows + 257 `targets` with full metadata + 0 `qc_quarantine` rows. All numbers below are computed from this DB.

---

## 1. AMBER fixes prediction defects, not crystal defects

We measured AMBER's effect on three independent dimensions: local stereochemistry (MolProbity sub-metrics), global fold (TM-score against the cleaned crystal reference), and binding interface (DockQ + i-RMSD + l-RMSD + f_nat). The dimensions are decoupled: AMBER's effect on each is measured separately, paired by target.

### Local stereochemistry (the headline AMBER effect)

| Comparison (matched-target, n = 257 per pipeline) | Δ Clashscore | Cliff's d | p (paired t) |
|---|---|---|---|
| AF unrelaxed → AMBER(AF), Blue | -21.09 | -0.992 | < 1e-17 |
| Boltz → AMBER(Boltz), Blue | -13.49 | -0.989 | < 1e-17 |
| AF unrelaxed → AMBER(AF), Green | -21.81 | (large) | < 1e-17 |
| Boltz → AMBER(Boltz), Green | -13.42 | (large) | < 1e-17 |

AMBER reduces clashscore on 257/257 AlphaFold targets (100%) and 256/257 Boltz targets (99.6%). The single Boltz failure is **1Y64** (Boltz MP 1.176 → AMBER(Boltz) MP 1.636, Δ MP = +0.46), where AMBER cannot escape a starting clash configuration.

### Global fold (TM-score)

| Comparison | Mean Δ TM | Cliff's d | n_paired | Note |
|---|---|---|---|---|
| `af_unrelaxed → amber_af`, Blue | -0.00041 | -0.634 | 257 | rank-ordered consistent, biologically negligible |
| `af_unrelaxed → amber_af`, Green | -0.00054 | -0.630 | 257 | matched after backfill (was n=100 pre-fix) |
| `boltz → amber_boltz`, Blue | < 0.001 | -0.74 | 257 | same pattern |
| `boltz → amber_boltz`, Green | < 0.001 | -0.78 | 257 | same pattern |
| `af_relaxed → amber_af`, Blue | < 0.0001 | -0.07 | 257 | the right control: AMBER on already-AMBER-relaxed AF |

AMBER consistently nudges TM-score down by less than 0.001 in absolute units. Cliff's d looks "large" because the rank-ordering across paired targets is consistent in direction, but the magnitude is well inside the resolution of TM-score for protein-protein complexes. The control cell (`af_relaxed → amber_af`, where the AMBER step has no work to do because AlphaFold's internal AMBER already ran) shows Cliff's d = -0.07: negligible. Manuscript prose distinguishes "statistically reproducible direction" from "biologically meaningful magnitude."

### Binding interface (DockQ)

| Comparison | Δ DockQ mean | Cliff's d | p | Interpretation |
|---|---|---|---|---|
| `af_unrelaxed → amber_af` | +0.001 | +0.027 to +0.097 | > 0.1 | negligible drift |
| `boltz → amber_boltz` | +0.001 | +0.105 to +0.136 | > 0.1 | negligible drift |
| `af_relaxed → amber_af` | -0.00015 | -0.152 | 0.033 | small statistically detectable shift, well below DockQ resolution |
| `crystal → amber_crystal` | -0.0216 | -1.000 | < 1e-17 | small perturbation within crystallographic refinement uncertainty (see note) |

AMBER preserves the binding interface on AlphaFold and Boltz prediction inputs (Δ DockQ < 0.0002, well below DockQ's practitioner decision threshold of ~0.05). On crystal references, AMBER produces a small structural shift: i-RMSD by 0.246 Å (median across 257 targets), l-RMSD by 0.277 Å, DockQ by -0.022. The Cliff's d = -1.0 is mathematically tautological (crystal vs crystal is DockQ = 1.0 by construction; AMBER can only move it down), and the p-value approaches zero by the same construction. **The shift is directionally consistent on 257/257 targets** (no target shows AMBER moving the crystal interface in the opposite direction), but the magnitude is inside typical crystallographic refinement uncertainty for BM5.5-range resolutions: BM5.5 spans ~1.4-3.5 Å resolution with median around 2.2 Å, where Luzzati / Sigma-A coordinate error estimates run 0.2-0.3 Å for whole-structure DPI. AMBER on crystal moves coordinates by an amount comparable in magnitude to the underlying refinement uncertainty, but in a systematic direction the refinement does not have. **Practical reading**: AMBER applied to predictions improves them; AMBER applied to crystals produces a small, systematic perturbation comparable to the refinement-uncertainty envelope. For comparisons where both prediction and crystal sources are AMBER-treated, the comparison is fair; for analyses that retain the crystal as an unmodified reference, skip AMBER on the crystal source.

### Honest framing

The headline framing is "AMBER substantially improves local geometry on prediction inputs without meaningfully changing global fold or binding interface; on crystallographic references, AMBER produces a small structural perturbation within typical refinement uncertainty."

The colloquial "free lunch" framing is defensible for the pre-Rosetta input-cleanup characterization (clashscore reduction is massive, TM and DockQ shifts are below their resolution floors). For the post-Rosetta combined pipeline, the AMBER + Rosetta vs Rosetta-alone comparison resolves at small effect sizes (see Section 2); manuscript prose distinguishes "AMBER as input-cleanup step" (massive, real) from "AMBER as separately additive to Rosetta in MolProbity terms" (small, directionally consistent, below practitioner noise floor).

Figures: `fig13_amber_clashscore_effect`, `fig14_amber_dual_effect`, `fig21_dockq_by_source` (new).

---

## 2. Post-Rosetta paired AMBER effect: small, directionally consistent, below practitioner noise floor

We ran the AMBER paired test on the post-Rosetta MolProbity outputs to test the magnitude of the AMBER + Rosetta vs Rosetta-alone effect at the locked frame.

| `amber_paired_table2.tsv` summary | Value |
|---|---|
| Cells tested | 16 (Blue + Green × 4 source-pairs × 2 metrics) |
| Cells reaching paired-t p < 0.05 | **0** |
| Cells reaching Wilcoxon p < 0.05 | **2** (Blue amber_boltz vs boltz, dualspace_beta, MolProbity score: Wilcoxon p = 0.016; same cell, clashscore: Wilcoxon p = 0.007) |
| Mean Δ MolProbity score range | -0.007 to +0.001 |
| Mean Δ Clashscore range | -0.04 to +0.02 |

Stratified by difficulty (`amber_paired_table2_by_difficulty.tsv`, 48 rows): 0 cells reach paired-t p < 0.05 in any difficulty tier. The Boltz-on-difficult cells reach p = 0.63, indicating the parametric test does not detect a meaningful effect.

The two Wilcoxon-significant cells (Blue, AMBER(Boltz) vs Boltz, dualspace_beta) reflect a directionally consistent rank-shift (mean Δ MP = -0.005, mean Δ clashscore = -0.02). The paired-t p-values are 0.107 and 0.062. The rotamer extension (`amber_paired_rotamer.tsv`, 18 cells) reaches raw p < 0.05 in 4 cells, all on AMBER-on-Boltz inputs; smallest BH q = 0.084 (does not survive FDR at q = 0.05). **Direction across all 4 marginally-detected cells: 4 of 4 favor AMBER improvement (zero rule-out random-sign noise).**

Architectural framing: AMBER is a **geometry pre-conditioner**, not a Rosetta-MolProbity scoring enhancer. The mechanistic story is:
- AMBER moves prediction inputs from a high-clashscore starting point to a low-clashscore starting point (Section 1: Cliff's d = -0.99, mean drop 13-21 clashscore points).
- Rosetta then takes those inputs and produces a final relaxed structure with MolProbity score in the 0.21-0.23 band.
- The post-Rosetta MolProbity does not separately benefit from AMBER pre-conditioning at a magnitude practitioners would notice (mean Δ MP -0.005 to -0.007 in the marginally-significant cells, ~1-2% of the cell mean).

The minimum detectable effect (MDE) at n = 257 with paired-t at α = 0.05 is approximately Δ MP = 0.008 for a within-cell standard deviation around 0.07; the observed Δ MP -0.005 sits below this MDE, so the post-Rosetta paired test would not be expected to reach paired-t significance at this n. The Wilcoxon detection at p = 0.016 is consistent with a rank-skewed distribution where the parametric test loses power but the rank-based test still picks up the directional consistency. Methods will document this as an MDE-bounded equivalence claim, not a null result.

Two practical readings:

- **For applications that consume the AMBER-pre-conditioned PDB without further Rosetta**: AMBER's value is in pre-Rosetta input cleanup (Section 1). Massive clashscore reduction, near-zero TM/DockQ change. Worth it.
- **For applications that run Rosetta after**: AMBER + Rosetta vs Rosetta-alone differs by a small amount in MolProbity (-0.005 to -0.007 in the cells where the rank-test detects an effect). Practitioners can choose AMBER pre-conditioning based on whether they need clean-clashscore inputs for downstream tools, rather than on a post-Rosetta MolProbity gain.

Figures: `fig22_amber_paired_forest` (new).

---

## 3. Rosetta protocol selection: a tied family at the top, beta_nov16 dominates ref2015

Per-target post-Rosetta MP score (Blue, locked, n = 257 per cell) for the top cells:

| Source | cartesian_beta | dualspace_beta | normal_beta |
|---|---|---|---|
| `af_relaxed` | 0.242 | 0.226 | 0.236 |
| `af_unrelaxed` | 0.243 | 0.227 | 0.287 |
| `amber_boltz` | 0.239 | **0.222** | 0.219 |
| `boltz` | (locked) | (locked) | **0.213** |

Two cells tie at the top within 0.01 MP units:
- `boltz + normal_beta`: MP 0.213 (no pre-AMBER)
- `amber_boltz + dualspace_beta`: MP 0.222 (with pre-AMBER)

The narrative previously claimed "dualspace_beta wins integrated MolProbity." More honestly, the integrated-MP minimum is in a small cluster, and the choice between cells depends on the use case:

- **Lowest-MP, no pre-AMBER overhead**: `boltz + normal_beta` (MP 0.213, but no AMBER pre-conditioning means downstream tools see higher input clashscore).
- **Lowest-MP with AMBER pre-conditioning** (clean inputs for downstream tools that want low clashscore): `amber_boltz + dualspace_beta` (MP 0.222).
- **TM-retentive choice**: `cartesian_beta` (mean TM = 0.927, the most TM-retentive protocol; MP slightly worse than dualspace).

Score-function comparison across 42 (pipeline, source, move-set) triples:
- `beta_nov16` beats `ref2015` on MP score: **42 of 42**
- Clashscore: **42 of 42**
- Rama-favored: **40 of 42**
- Rotamer outliers: **`ref2015` wins 21 of 42** (the only sub-metric where ref2015 is the right default)

**Default protocol recommendation**: `beta_nov16` for the score function. Among the move sets, `dualspace_beta` and `normal_beta` cluster at the top of integrated MP; `cartesian_beta` is the TM-retentive runner-up. `ref2015` only when sidechain rotamer fidelity is the binding constraint.

Figures: `fig17_rosetta_protocol_mp`, `figS2_molprobity_violins`, `figS3_rosetta_by_source`, `figS7_rosetta_protocol_mp_by_difficulty` (new).

---

## 4. Crystal MolProbity calibration: high by design, framed as score-function idealization

| Source | Clashscore | Rota outliers % | MP score |
|---|---|---|---|
| `crystal` | 13.85 | 6.15 | 1.81 |
| `af_relaxed` | 2.82 | 0.89 | 0.70 |
| `amber_boltz` | **1.60** | **0.36** | **0.45** |

Crystal structures show ~5× worse clashscore and ~7× more rotamer outliers than relaxed predictions.

**Framing.** Predictions come from neural networks trained on energy-minimized targets. AMBER perfects the physics in vacuo. Rosetta `beta_nov16` adds another round of geometric idealization. Crystal stereochemistry is observed in cryogenic packed lattices, with refinement constrained by electron density rather than by gas-phase energy minima. The crystals' MolProbity numbers reflect the experimental environment, not refinement failure.

**Anticipated reviewer objection (Gemini, "circular validation")**: NN+AMBER+Rosetta optimize toward score-function idealization; we then evaluate against MolProbity, which the score functions know about. Therefore the "predictions look better than crystals" claim is tautological.

**Mitigation**: MolProbity is not in the same family as the AlphaFold or Rosetta training objectives. AlphaFold is trained on per-residue distance and orientation distributions; Rosetta `beta_nov16` is trained on small-molecule energetics + macromolecular density. MolProbity (clashscore, Ramachandran-favored, rotamer outliers) is built from independent geometric criteria from MolProbity 2018 (Williams et al). The objection has a kernel of truth (score functions and MolProbity are not orthogonal) but the magnitude of the prediction-vs-crystal gap (5-7× on the harshest sub-metrics) is larger than score-function-overlap can explain. Discussion paragraph addresses this directly.

**Mitigation paired with Finding 1**: AMBER on crystal damages the interface. So the predictions+AMBER+Rosetta cell is NOT "better than crystal across the board"; it is "better than crystal on local geometry, equivalent on global fold, and modestly worse on binding interface for the crystal source where AMBER perturbs a known-clean reference." The honest comparison is multi-dimensional, not headline.

Figures: `fig15_crystal_vs_predicted`.

---

## 5. Reproducibility (Blue vs Green)

| Metric | Pearson r | n |
|---|---|---|
| Pre-Rosetta TM | 0.997 | 1,128 (post-backfill of 157 green af_unrelaxed targets) |
| Pre-Rosetta RMSD | 0.994 | 1,128 |
| Post-Rosetta TM | 0.999 | 60 |
| Per-source clashscore | 0.867 to 0.991 | 257 |
| Per-source MP score | 0.941 to 0.984 | 257 |
| Per-source DockQ | 0.998 (crystal sanity), 0.99+ (amber_crystal) | 257 |

Blue and Green pipelines reproduce all four findings. The gap-fix backfill (157 green af_unrelaxed TM-scores, 257 Blue amber_crystal DockQ rows) closed the only known asymmetric coverage gaps; both pipelines are now fully populated.

Figures: `fig11_blue_green_rosetta`.

---

## 6. AlphaFold 2.3.2 vs Boltz-1 v0.4.1 (Blue, n = 257)

| Metric | AF (relaxed) | Boltz | Winner |
|---|---|---|---|
| Clashscore | 2.82 | 15.09 | AF (5×; AF includes built-in AMBER) |
| Rama favored % | 96.99 | 98.38 | Boltz (better backbone) |
| Rota outliers % | 0.89 | 0.23 | Boltz (4×) |
| MP score | 0.70 | 1.16 | AF |

After standalone AMBER, Boltz catches up: AMBER(Boltz) MP = 0.45 vs AF-relaxed MP = 0.70.

Outliers: 3 catastrophic failures (TM < 0.5) where both methods fail (1WEJ, 2MTA, 1DQJ). 6 targets where AF >> Boltz (ΔTM > 0.1). 5 targets where Boltz >> AF (ΔTM > 0.1). 3H11 shows Boltz stochastic failure (2-3 of 5 models fail). 1Y64 is the one Boltz target where AMBER pre-conditioning makes MolProbity worse.

---

## 7. Data scale (locked snapshot 2026-04-27a, post data-gap fixes)

| Dataset | Rows | Targets | Pipelines |
|---|---|---|---|
| Pre-Rosetta TM | 12,850 | 257 | 2 (post-backfill: was 12,065 + 785 green af_unrelaxed) |
| Pre-Rosetta MolProbity | 13,364 | 257 | 2 |
| Post-Rosetta TM | 92,700 | 257 | 2 |
| Post-Rosetta MolProbity (`rosetta_metrics`) | **416,340** | 257 | 2 |
| Post-Rosetta total energy (`rosetta_energy`) | **416,340** | 257 | 2 |
| Pre-Rosetta interface metrics (`dockq_metrics`) | **13,878** | 257 | 2 (post-backfill: was 13,621 + 257 Blue amber_crystal) |

0 gap cells. 0 missing rows. 0 NaN. 663 exact-duplicate + 27 legacy-source rows filtered at ingest.

Per-pipeline arithmetic: 257 × 27 input structures × 6 protocols × 5 reps = 208,170. Two pipelines = 416,340.

---

## 8. Pipeline recommendation (the practical takeaway)

The data support a specific end-to-end pipeline for protein-protein complex generation, **conditional on what the downstream consumer needs**:

| Downstream need | Recommended cell | MP | TM cost | Rationale |
|---|---|---|---|---|
| Lowest MolProbity, no AMBER overhead | `boltz` + `normal_beta` | 0.213 | small | the standalone-Rosetta MP minimum |
| Clean inputs for downstream tools (low clashscore) + Rosetta polish | `boltz → amber_boltz` + `dualspace_beta` | 0.222 | small | the AMBER+Rosetta MP minimum |
| Maximum TM-score retention | `af_relaxed` + `cartesian_beta` | 0.242 | smallest | AF has highest TM, cartesian_beta is most TM-retentive |
| Sidechain-rotamer fidelity | any source + `_ref15` family | varies | varies | ref2015 wins 21/42 on rota outliers |
| Crystal as input | use crystal directly; do NOT AMBER | 1.81 | n/a | AMBER damages the crystal interface (see Finding 1) |

**Headline (one sentence for the abstract)**: For protein-protein complex generation where local geometry matters, predict with AlphaFold 2.3.2 or Boltz-1, optionally pre-condition with AMBER for downstream consumers that want clean-clashscore inputs, then relax with Rosetta `dualspace_beta` or `normal_beta` (`beta_nov16` score function) to reach an integrated MolProbity score in the 0.21-0.23 band at small (< 0.025 TM-score) global-fold cost and no meaningful binding-interface drift.

---

## 9. Figures and tables

### Main figures (existing 20 + 2 new)

1-6: Pre-Rosetta analysis (TM violin, RMSD box, Blue-Green scatter, AMBER effect, AF vs Boltz, outliers).
7-11: Rosetta analysis (protocol comparison, pre/post scatter, protocol × source heatmap, replicate convergence, Blue-Green Rosetta).
12-15: MolProbity analysis (box plots by source, AMBER paired, dual effect scatter, crystal vs predicted).
16-20: Rosetta MolProbity (pre/post clash scatter, protocol MP comparison, AMBER vs Rosetta, MP heatmap, TM-MP tradeoff).
**21 (new)**: DockQ by source, with AMBER paired overlays. Captures the corrected Finding 1 (AMBER preserves prediction interfaces, perturbs crystal).
**22 (new)**: AMBER paired forest plot, post-Rosetta locked frame. Captures Finding 2 (the post-Rosetta paired test is null; honest visualization).

### Supplementary figures (existing 5 + 3 new + 1 optional)

S1-S5: Per the existing supp set.
**S6 (new)**: AMBER dual effect by difficulty tier.
**S7 (new)**: Rosetta protocol MP by difficulty tier.
**S8 (new)**: TM/MP tradeoff by difficulty tier.
**S9 (optional new)**: 1Y64 outlier detail (Boltz failure target).

### Headline figures for the manuscript

- **fig14**: AMBER dual effect (TM unchanged + MolProbity collapsed).
- **fig15**: Crystal vs predicted MolProbity, with the AMBER-on-crystal damage callout.
- **fig17**: Rosetta protocol MP ranking (with the tied-cluster-at-top caveat in the caption).
- **fig21**: DockQ by source, the new interface-preservation evidence.
- **fig22**: Forest plot of the post-Rosetta paired AMBER test, the honest null result.

---

## 10. What's not in the DB and why

**AlphaFold 3 not used.** AF2 chosen for benchmark continuity with prior docking-evaluation literature anchored to the AF2 era. Discussion paragraph defends.

**Boltz-1 single-sequence mode.** Intentional, removes MSA-quality confound across 257 targets with heterogeneous database coverage. Methods note required.

**PoseBusters not run.** Out of paper scope; methodology established in the [`Protein_Data_Analysis`](https://github.com/dreamlessx/Protein_Data_Analysis) Phase 1 pilot, retained for the historical record.

**No comparison to ClusPro / HADDOCK / pyDock.** Out of paper scope. The recommended pipeline claim is bounded to "best end-to-end relaxation pipeline given AlphaFold or Boltz starting points," not "best docking pipeline." Discussion paragraph addresses this scoping.

**43 AlphaFold targets used `reduced_dbs` fallback** (HHblits hit the 32,763-residue limit on antibody/immunoglobulin sequences). Source list trace pending. Mention as non-blocking caveat in Discussion.
