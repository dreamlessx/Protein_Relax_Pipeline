# Red Analysis: Paper Findings Summary
## BM5.5 Protein Structure Relaxation Benchmark

**Last updated: 2026-03-22**
**Status: Phase 2 COMPLETE (257 targets, 93K rows), Phase 4 expanding (160 targets, 85K rows)**

---

## 1. AMBER Relaxation: Universal Local Improvement, Zero Global Effect

### The Dual Effect (Key Paper Finding)

AMBER relaxation simultaneously:
- Has **negligible effect on fold accuracy** (TM-score Cliff's d = -0.01)
- Has **massive effect on local geometry** (Clashscore Cliff's d = -0.99)

This is not a tradeoff — it's a free lunch. AMBER fixes what TM-score can't see.

### Numbers (257 targets, Blue pipeline)

| Comparison | ΔTM-score | TM Effect Size | ΔClashscore | Clash Effect Size |
|------------|-----------|-----------------|-------------|-------------------|
| AF unrelaxed → AMBER(AF) | -0.0004 | d=-0.012 (negligible) | -21.1 | d=-0.992 (large) |
| Boltz → AMBER(Boltz) | -0.0009 | d=-0.010 (negligible) | -13.5 | d=-0.989 (large) |

### Universality
- AMBER improves MolProbity for **257/257 AF targets (100%)**
- AMBER improves MolProbity for **256/257 Boltz targets (99.6%)**
- This is not sampling bias — it's a physical law: energy minimization always reduces steric clashes

### Best Source for Local Geometry
- AMBER(Boltz) achieves best MolProbity score for **77% of targets**
- AF relaxed (built-in AMBER) is already good — standalone AMBER adds negligible improvement (2.82 → 2.82 clashscore)
- The big win is AMBER on Boltz: clashscore drops from 15.1 to 1.6

### Why This Matters for the Paper
This finding challenges the assumption that AMBER relaxation is optional. While AMBER doesn't help if your metric is TM-score, it dramatically improves the physics. For downstream applications (docking, MD simulation, drug design), local geometry matters. AMBER should be standard practice.

---

## 2. Crystal Structures Have the Worst MolProbity

### A Provocative Finding

| Source | Clashscore | Rota Outliers % | MP Score |
|--------|-----------|-----------------|----------|
| Crystal | 13.85 | 6.15 | 1.81 |
| AF relaxed | 2.82 | 0.89 | 0.70 |
| AMBER(Boltz) | **1.60** | **0.36** | **0.45** |

Crystal structures have **5x worse clashscore** and **7x more rotamer outliers** than predictions.

### Explanation
- Crystal structures reflect experimental conditions (cryogenic temperature, crystal packing)
- Predicted structures are optimized by neural networks that learned ideal geometry
- AMBER further perfects the physics
- This doesn't mean predictions are "better" — crystal is ground truth for global fold
- But for local geometry, predictions with AMBER are superior

### Paper Narrative
Frame this carefully: predictions beat crystal on MolProbity because they're computationally idealized. Crystal structures capture real experimental constraints (contacts, disorder, alternative conformations). This is expected, not surprising — but quantifying it across 257 targets is novel.

---

## 3. Rosetta Relaxation: Consistent TM-score Degradation

### All 6 Protocols Degrade TM-score (Medium Effect)

| Source | ΔTM-score | Effect Size | p-value | Improved/Degraded |
|--------|-----------|-------------|---------|-------------------|
| AF relaxed | -0.019 | d=-0.36 (medium) | 8e-17 | 2/97 |
| AF unrelaxed | -0.025 | d=-0.43 (medium) | 7e-10 | 1/53 |
| Boltz | -0.024 | d=-0.36 (medium) | 4e-06 | 0/19 |

### Protocol Ranking
Friedman test: p < 1e-33 (highly significant protocol differences)

| Protocol | Blue TM | Green TM |
|----------|---------|----------|
| Cartesian β | 0.932 | 0.932 |
| Normal β | 0.931 | 0.933 |
| Normal REF15 | 0.930 | 0.932 |
| Cartesian REF15 | 0.929 | 0.929 |
| Dualspace β | 0.926 | 0.926 |
| Dualspace REF15 | 0.922 | 0.924 |

**Top tier**: cartesian_beta ≈ normal_beta (statistically indistinguishable)
**Bottom tier**: dualspace protocols (significantly worse)

### Phase 4 Robust Results (100 targets, 43K rows, both pipelines)

**Rosetta DRAMATICALLY improves MolProbity — MORE than AMBER.**

| Protocol | Clashscore | MP Score | Rota Outliers % | Rama Favored % |
|----------|-----------|----------|-----------------|----------------|
| dualspace_beta | **0.650** | **0.215** | 0.015 | 98.66 |
| cartesian_beta | 0.673 | 0.223 | 0.017 | 98.61 |
| normal_beta | 0.685 | 0.242 | 0.028 | 98.21 |
| normal_ref15 | 1.046 | 0.335 | 0.026 | 98.04 |
| dualspace_ref15 | 1.294 | 0.358 | 0.016 | 98.57 |
| cartesian_ref15 | 1.416 | 0.383 | 0.020 | 98.48 |

**Key insight**: Beta energy function is the primary differentiator (0.65-0.69 vs 1.0-1.4).
Sampling method (cartesian/dualspace/normal) matters less within same energy function.

**Pre vs Post Rosetta (Blue, 99 targets):**
- AF unrelaxed → Rosetta: clashscore 22.9 → 1.03 (**96% reduction**, d=-1.00)
- Boltz → Rosetta: clashscore 12.5 → 0.87 (**93% reduction**, d=-1.00)
- Crystal → Rosetta: clashscore 16.2 → 0.95 (**94% reduction**, d=-1.00)
- AF relaxed → Rosetta: clashscore 2.66 → 0.97 (**64% reduction**, d=-0.58)
- Rosetta eliminates rotamer outliers: 0.80% → 0.02% (d=-0.96)

**Rosetta vs AMBER direct comparison (99 targets):**
- AMBER(AF) clashscore: 2.66 → Rosetta(AF): 0.97 (Δ=-1.69, **Rosetta 64% better**)
- AMBER(AF) MP score: 0.66 → Rosetta(AF): 0.29 (Δ=-0.37)
- AMBER(AF) rota outliers: 0.83% → Rosetta(AF): 0.02% (Δ=-0.81)

**Per-target best protocol (clashscore):**
- normal_beta: 47% | dualspace_beta: 28% | cartesian_beta: 21%
- Beta protocols collectively win 96% of targets

**Blue-Green agreement**: Protocol rankings match perfectly between pipelines.

**The verdict (Outcome A, confirmed at scale):** Rosetta improves MolProbity MORE than AMBER, but at a TM-score cost of ~0.02. The beta energy function is the critical factor. For applications requiring perfect local geometry (docking, drug design), Rosetta post-processing adds significant value.

**Recommended protocol: cartesian_beta** — best TM retention (Phase 2) with near-best MolProbity (Phase 4). While normal_beta wins slightly more targets on MolProbity, cartesian_beta provides the best overall tradeoff.

---

## 4. Reproducibility (Blue-Green Agreement)

### Everything Reproduces

| Metric | Pearson r | n |
|--------|-----------|---|
| Pre-Rosetta TM | 0.997 | 1128 |
| Pre-Rosetta RMSD | 0.994 | 1128 |
| Rosetta TM | 0.988-0.992 | 96-105 |
| Clashscore | 0.867-0.991 | 100-257 |
| MP Score | 0.941-0.984 | 100-257 |

Blue and Green are independent implementations of the same protocol. Their agreement validates:
1. The protocol specification is unambiguous
2. Results are not implementation artifacts
3. Stochastic variation (AF/Boltz/Rosetta seeds) averages out across 257 targets

---

## 5. Prediction Method Comparison

### AF2 vs Boltz-1

| Metric | AF2 (relaxed) | Boltz-1 | Winner |
|--------|--------------|---------|--------|
| TM-score | 0.943 | 0.939 | AF2 (54% of targets) |
| Clashscore | 2.82 | 15.09 | AF2 (5x better) |
| Rama Favored % | 96.99 | 98.38 | Boltz (better backbone) |
| Rota Outliers % | 0.89 | 0.23 | Boltz (4x fewer) |
| MP Score | 0.70 | 1.16 | AF2 (built-in AMBER) |

**Nuance**: Boltz has better backbone geometry but worse clashscore. After AMBER, Boltz catches up:
- AMBER(Boltz) MP Score: **0.45** vs AF relaxed: 0.70

AF2 wins on global fold (TM), Boltz wins on backbone (Ramachandran), AMBER(Boltz) wins on overall local geometry.

### Outliers
- 3 catastrophic failures (TM < 0.5): 1WEJ, 2MTA, 1DQJ — all methods fail
- 6 targets where AF >> Boltz (ΔTM > 0.1)
- 5 targets where Boltz >> AF (ΔTM > 0.1)
- 3H11: Boltz stochastic failure (2-3 of 5 models completely wrong)

---

## 6. Data Scale

| Dataset | Rows | Targets | Pipelines | Sources |
|---------|------|---------|-----------|---------|
| Pre-Rosetta TM | 12,065 | 257 | 2 | 5 |
| MolProbity | 12,539 | 257 | 2 | 6 |
| Rosetta TM | 92,700 | 257 | 2 | 6 |
| Rosetta MolProbity | 84,825 | 160 | 2 | 6 |
| **Total** | **~202K** | **257** | **2** | **6** |

---

## 7. Figures for Paper

### Main Figures (20)
1-6: Pre-Rosetta analysis (TM violin, RMSD box, Blue-Green scatter, AMBER effect, AF vs Boltz, outliers)
7-11: Rosetta analysis (protocol comparison, pre/post scatter, protocol×source heatmap, convergence, Blue-Green Rosetta)
12-15: MolProbity analysis (box plots by source, AMBER paired, **dual effect scatter**, crystal vs predicted)
16-20: Rosetta MolProbity (pre/post clash scatter, protocol MP comparison, **AMBER vs Rosetta**, MP heatmap, **tradeoff plot**)

### Supplementary Figures (5)
S1: All 7 AMBER MolProbity metrics
S2: MolProbity violin plots
S3: Rosetta TM by input type
S4: MolProbity metric correlation matrix
S5: AF vs Boltz MolProbity scatter

### Key Figures
- **fig14 (AMBER Dual Effect)**: Two-panel scatter showing TM (no change) vs MolProbity (improvement)
- **fig18 (AMBER vs Rosetta)**: Bar chart comparing all methods' MolProbity — shows Rosetta beats AMBER
- **fig20 (The Tradeoff)**: ΔTM vs ΔClashscore scatter — each target as a point, AMBER as star reference

---

## 8. Remaining Work

1. **Phase 2**: 257 targets (92.7K rows) — ALL TARGETS COMPLETE
2. **Phase 4**: 160 targets (84.8K rows, both pipelines) — 16 workers expanding coverage
3. **Rosetta relaxation**: Blue 251/257, Green 257/257 COMPLETE
4. **Final analysis**: Phase 4 workers will reach 200+ targets within 24h
5. **Write paper**: All findings documented, 25 figures + 6 LaTeX tables + 11 TSVs ready
6. **PowerPoint**: 16-slide presentation generated (`BM55_Relaxation_Benchmark.pptx`)
7. **GitHub**: `data-analysis` branch pushed (commit 0b90e30f)
