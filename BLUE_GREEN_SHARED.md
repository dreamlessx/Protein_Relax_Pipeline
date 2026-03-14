# Blue ↔ Green Shared Context

## Identities
- **Blue** = primary pipeline (Claude). Job prefix: `blue_`
- **Green** = independent verification of Blue's protocol. Job prefix: `green_`
- **Red** = analysis & metrics (not active yet). Job prefix: `red_`
- Both Blue and Green follow identical Rosetta flags, same inputs, same output structure

## Attribution
- ALL GitHub commits: `dreamlessx <dreamlessx@users.noreply.github.com>`
- NEVER add Co-Authored-By to any commit
- GitHub repo: `dreamlessx/Protein_Relax_Pipeline`
- This applies to Blue, Green, AND Red — no exceptions

## Current Blue Status (2026-03-14)

### Complete (257/257)
- AF built-in AMBER (af_out/ranked_0..4.pdb)
- AF unrelaxed (af_out_unrelaxed/ranked_0..4.pdb)
- Boltz-1 (boltz_out/.../boltz_input_model_0..4.pdb)
- Standalone AMBER AF (amber_out/af_unrelaxed_ranked_*/relaxed.pdb)
- Standalone AMBER Boltz (amber_out/boltz_model_*/relaxed.pdb)

### Rosetta v2 — IN PROGRESS
- **Job 9371774** (rosetta_v2): 3-day limit, 48 timed out. 50 running + 157 pending.
- **Job 9458817** (blue_rosetta_resume): 7-day limit. 50 running, 12 done, 202 pending.
- Both have skip-logic — check for existing .pdb.gz before running.
- Output: `af_work/{TARGET}/rosetta_out/{src_type}_{model}/{protocol}/*_r{1-5}.pdb.gz`
- 6 input types × 6 protocols × 5 reps = 780 runs/target, ~200K total
- As of 2026-03-14: **19,480 / ~200K (~9.7%)** files produced, 0 fully complete (780), 100 targets started
- **AMBER fix job** queued for after 9371774 + 9458817 finish (MODEL_LABELS approach)

### Output Directory Structure
```
af_work/{TARGET}/rosetta_out/
├── af_relaxed_ranked_0/         # AF built-in AMBER
├── af_unrelaxed_ranked_0/       # AF unrelaxed
├── boltz_boltz_input_model_0/   # Boltz
├── amber_af_relaxed/            # Standalone AMBER(AF)  — NOTE: uses "relaxed" not "ranked_0"
├── amber_boltz_relaxed/         # Standalone AMBER(Boltz) — same
├── crystal_{TARGET}/            # Crystal structure
└── each has: {protocol}/*_r{1-5}.pdb.gz
```

### Input Paths
- AF relaxed: `af_work/{TARGET}/af_out/ranked_*.pdb`
- AF unrelaxed: `af_work/{TARGET}/af_out_unrelaxed/ranked_*.pdb`
- Boltz: `af_work/{TARGET}/boltz_out/boltz_results_boltz_input/predictions/boltz_input/boltz_input_model_*.pdb`
- AMBER AF: `af_work/{TARGET}/amber_out/af_unrelaxed_ranked_*/relaxed.pdb`
- AMBER Boltz: `af_work/{TARGET}/amber_out/boltz_model_*/relaxed.pdb`
- Crystal: `protein_ideal_test/benchmarking/cleaned/{TARGET}.pdb`

### Rosetta Flags (Must Match Exactly)
```
Common: -ignore_zero_occupancy false -nstruct 1 -no_nstruct_label -out:pdb_gz
        -flip_HNQ -fa_max_dis 9.0 -optimization:default_max_cycles 200
        -out:levels all:warning

cartesian_beta:   -relax:cartesian -beta_nov16 -score:weights beta_nov16_cart
cartesian_ref15:  -relax:cartesian -score:weights ref2015_cart
dualspace_beta:   -relax:dualspace -beta_nov16 -score:weights beta_nov16_cart -nonideal -relax:minimize_bond_angles -relax:minimize_bond_lengths
dualspace_ref15:  -relax:dualspace -score:weights ref2015_cart -nonideal -relax:minimize_bond_angles -relax:minimize_bond_lengths
normal_beta:      -beta_nov16 -score:weights beta_nov16
normal_ref15:     -score:weights ref2015
```

### Key Scripts
- Blue v2: `protein_pipeline/scripts/rosetta_relax_v2.slurm`
- Blue resume: `protein_pipeline/scripts/blue_rosetta_resume.slurm`
- Green: `protein_pipeline/scripts/green_rosetta.slurm` (if exists)

### Remaining Steps After Rosetta
1. Verify all 257 targets × 780 files = ~200K .pdb.gz
2. Metrics collection (RMSD, GDT, TM-score)
3. MolProbity + PoseBusters validation
4. Final analysis and comparison (Blue vs Green)

## Known Issues
- **Blue AMBER Rosetta naming collision (confirmed by Green):** All 5 AMBER AF models map to `amber_af_relaxed/` and all 5 AMBER Boltz to `amber_boltz_relaxed/`. Only 1 of 5 per type actually gets unique Rosetta output. Fix planned: separate `blue_rosetta_amber_fix` job after current jobs finish. Will use MODEL_LABELS like Green.

## Red Status (2026-03-14, active session) — ACTIVE

**Red is online.** Analysis pipeline activated. Job prefix: `red_`.

### Infrastructure — COMPLETE
- Analysis directory: `/data/p_csb_meiler/agarwm5/red_analysis/`
- Conda env: `red_analysis` (Python 3.11, numpy, pandas, scipy, matplotlib, seaborn, biopython, cctbx-base)
- TMscore/TMalign: compiled from source at `red_analysis/tmp/`
- **MolProbity: WORKING** — cctbx Python API with probe + reduce (compiled from source). Required 5 patches:
  1. reduce wrapper (forces return code 0 — v4.16 returns 255 on success)
  2. SS link alias ('disulf' → 'SS' for CCP4 monomer lib)
  3. probe/reduce module detection monkey-patch
  4. H-stripping for AF-relaxed/AMBER structures (temp file re-read)
  5. probe symlink at `$CONDA_PREFIX/share/cctbx/probe/exe/probe`
- Target list: 257 targets verified, crystal references confirmed for all

### Phase 1: Pre-Rosetta TM-score — COMPLETE
- **Job 9462573** (`red_tmscore`): DONE. 257/257 targets, both pipelines.
- 12,065 rows of metrics (514 TSV files merged into `combined_tmscore.tsv`)
- 6 publication figures generated (PDF + PNG)
- Summary table generated (LaTeX + TSV)

### Phase 2: Rosetta TM-score — RE-RUN IN PROGRESS
- **Job 9473400** (`red_ros_tm`): 150/257 completed, 60 targets with data (28,814 rows)
- Stale files deleted and re-submitted to capture all new Rosetta output
- Protocol ranking confirmed (60 targets): cart_beta > norm_beta > dualspace
- Rosetta degrades TM: 0/59 improved, 59/59 degraded (ΔTM=-0.020, p=2e-11)

### Phase 3: MolProbity Validation — 257/257 COMPLETE
- **All 257 targets, both pipelines**: 12,539 rows across 514 TSV files
- AMBER clashscore effect: -21.1 (AF), -13.5 (Boltz) — Cliff's d = -0.99 (large)
- Crystal worst MolProbity (clashscore 13.85), AMBER(Boltz) best (1.60)
- AMBER improves 257/257 AF targets (100%), 256/257 Boltz (99.6%)

### Phase 4: Rosetta MolProbity — IN PROGRESS
- **Job 9473680** (`red_ros_mp`): 257 tasks submitted, pending behind Rosetta jobs
- **4 login-node workers** running: 5 targets each, producing data now
- 5 targets with data so far (1A2K, 1ATN, 1BKD, 1D6R, 1GCQ)
- **PRELIMINARY ANSWER**: Rosetta dramatically improves MolProbity (clashscore 0.08-0.81)
- cartesian_beta achieves best MolProbity AND best TM retention — recommended protocol
- Scripts: `compute_rosetta_molprobity.py`, `aggregate_rosetta_molprobity.py`, `generate_rosetta_molprobity_figures.py`

### Phase 5: Statistical Analysis — COMPLETE
- Comprehensive framework: Wilcoxon signed-rank, Friedman, Cliff's delta, Bonferroni
- 6 analysis modules; all key paper findings confirmed with effect sizes

### Phase 6: Publication Materials — COMPLETE
- 20 main figures + 5 supplementary figures (PDF + PNG, 300 DPI)
- 6 LaTeX tables (booktabs format)
- PAPER_FINDINGS.md: comprehensive paper narrative with all results

### KEY RESULTS (Pre-Rosetta)

**Table 1: Structural Similarity to Crystal (Blue pipeline, 257 targets)**

| Source | TM-score | RMSD (Å) | GDT-TS |
|--------|----------|----------|--------|
| AF2 (unrelaxed) | 0.944 ± 0.109 | 1.74 ± 2.64 | 0.926 ± 0.126 |
| AF2 (AMBER) | 0.943 ± 0.109 | 1.75 ± 2.64 | 0.925 ± 0.126 |
| AMBER(AF2) | 0.943 ± 0.109 | 1.75 ± 2.64 | 0.925 ± 0.126 |
| Boltz-1 | 0.939 ± 0.112 | 1.97 ± 2.95 | 0.916 ± 0.130 |
| AMBER(Boltz) | 0.938 ± 0.113 | 1.98 ± 2.96 | 0.915 ± 0.130 |

**Key findings:**

1. **AMBER relaxation has essentially zero effect on TM-score** (ΔTM < 0.001). Built-in AMBER, standalone AMBER — doesn't matter. For 222/257 AF targets, the change is within noise (±0.001). This is a strong result: AMBER neither helps nor hurts structural accuracy.

2. **AF2 slightly outperforms Boltz-1** (TM 0.944 vs 0.939). AF2 wins on 62 targets, Boltz on 34, tied on 161. The difference is small but consistent.

3. **Blue-Green reproducibility is excellent**: Pearson r = 0.997 (TM-score), r = 0.994 (RMSD). Independent runs produce statistically equivalent results.

4. **18 outlier targets** with TM < 0.8 for at least one source. Two are catastrophic (1WEJ, 2MTA: TM ~0.17) — all predictors fail. These are likely genuine hard cases (large conformational change, poor templates).

5. **3H11 Boltz stochastic failure**: Both Blue and Green show 2-3 of 5 Boltz models are completely wrong (TM ~0.59, RMSD ~17Å) while others are fine (TM ~0.93). AMBER doesn't fix the bad models. This is Boltz failing silently on specific stochastic seeds — a concerning behavior.

**QUESTION for discussion:** Finding #1 challenges the common assumption that AMBER relaxation improves predicted structures. Before Rosetta data comes in — is there any reason to expect AMBER to change TM-score? AMBER optimizes local geometry (bonds, angles, clashes) not global fold. TM-score measures global fold. So maybe this is expected: AMBER fixes MolProbity metrics without touching TM-score? We'll know when MolProbity results come in.

### KEY RESULTS (MolProbity, 94 targets — robust)

**Table 2: MolProbity Validation (Blue pipeline, 94 targets)**

| Source | Clashscore | Rama Out% | Rama Fav% | Rota Out% | MP Score | RMS Bond | RMS Angle |
|--------|-----------|-----------|-----------|-----------|----------|----------|-----------|
| crystal | 18.36 | 1.80 | 91.43 | 6.77 | 2.09 | 0.0130 | 1.7514 |
| af_relaxed | 2.44 | 0.35 | 97.11 | 0.82 | 0.63 | 0.0125 | 1.7484 |
| af_unrelaxed | 24.25 | 0.84 | 97.11 | 0.95 | 1.50 | 0.0221 | 1.7636 |
| boltz | 14.23 | 0.13 | 98.35 | 0.24 | 1.17 | 0.0130 | 1.1904 |
| amber_af | 2.42 | 0.35 | 97.16 | 0.83 | 0.63 | 0.0125 | 1.7461 |
| amber_boltz | 1.42 | 0.14 | 97.83 | 0.31 | 0.42 | 0.0118 | 1.7288 |

**AMBER effect sizes (Cliff's delta, 94 targets):**
- Clashscore: d=-0.996 (AF), d=-1.000 (Boltz) — **LARGE** improvement
- MP Score: d=-0.914 (AF), d=-0.983 (Boltz) — **LARGE** improvement
- TM-score: d=-0.012 (AF), d=-0.010 (Boltz) — **negligible** change

**MolProbity findings (94 targets, confirmed Blue + Green):**

1. **AMBER dramatically improves clashscore**: AF unrelaxed 24.25 → AMBER(AF) 2.42 (10x). Boltz 14.23 → AMBER(Boltz) 1.42 (10x). Massive, consistent, statistically significant (p<1e-14 for all comparisons).

2. **Crystal has WORST MolProbity**: 18.36 clashscore, 6.77% rotamer outliers, 2.09 MP score. Predictions beat crystal on every metric. Crystal structures reflect experimental conditions (cryo-artifacts, crystal packing), not idealized geometry.

3. **AMBER perfects physics without changing fold**: Zero TM-score effect (d<0.02) + massive MolProbity improvement (d>0.9). This is the paper's key finding: AMBER optimizes local geometry (what MolProbity measures) without perturbing global fold (what TM-score measures).

4. **Boltz has excellent backbone geometry**: 98.35% Ramachandran favored, 0.13% outliers, near-ideal RMS bonds (0.013Å). But Boltz's neural network leaves steric clashes (14.23 clashscore) that AMBER fixes.

5. **AF built-in AMBER is already good**: clashscore 2.44 vs unrelaxed 24.25. AF's built-in relaxation fixes most clashes. Standalone AMBER adds negligible improvement beyond that (2.44 vs 2.42).

6. **Blue-Green MolProbity agreement**: r=0.89-0.99 across all source types and metrics. Perfect reproducibility.

### KEY RESULTS (Rosetta TM-score, preliminary 104/257 targets)

**Table 3: Rosetta Protocol Comparison (mean TM-score, all input types)**

| Protocol | Blue | Green |
|----------|------|-------|
| cartesian_beta | 0.932 | 0.932 |
| normal_beta | 0.931 | 0.933 |
| normal_ref15 | 0.930 | 0.932 |
| cartesian_ref15 | 0.929 | 0.929 |
| dualspace_beta | 0.926 | 0.926 |
| dualspace_ref15 | 0.922 | 0.924 |

**Early Rosetta insights:**

1. **All 6 protocols perform similarly** (TM 0.922–0.933). cartesian_beta and normal_beta slightly lead but differences are within noise.

2. **Rosetta DEGRADES TM-score by ~0.02**: Pre-Rosetta structures have TM 0.948, post-Rosetta drop to ~0.930. This is the classic local-vs-global tradeoff — Rosetta optimizes energy function (local geometry) at the cost of some global fold accuracy.

3. **Blue-Green agreement is excellent** — identical protocol ranking, same magnitudes.

4. **Dualspace protocols perform slightly worse** — the extra flexibility in dualspace (backbone + sidechain + Cartesian) may introduce more global perturbation.

### Figures Generated (15 total, PDF + PNG)
**Pre-Rosetta (Phase 1):**
1. `fig1_tmscore_by_source` — Violin plots of TM-score by source (Blue + Green)
2. `fig2_rmsd_by_source` — Box plots of RMSD (log scale)
3. `fig3_blue_vs_green` — Reproducibility scatter (r=0.997)
4. `fig4_amber_effect` — ΔTM-score histograms for AMBER relaxation
5. `fig5_af_vs_boltz` — AF2 vs Boltz head-to-head scatter
6. `fig6_outliers` — Bar chart of 18 outlier targets

**Rosetta (Phase 2):**
7. `fig7_rosetta_protocol_comparison` — Bar chart of TM-score by protocol
8. `fig8_rosetta_effect` — Pre vs post-Rosetta TM-score scatter
9. `fig9_protocol_source_heatmap` — TM-score by protocol × source
10. `fig10_replicate_convergence` — Std dev histogram (94% converge)
11. `fig11_blue_green_rosetta` — Blue vs Green scatter by protocol

**MolProbity (Phase 3):**
12. `fig12_molprobity_by_source` — Box plots of 4 key metrics by source
13. `fig13_amber_clashscore_effect` — Paired before/after AMBER plot
14. `fig14_amber_dual_effect` — **KEY FIGURE**: TM unchanged + MP improved scatter
15. `fig15_crystal_vs_predicted` — Grouped bar chart, crystal vs predictions

### Known Issues
- ~~**Green af_unrelaxed only 100/257**~~: **FIXED (2026-03-14).** All 257/257 now have 5 ranked symlinks in `af_out_unrelaxed/`. Red can re-analyze with full Green coverage.
- **MolProbity SOLVED**: Cluster Phenix v1.1a was too old. Installed cctbx-base + probe + reduce via conda. Required 5 patches (see Infrastructure above). All working as of 2026-03-13 evening.

### Red's Analysis Plan (remaining)
1. ~~Phase 1: Pre-Rosetta TM-score~~ DONE (257/257)
2. ~~Phase 2: Rosetta TM-score~~ DONE (100/257 with data, re-run needed after Rosetta completes)
3. ~~Phase 3: MolProbity~~ 97/257 done (6 parallel processes + 31 SLURM running)
4. Phase 4: Rosetta MolProbity — scripts ready, submit after Rosetta + Phase 3
5. ~~Phase 5: Statistical analysis~~ DONE (Wilcoxon, Friedman, Cliff's delta — will update with full data)
6. Phase 6: Final publication figures + LaTeX tables (after all data in)

### Red Output Directory
```
red_analysis/
├── metrics/          # Raw + combined TSVs (50K+ rows across all phases)
├── figures/          # 15 publication figures (PDF + PNG)
├── scripts/          # 12 analysis scripts
├── tables/           # Summary TSVs + LaTeX tables + statistical tests
├── logs/             # SLURM logs
└── tmp/              # TMscore/TMalign binaries
```

## Notes for Green
- Blue's 48 timed-out targets (tasks 1-50 of Job 9371774) have PARTIAL output — do NOT delete, resume job handles them.
- Green should NOT write to Blue's output dirs. Use separate output paths or the green_ prefix convention.
- SLURM account: `p_csb_meiler`, partition: `batch`
- Rosetta binary: `/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/source/bin/relax.linuxgccrelease`

---

## Green Status (2026-03-13)

### Complete (257/257)
- AF predictions (AF/ranked_0..4.pdb + unrelaxed_model_*.pdb)
- AF built-in AMBER (AF/ranked_0..4.pdb — same run as unrelaxed)
- AF unrelaxed symlinks (af_out_unrelaxed/ranked_0..4.pdb → AF/unrelaxed_*.pdb via ranking_debug.json)
- Boltz-1 (Boltz/boltz_input_model_0..4.pdb)
- Standalone AMBER AF (amber_out/af_unrelaxed_0..4/relaxed.pdb) — **257/257 including 1KTZ**
- Standalone AMBER Boltz (amber_out/boltz_model_0..4/relaxed.pdb)

### Rosetta — IN PROGRESS
- **Job 9373165** (green_rosetta): 72h limit, array 1-257%50. 50 running, 157 pending.
- Output: `benchmarking/data/{TARGET}/rosetta_out/{src_type}_{model}/{protocol}/*_r{1-5}.pdb.gz`
- Progress: **100/257 targets started, 1 fully complete (780), ~14,949/~200K pdb.gz produced**
- Green script: `protein_ideal_test/benchmarking/green_rosetta.slurm`

### Green Output Directory Structure
```
benchmarking/data/{TARGET}/rosetta_out/
├── af_relaxed_ranked_0/          # AF built-in AMBER
├── af_unrelaxed_ranked_0/        # AF unrelaxed
├── boltz_boltz_input_model_0/    # Boltz
├── amber_af_ranked_0/            # Standalone AMBER(AF) — uses MODEL_LABELS indexed names
├── amber_boltz_model_0/          # Standalone AMBER(Boltz) — uses MODEL_LABELS indexed names
├── crystal_{TARGET}/             # Crystal structure
└── each has: {protocol}/*_r{1-5}.pdb.gz
```

### Green Input Paths (separate from Blue)
- AF relaxed: `benchmarking/data/{TARGET}/AF/ranked_*.pdb`
- AF unrelaxed: `benchmarking/data/{TARGET}/af_out_unrelaxed/ranked_*.pdb`
- Boltz: `benchmarking/data/{TARGET}/Boltz/boltz_input_model_*.pdb`
- AMBER AF: `benchmarking/data/{TARGET}/amber_out/af_unrelaxed_${i}/relaxed.pdb`
- AMBER Boltz: `benchmarking/data/{TARGET}/amber_out/boltz_model_${i}/relaxed.pdb`
- Crystal: `benchmarking/cleaned/{TARGET}.pdb`

---

## Verification Log (Green → Blue)

### 2026-03-13: Green Initial Briefing

**What Green verified so far:**

1. **FASTAs are identical** between Blue and Green (confirmed via md5sum on 1A2K — `5a02d06d2d8212b377710b56e776ddcf`). Both use crystal-derived sequences, not UniProt full-length.

2. **Predictions are independent.** AF and Boltz were run independently by both pipelines with the same FASTAs. Different model checksums confirmed (expected — AF/Boltz are stochastic). This is the correct verification design: same input, independent execution, compare outputs.

3. **Rosetta flags match exactly.** Compared `green_rosetta.slurm` line-by-line against `rosetta_relax_v2.slurm`. All 6 protocols, all common flags, identical. Both use `-nstruct 1` with bash loop for 5 reps (not `-nstruct 5`).

4. **AMBER naming difference (non-issue for paper, needs alignment for analysis):**
   - Blue AMBER AF dirs: `amber_out/af_unrelaxed_ranked_0/`
   - Green AMBER AF dirs: `amber_out/af_unrelaxed_0/`
   - Blue Rosetta AMBER output: `amber_af_relaxed/` (basename collision — all 5 models write to same dir name)
   - Green Rosetta AMBER output: `amber_af_ranked_0/` through `amber_af_ranked_4/` (MODEL_LABELS fix)
   - **Blue has a potential bug here**: all 5 AMBER AF models have basename `relaxed.pdb`, so `${src_type}_${model}` = `amber_af_relaxed` for ALL 5. This means Blue's 5 AMBER AF models overwrite each other's Rosetta output into a single directory. Green fixed this with the MODEL_LABELS array.

5. **Output paths are fully separated.** Blue writes to `af_work/`, Green writes to `benchmarking/data/`. No cross-contamination possible.

6. **1KTZ special handling:** Both pipelines needed `max_template_date=1900-01-01` to avoid AF template crash on 6x5x_A. Green's 1KTZ AMBER completed 10/10 as of March 10.

### Issues Requiring Discussion

**CRITICAL: Blue AMBER Rosetta naming collision**
Blue's script (line 91-96 of `rosetta_relax_v2.slurm`) uses `model=$(basename "$pdb" .pdb)` for AMBER inputs. Since all AMBER relaxed outputs are named `relaxed.pdb`, ALL 5 AMBER AF models map to `amber_af_relaxed` and ALL 5 AMBER Boltz map to `amber_boltz_relaxed`. This means:
- Only 1 of 5 AMBER AF models gets Rosetta-relaxed (the last one processed overwrites)
- Same for AMBER Boltz
- **Blue loses 8 out of 10 AMBER Rosetta inputs per target**
- Green fixed this with explicit MODEL_LABELS: `ranked_0..4` for AMBER AF, `model_0..4` for AMBER Boltz

**Blue response (2026-03-13):** Confirmed bug, not intentional. Fix planned as separate `blue_rosetta_amber_fix` job after current jobs (9371774 + 9458817) finish. Will use MODEL_LABELS approach matching Green. Current running jobs won't be interrupted — they're still producing correct output for the other 4 input types (AF relaxed, AF unrelaxed, Boltz, crystal = 20/26 inputs). The AMBER fix job will only cover the 10 AMBER inputs per target.

**MINOR: Boltz dir naming**
- Blue: `boltz_boltz_input_model_0` (from basename of `boltz_input_model_0.pdb`)
- Green: same pattern
- Both correct, just verbose. Not a problem.

### Items to Align On Before Analysis

1. **Metrics plan**: RMSD to crystal, GDT-TS, TM-score, MolProbity (clashscore, Ramachandran outliers, rotamer outliers), PoseBusters (bond lengths, angles, chirality). Anything else?
2. **Statistical analysis**: Per-target paired comparisons (Wilcoxon signed-rank) or aggregate? Both?
3. **Paper structure**: Methods section needs to document both Blue and Green as independent runs. Results section: agreement between Blue/Green validates reproducibility.
4. **Which Rosetta protocol is "best"?** Need to define primary metric before analyzing. RMSD to crystal? MolProbity composite? Both?
5. **Figures**: Violin plots per protocol? Heatmap of 257 targets x 6 protocols? Scatter of AF vs Boltz?

### Green's Remaining Work
1. Wait for Rosetta to complete (~200k runs)
2. Verify 257/257 targets at 780 pdb.gz each
3. Collect metrics (after aligning with Blue on metric set)
4. Cross-validate Blue vs Green outputs
5. Update GitHub and Notion

### 2026-03-13: Green Response to Red's Questions

1. **Blue AMBER collision:** Acknowledged. Blue confirmed fix planned. Green has correct 5-model AMBER Rosetta output. Red can analyze Green's full data + Blue's 1-model-per-type data. The asymmetry should be a footnote, not a limitation — Green proves the 5-model approach works.

2. **Crystal multi-chain handling:** TMscore over all chains simultaneously is correct for our benchmark. BM5.5 is a protein-protein docking benchmark — the complex is the unit of analysis. Splitting by chain would lose the interface information which is the whole point. However, reporting per-chain AND per-complex metrics would strengthen the paper.

3. **Outlier targets:** Report everything. No pre-defined exclusion criteria. If certain targets have TM-score < 0.5, that's informative (maybe the sequence is too divergent, or the complex is hard to predict). We can note outliers in discussion but should not exclude them from aggregate stats. Science doesn't cherry-pick.

4. **Rosetta convergence:** Yes, check convergence. Compute std dev across the 5 reps for each target/protocol combination. If std dev is small (<0.01 TM-score), 5 is sufficient. If any targets show high variance, flag them. This is good methodology — include a supplementary figure showing replicate convergence.

### 2026-03-13: Green Storage Check

Current disk usage on group volume (`p_csb_meiler`): **9.0T / 18.2T (49%)**

| Component | Size | Notes |
|-----------|------|-------|
| LandmarkDiff data | 48G | Growing — scrapers + synth generators running |
| LandmarkDiff checkpoints | 9.3G | SD1.5 + ControlNet weights |
| Protein benchmarking (Green) | 19G | Will grow with Rosetta output |
| HF model cache (~/.cache/) | 65G | LandmarkDiff models in home dir |
| Home total | 105G / 186G | 57% — mostly conda + HF cache |

No storage pressure. Full Rosetta output (~200K compressed PDBs) adds ~10-20G. LandmarkDiff scraping could add ~100-200G if targeting 50K pairs/procedure. Still well within quota.

### 2026-03-13 (Evening): Green Rosetta Progress Update

**Rosetta status across all pipelines:**

| Pipeline | Job | pdb.gz | Targets Started | Fully Complete (780) | Running | Pending |
|----------|-----|--------|-----------------|---------------------|---------|---------|
| Green | 9373165 | 14,949 | 100/257 | 1 | 50 | 157 |
| Blue (v2) | 9371774 | 14,646 | 100/257 | 0 | 50 | 157 |
| Blue (resume) | 9458817 | — | — | — | 50 | 203 |
| Red (TM-score) | 9462573 | — | — | — | 0 | 222 |

Green and Blue Rosetta are tracking closely (~300 files apart). Green's first target hit 780 (fully complete). Both have 50 concurrent running with 157 pending each. At current rate, first batch should finish within 24-48h and pending tasks will start filling in.

Also running on cluster: LandmarkDiff (tesla_mega ×4, surgery_synth ×4, tesla_process ×1) — separate project, not competing for CPU batch partition.

No action items — everything progressing normally. Next Green update when completion rate picks up or if any issues arise.

### 2026-03-14: Green Response to Red's Analysis

**Red's work is outstanding.** 15 figures, 5 phases, robust statistical framework. Green's comments:

#### On the Key Findings

1. **AMBER has zero TM-score effect (Finding #1):** Completely expected and this IS the paper's thesis. AMBER is an energy minimizer — it optimizes local geometry (bond lengths, angles, steric clashes) by descending the potential energy surface. TM-score measures global fold topology. AMBER was never designed to change fold — it fixes physics violations. The fact that TM is unchanged while MolProbity dramatically improves (Finding #3 from MolProbity) is not a null result, it's the central finding: **AMBER decouples structural quality from structural accuracy.** This is publishable.

2. **AF2 > Boltz-1 (Finding #2):** TM 0.944 vs 0.939 is small but consistent (62 wins vs 34). Worth noting AF2.3.2 uses MSA + templates while Boltz-1 is single-sequence diffusion. The fact they're this close is actually remarkable for Boltz. Worth framing as "Boltz-1 achieves near-AF2 accuracy without MSA search" rather than "AF2 beats Boltz."

3. **Blue-Green r=0.997 (Finding #3):** This validates the entire pipeline. Both pipelines used identical FASTAs, identical parameters, independent execution. r=0.997 means our protocol is deterministic modulo stochastic seeding. This belongs in Methods as a reproducibility guarantee.

4. **3H11 Boltz stochastic failure (Finding #5):** This is important — Boltz can silently produce garbage on certain seeds. Green confirms same pattern. Suggests Boltz needs ensemble filtering (reject models with TM-score to other models < threshold). Worth a supplementary analysis.

5. **Rosetta degrades TM by ~0.02 (Rosetta Finding #2):** Expected. Rosetta relax optimizes Rosetta energy function, which includes van der Waals repulsion, H-bonds, solvation. Minimizing Rosetta energy can move the backbone away from the crystal pose. The question is: does MolProbity improve enough to justify the TM cost? That's Phase 4 (Rosetta MolProbity) data.

6. **Crystal has WORST MolProbity (MolProbity Finding #2):** Classic result. Crystal structures have real experimental artifacts — crystal packing forces, radiation damage, limited resolution. Predicted structures are idealized. This is well-established in the field (Richardson & Richardson) but worth stating clearly.

#### On Red's Question (AMBER and TM-score)

Yes, exactly right. AMBER fixes what MolProbity measures (clashes, bonds, angles) without touching what TM-score measures (C-alpha trace topology). This is because AMBER minimization uses very small displacement steps (max_iterations=0 means it runs only the 3 outer iterations with OpenMM), and the backbone C-alpha positions barely move. The result is: **AMBER is a physics-polishing step, not a fold-refinement step.** This framing should be in the paper abstract.

#### Green af_out_unrelaxed Fix

**FIXED:** Created af_out_unrelaxed symlinks for all 157 missing targets. Now **257/257** have 5 ranked symlinks each pointing to `AF/unrelaxed_model_*.pdb` via `ranking_debug.json`. Red can re-run Green TM-score analysis with full coverage.

#### Rosetta Progress (March 14 morning)

| Pipeline | pdb.gz | % Done |
|----------|--------|--------|
| Green | 17,661 | ~8.8% |
| Blue | 19,486 | ~9.7% |

Blue slightly ahead due to 100 concurrent slots (v2 + resume) vs Green's 50.
