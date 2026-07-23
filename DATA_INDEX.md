# Data Index

**"I want X. Where do I get it?"** Every artifact in the BM5.5 docking-relaxation benchmark, where it lives, and one command to fetch it.

## The story in one paragraph

We benchmark structure-prediction and relaxation pipelines on the 257-complex BM5.5 protein-protein docking dataset. AlphaFold 2.3.2 (5 ranked + 5 unrelaxed) and Boltz-1 v0.4.1 (5 single-sequence) predictions feed a relaxation matrix: AMBER as a pre-conditioning step, then Rosetta 3.15 under 6 protocols (cartesian / dualspace / normal × beta_nov16 / ref2015) at 5 replicates each. Two independent pipelines (Blue, primary; Green, matched-parameters re-run) produce 416,340 Rosetta MolProbity rows, 416,340 Rosetta energy rows, 105,550 TM scores, 13,364 pre-Rosetta MolProbity rows, 13,878 DockQ rows, and 257 targets with full metadata in a single locked SQLite DB at snapshot 2026-04-27a (qc_status = pass). Five findings emerge: AMBER fixes prediction defects without meaningfully changing global fold or interface; AMBER on crystal produces a small perturbation comparable to refinement uncertainty; post-Rosetta paired AMBER effect is small and below practitioner noise floor (16 of 16 cells below MDE); dualspace_beta and normal_beta tie at top of integrated MolProbity (within 0.01 MP units), with beta_nov16 dominating ref2015 on 40-42 of 42 triples; crystal MolProbity reflects score-function idealization rather than crystallographer failure.

Last verified: 2026-04-28. Every fact table at 100% coverage; no audit-log padding.

---

## Cheat sheet

```bash
cd ~/proteins-workspace
make sync-all          # everything from ACCRE: figures + tables + raw metrics + DB + pptx
make sync-figures      # 64 PNG + 64 PDF (37 MB)
make sync-tables       # 16 TSV + 17 LaTeX aggregate tables (1.8 MB) → aggregates/
make sync-metrics      # 5 raw row-level TSVs (~100 MB) → aggregates/raw/
make sync-db           # bm55.sqlite via sqlite3 .backup (~225 MB) → db/
make sync-pptx         # 3 BM55 .pptx → briefs/pptx/
make accre-status      # squeue + recent ACCRE log tail
```

All sync scripts are atomic (`rsync --partial --delay-updates`) and idempotent.

---

## "I want X" matrix

### Final figures (paper-ready PNG + PDF)
- **Local:** `figures/` (after `make sync-figures`)
- **ACCRE source:** `accre:/data/p_csb_meiler/agarwm5/red_analysis/figures/`
- **Repo:** `repos/Protein_Relax_Pipeline/red_analysis/figures/` (snapshot, may lag ACCRE)
- **Inventory:** 64 PNG + 64 PDF. fig1-fig22 main + figS1-figS5 supplementary + audit_*. Most have `_blue` and `_green` variants.
- **Naming convention:** `fig<N>_<name>.png` (combined), `fig<N>_<name>_blue.png` / `_green.png` (per-pipeline), `figS<N>_<name>.png` (supp), `audit_*` (not paper figs).
- **Headline figs:** fig4 (AMBER dual effect), fig14 (AMBER PI slide), fig15 (crystal worst MP), fig17 (Rosetta protocol MP), fig18 (AMBER vs Rosetta), fig20 (TM/MP tradeoff), figS2 (violins), figS3 (per-source).

### Aggregate tables (TSV + LaTeX, paper-ready)
- **Local:** `aggregates/` (after `make sync-tables`)
- **ACCRE source:** `accre:/data/p_csb_meiler/agarwm5/red_analysis/tables/`
- **Repo:** `repos/Protein_Relax_Pipeline/red_analysis/tables/` (committed)
- **Inventory:** 16 TSV + 17 TeX. Includes `amber_paired_table2.tsv` (16 cells main paired test), `amber_paired_rotamer.tsv` (18 cells rotamer extension), `per_target_mp_means.tsv`, `pairwise_wilcoxon.tsv`, `summary_by_source.tsv`, `rosetta_protocol_summary.tsv`, plus LaTeX `table1_*.tex` through `table6_*.tex`.

### Raw row-level metrics (per-structure measurements)
- **Local:** `aggregates/raw/` (after `make sync-metrics`)
- **ACCRE source:** `accre:/data/p_csb_meiler/agarwm5/red_analysis/metrics/`
- **Inventory:**
  | File | Rows | What |
  |---|---|---|
  | `combined_molprobity.tsv` | 13,344 | Pre-Rosetta MP (initials: AF / Boltz / AMBER variants / crystal MolProbity before Rosetta). 20 Blue crystal rows backfilled in DB from Green source via Stage C. |
  | `combined_tmscore.tsv` | 12,850 | Pre-Rosetta TM-score (initials vs crystal); post green af_unrelaxed backfill |
  | `combined_rosetta_molprobity.tsv` | 418,060 raw / 416,340 in lock | Post-Rosetta MP (the paper's primary fact table) |
  | `combined_rosetta_tmscore.tsv` | 92,700 | Post-Rosetta TM-score |
  | `combined_rosetta_energy.tsv` | 416,370 raw / 416,340 in DB | Rosetta total_score + per_residue_energy. 30 legacy src_type rows filtered at ingest. |
- **Schema (header rows):** `target | pipeline | source | model_idx | clashscore | rama_outliers | rama_favored | rota_outliers | molprobity_score | cbeta_outliers | rms_bonds | rms_angles` (MP) and `target | pipeline | source | model_idx | rmsd | tmscore | gdtts | gdtha | aligned_len | seq_len` (TM).

### SQLite DB (queryable, locked snapshot)
- **Local:** `db/bm55.sqlite` (after `make sync-db`, ~225 MB VACUUMed)
- **ACCRE source:** `accre:/data/p_csb_meiler/agarwm5/red_analysis/db/bm55.sqlite`
- **Schema:** `repos/Protein_Relax_Pipeline/db/sql/schema.sql`
- **Builders:** `repos/Protein_Relax_Pipeline/db/scripts/build_db.py` (locks `rosetta_metrics` from per-target TSVs) + `build_db_supplements.py` (idempotent companion that loads the supplemental tables)
- **Release:** [`db-2026-04-27a-supp`](https://github.com/dreamlessx/Protein_Relax_Pipeline/releases/tag/db-2026-04-27a-supp) on the primary repo: 6 assets (DB + 5 raw TSVs)
- **Snapshot:** 2026-04-27a (qc_status=pass)
- **Tables (current state):**
  | Table | Rows | Notes |
  |---|---|---|
  | `rosetta_metrics` | **416,340** | Locked. Authoritative post-Rosetta MP. |
  | `targets` | 257 | All metadata populated: difficulty (162R/60M/35D from zlab), category (8 zlab codes: AA/AS/EI/ER/ES/OG/OR/OX), n_chains (605 total), n_residues (122,966 total). 4 non-standard flagged (BAAD, BOYV, BP57, CP57) with `parent_pdb_id` (3AAD_A:B, 1OYV_B:I, 3P57_AB:P, 3P57_CD:P). |
  | `pipelines` | 2 | blue, green |
  | `sources` | 7 | af_relaxed, af_unrelaxed, amber_af, amber_boltz, amber_crystal, boltz, crystal |
  | `protocols` | 6 | cartesian/dualspace/normal × beta_nov16/ref2015 |
  | `prerosetta_metrics` | **13,364** | 13,344 from `combined_molprobity.tsv` + 20 Stage C backfilled from Green crystal MP (PDB-identical, MolProbity deterministic). |
  | `tm_scores` | **105,550** | 12,850 pre-Rosetta (post green af_unrelaxed backfill of 157 missing targets) + 92,700 post-Rosetta. PK uses `is_post_rosetta` flag; pre-Rosetta has NULL protocol_id and rep. |
  | `rosetta_energy` | **416,340** | Total_score + per_residue_energy per (target, pipeline, src_type, protocol, rep) cell. 100% coverage of rosetta_metrics, 0 orphans, 0 gaps. Extracted via patched `extract_rosetta_energy.py` (relax.fasc + sidecar score_*.sc + PDB POSE_ENERGIES_TABLE fallback). |
  | `dockq_metrics` | **13,878** | DockQ + i-RMSD + l-RMSD + f_nat per (target, pipeline, src_type) cell on the 27 input structures. 100% symmetric across Blue and Green after 257-row Blue amber_crystal backfill. |
  | `qc_quarantine` | **0** | Empty; the snapshot is fully consistent. |
  | `build_runs` | 1 | Single locked snapshot 2026-04-27a, qc_status=pass. |
  | views: `v_cell_summary`, `v_per_target_mp_means` | n/a | Pre-built convenience views |

  **Status (2026-04-29):** Every fact table at 100% coverage post data-gap fixes (157 green af_unrelaxed TM rows + 257 Blue amber_crystal DockQ rows backfilled; 1Y64 Boltz outlier identified). Loader: `db/scripts/build_db_supplements.py` (idempotent, applies schema migrations for `rosetta_energy` + `dockq_metrics` + `targets.parent_pdb_id` on first run; wipe-and-reload on each run). Single citation point for the manuscript: snapshot 2026-04-27a.

### Per-target output (Blue pipeline, individual target metrics)
- **ACCRE source:** `accre:/data/p_csb_meiler/agarwm5/red_analysis/metrics/molprobity_blue_<TARGET>.tsv` and `molprobity_green_<TARGET>.tsv` (one per target × per pipeline)
- **Use case:** debugging a single target. Not synced by default (would be 514 small files); pull individually via `ssh accre 'cat /data/p_csb_meiler/agarwm5/red_analysis/metrics/molprobity_blue_1ACB.tsv'`.

### PowerPoint slides (PI-ready)
- **Local:** `briefs/pptx/` (after `make sync-pptx`)
- **ACCRE source:** `accre:/data/p_csb_meiler/agarwm5/red_analysis/BM55_Relaxation_Benchmark_*.pptx`
- **Inventory:** Blue (16 slides), Green (16 slides), combined.

### Crystal PDBs (the BM5.5 originals)
- **Repo (canonical):** `repos/Protein_Relax_Pipeline/cleaned/` (257 cleaned crystal PDBs, homo-multimer dedup applied)
- **ACCRE Green source:** `accre:/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/cleaned/`
- **FASTAs:** `repos/Protein_Relax_Pipeline/fasta/` (257 crystal-derived FASTAs)
- **Note:** Blue's original `protein_pipeline/cleaned/` has been removed from ACCRE; Green's cleaned crystal PDBs are byte-identical to Blue's per empirical verification on 1A2K, 7CEI, 1ACB.

### AlphaFold + Boltz predictions (the upstream models)
- **ACCRE source:** `accre:/data/p_csb_meiler/agarwm5/protein_pipeline/predictions/alphafold/` and `predictions/boltz/`
- **Per target:** 5 ranked AF (AMBER-relaxed) + 5 unrelaxed AF + 5 Boltz
- **Total upstream input structures fed into Rosetta:** 6,939 = 257 crystal + 1,285 AF ranked + 1,285 AF unrelaxed + 1,285 Boltz + 257 amber_crystal + 1,285 amber_af + 1,285 amber_boltz

### Rosetta outputs (the 416,340 PDBs)
- **ACCRE source:** `accre:/data/p_csb_meiler/agarwm5/protein_pipeline/{blue,green}_relax/...` (mirrored under protein_ideal_test/)
- **Per target per pipeline:** 27 input models × 6 protocols × 5 reps = 810 .pdb.gz files
- **Not synced locally by default** (208,170 × 2 pipelines = 416,340 small files would take hours).
- **Re-fetch a single (target, source, protocol) cell:** `rsync accre:/data/p_csb_meiler/agarwm5/protein_pipeline/blue_relax/1ACB/<src_type>/<protocol>/r*.pdb.gz ./tmp/`

### Manuscript (drafted from scratch, post-deletion 2026-04-28)
- **Vault:** `Dreamless_Machine/03-Research/Protein/Docs/manuscript/` (v2 workspace; v1 cleared 2026-04-28)
- **Workspace:** `manuscript/` (workspace-side new sections, response-to-reviewer drafts, build artifacts)
- **Style:** Concise academic register with quantitative, data-led results. No iteration logs in frontmatter; use the `verifier_history` list.
- **Target:** Proteins: Structure, Function, and Bioinformatics (250-word abstract limit).

### PI briefs and weekly status
- **Workspace:** `briefs/<YYYY-MM-DD>_<topic>.md`
- **Vault PI brief reference:** `pi_review_2026-04-26/PI_REVIEW_BRIEF.md` (DELIVERED, post-meeting note added)
- **Template:** model new briefs on the canonical PI brief structure (1-line narrative → 3 headline findings → method defensibility → open items → file index).

### Vault landmark docs
- `Plan.md`: strategic plan
- `Status.md`: live project status
- `Notes.md`, `Ideas.md`: Damien / Ary owned, free-form
- `Protein.md`: hub dashboard with dataview queries
- `Docs/{Architecture, Overview, Decisions, Glossary, References, Changelog, Status}.md`, landmark docs
- `Docs/audits/`, historical audit snapshots (citation audit, data completeness iter48, interim rosetta 2026-04-18)
- `pi_review_2026-04-26/`, PI brief + files index

---

## Repo map (GitHub, public)

| Repo | URL | Role |
|---|---|---|
| `Protein_Relax_Pipeline` | github.com/dreamlessx/Protein_Relax_Pipeline | BLUE pipeline + canonical analysis (`red_analysis/scripts/`) + DB schema/build + manuscript code. **Primary reference repo.** |
| `Protein_Ideal` | github.com/dreamlessx/Protein_Ideal | GREEN pipeline + dataset prep |
| `Protein_Data_Analysis` | github.com/dreamlessx/Protein_Data_Analysis | MolProbity validation toolchain (Phase 1 pilot, 20 proteins) |

All three are clean at upstream HEAD as of 2026-04-29. DB at 100% coverage: `rosetta_metrics` 416,340 + `prerosetta_metrics` 13,364 + `tm_scores` 105,550 + `rosetta_energy` 416,340 + `dockq_metrics` 13,878 + `targets` 257 (full metadata + parent_pdb_id) + `qc_quarantine` 0. Loader: `Protein_Relax_Pipeline/db/scripts/build_db_supplements.py`. Energy extractor: `Protein_Relax_Pipeline/red_analysis/scripts/extract_rosetta_energy.py`. DockQ extractor: `Protein_Relax_Pipeline/red_analysis/scripts/compute_dockq_inputs.py`. Live SQLite + 6 raw TSVs in the `db-2026-04-27a-supp` GitHub Release.

---

## Number reference (verified against locked DB)

```
257 BM5.5 complexes      = 162 rigid + 60 medium + 35 difficult (incl. 4 non-standard: BAAD, BOYV, BP57, CP57)
605 chains, 122,966 res  = post-strip across 257; 119 homo-multimers deduped to unique seqs; 41 his-tags removed
6,939 input structures   = 257 crystal + 1,285 AF ranked + 1,285 AF unrelaxed + 1,285 Boltz + 257 amber_crystal + 1,285 amber_af + 1,285 amber_boltz
416,340 Rosetta MP rows  = 5 × 77,100 (5-model sources) + 2 × 15,420 (1-model sources). Blue (208,170) + Green (208,170); the DB unifies under one snapshot.
208,170 per pipeline     = 257 × 810 outputs/target = 257 × 27 × 6 × 5
13,364 pre-Rosetta MP    = post-Stage-C. 13,344 from combined TSV + 20 Stage C Blue crystal backfill. 6,682 Blue + 6,682 Green; symmetric.
12,850 pre-Rosetta TM    = initials vs crystal TM-score (post green af_unrelaxed backfill of 157 missing targets)
92,700 post-Rosetta TM   = Rosetta outputs vs crystal TM-score
416,340 Rosetta energy   = total_score + per_residue_energy per Rosetta cell, 1:1 with rosetta_metrics
1,720 filtered at lock   = 663 exact-duplicate + 27 legacy-source + 1,030 secondary dedup (re-aggregation duplicates from per-target TSVs)
6,820 (Phase 1 pilot)    = 20 proteins × 341 structures (Protein_Data_Analysis only, pre-full-benchmark)
```

---

## Five findings (manuscript-ready, source-of-truth)

See `repos/Protein_Relax_Pipeline/red_analysis/PAPER_FINDINGS.md` for the canonical version with full numbers and figure pointers.

1. **AMBER fixes prediction defects without meaningfully changing global fold or interface.** Clashscore Cliff's d = -0.99 (mean drop 13-21 per target). TM mean Δ < 0.001 absolute; DockQ mean Δ < 0.0002 (below DockQ practitioner threshold). 257/257 AF + 256/257 Boltz improved. 1Y64 is the one Boltz outlier where AMBER raises MP by 0.46.
2. **AMBER on crystal: small structural perturbation comparable to crystallographic refinement uncertainty.** Median Δ i-RMSD 0.246 Å (within Luzzati 0.2-0.3 Å envelope at BM5.5 ~2.2 Å median resolution); Δ DockQ -0.022 below practitioner decision threshold. Cliff's d = -1.0 mathematically tautological. Skip AMBER on crystal sources held as unmodified reference.
3. **Post-Rosetta paired AMBER effect: small, directionally consistent, below practitioner noise floor.** 0/16 paired-t cells, 2/16 Wilcoxon (Blue amber_boltz dualspace_beta, mean Δ MP -0.005). Magnitude below n=257 paired-t MDE. AMBER as geometry pre-conditioner, not Rosetta-MP scoring enhancer.
4. **`dualspace_beta` and `normal_beta` tie at the top of integrated MolProbity** (within 0.01 MP units). beta_nov16 dominates ref2015 on 40-42 of 42 triples for MP/clash/Rama; ref2015 wins on rotamer outliers (21/42). cartesian_beta is the TM-retentive runner-up.
5. **Crystal MolProbity calibration.** Clashscore 13.85 vs predictions 0.68-2.82. Frame as score-function idealization (predictions+AMBER+Rosetta optimize toward score-function ideal; crystals reflect refinement constraints).
