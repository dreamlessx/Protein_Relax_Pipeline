# Data Index

**"I want X. Where do I get it?"** — every artifact in this project, where it lives, and one command to fetch it.

Last verified: 2026-04-28 (every fact table at 100% coverage; no audit-log padding). Data lock 2026-04-27a: 416,340 Rosetta MP + 13,364 prerosetta MP + 104,765 TM scores + **416,340 Rosetta energy (100%)** + 257 targets full metadata + 0 quarantine rows. qc_status=pass.

---

## Cheat sheet

```bash
cd ~/proteins-workspace
make sync-all          # everything from ACCRE: figures + tables + raw metrics + DB + pptx
make sync-figures      # 64 PNG + 64 PDF (37 MB)
make sync-tables       # 16 TSV + 17 LaTeX aggregate tables (1.8 MB) → aggregates/
make sync-metrics      # 5 raw row-level TSVs (~65 MB) → aggregates/raw/
make sync-db           # bm55.sqlite via sqlite3 .backup (115 MB) → db/
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
  | `combined_molprobity.tsv` | 13,344 | **Pre-Rosetta MP** (the initials: AF/Boltz/AMBER variants/crystal MolProbity before Rosetta) |
  | `combined_tmscore.tsv` | 12,065 | **Pre-Rosetta TM-score** (initials vs crystal) |
  | `combined_rosetta_molprobity.tsv` | 418,060 raw / 416,340 in lock | Post-Rosetta MP (the paper's primary table) |
  | `combined_rosetta_tmscore.tsv` | 92,700 | Post-Rosetta TM-score |
  | `combined_rosetta_energy.tsv` | 184,352 | Rosetta energy values |
- **Schema (header rows):** `target | pipeline | source | model_idx | clashscore | rama_outliers | rama_favored | rota_outliers | molprobity_score | cbeta_outliers | rms_bonds | rms_angles` (MP) and `target | pipeline | source | model_idx | rmsd | tmscore | gdtts | gdtha | aligned_len | seq_len` (TM).

### SQLite DB (queryable, locked snapshot)
- **Local:** `db/bm55.sqlite` (after `make sync-db`, 115 MB)
- **ACCRE source:** `accre:/data/p_csb_meiler/agarwm5/red_analysis/db/bm55.sqlite` (144 MB live)
- **Schema:** `repos/Protein_Relax_Pipeline/db/sql/schema.sql`
- **Builder:** `repos/Protein_Relax_Pipeline/db/scripts/build_db.py`
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
  | `tm_scores` | **104,765** | 12,065 pre-Rosetta + 92,700 post-Rosetta. PK uses `is_post_rosetta` flag; pre-Rosetta has NULL protocol_id and rep. |
  | `rosetta_energy` | **416,340** | Schema-extension table (added 2026-04-28). Total_score + per_residue_energy per (target, pipeline, src_type, protocol, rep) cell. **100% coverage of rosetta_metrics**, 0 orphans, 0 gaps. Extracted via patched `extract_rosetta_energy.py` (canonical `relax.fasc` + sidecar `score_*.sc` fillers + final fallback to `POSE_ENERGIES_TABLE` parsed directly from PDBs for any rep that wasn't otherwise covered). |
  | `qc_quarantine` | **0** | Empty; the snapshot is fully consistent. |
  | `build_runs` | 1 | Single locked snapshot 2026-04-27a, qc_status=pass, raw_manifest_hash recomputed against supplemental TSVs. |
  | views: `v_cell_summary`, `v_per_target_mp_means` | — | Pre-built convenience views |

  **Status:** Every fact table at 100% coverage as of 2026-04-28. Loader: `db/scripts/build_db_supplements.py` (idempotent, applies schema migration for `rosetta_energy` + `targets.parent_pdb_id` on first run, wipe-and-reload for energy table on every run). Energy extraction: `red_analysis/scripts/extract_rosetta_energy.py` (patched to handle sidecar `score_*.sc` fill-ins and PDB POSE_ENERGIES_TABLE fallback). Single citation point for the manuscript: snapshot 2026-04-27a.

### Per-target output (Blue pipeline, individual target metrics)
- **ACCRE source:** `accre:/data/p_csb_meiler/agarwm5/red_analysis/metrics/molprobity_blue_<TARGET>.tsv` and `molprobity_green_<TARGET>.tsv` (one per target × per pipeline)
- **Use case:** debugging a single target. Not synced by default (would be 514 small files); pull individually via `ssh accre 'cat /data/p_csb_meiler/agarwm5/red_analysis/metrics/molprobity_blue_1ACB.tsv'`.

### PowerPoint slides (PI-ready)
- **Local:** `briefs/pptx/` (after `make sync-pptx`)
- **ACCRE source:** `accre:/data/p_csb_meiler/agarwm5/red_analysis/BM55_Relaxation_Benchmark_*.pptx`
- **Inventory:** Blue (16 slides), Green (16 slides), combined.

### Crystal PDBs (the BM5.5 originals)
- **ACCRE source:** `accre:/data/p_csb_meiler/agarwm5/protein_pipeline/cleaned/<TARGET>.pdb` and `accre:/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/cleaned/`
- **Repo:** `repos/Protein_Relax_Pipeline/cleaned/` (257 cleaned crystal PDBs)
- **FASTAs:** `repos/Protein_Relax_Pipeline/fasta/` (257 crystal-derived FASTAs)

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
- **Vault:** `Dreamless_Machine/03-Research/Protein/Docs/manuscript/` (subdir to be created when the proteins agent generates fresh v2 drafts; v1 cleared 2026-04-28)
- **Workspace:** `manuscript/` (workspace-side new sections, response-to-reviewer drafts, build artifacts)
- **Style:** Academic register from `~/.claude/VOICE.md` (Bair pattern). No iter logs in frontmatter — use `verifier_history` list.
- **Target:** Proteins: Structure, Function, and Bioinformatics (250-word abstract limit).

### PI briefs and weekly status
- **Workspace:** `briefs/<YYYY-MM-DD>_<topic>.md`
- **Vault PI brief reference:** `pi_review_2026-04-26/PI_REVIEW_BRIEF.md` (DELIVERED, post-meeting note added)
- **Template:** model new briefs on the canonical PI brief structure (1-line narrative → 3 headline findings → method defensibility → open items → file index).

### Vault landmark docs
- `Plan.md` — strategic plan (operator-editable)
- `Status.md` — live status (operator refreshes)
- `Notes.md`, `Ideas.md` — Damien / Ary owned, free-form
- `Protein.md` — hub dashboard with dataview queries
- `Docs/{Architecture, Overview, Decisions, Glossary, References, Changelog, Status}.md` — landmark docs
- `Docs/audits/` — historical audit snapshots (citation audit, data completeness iter48, interim rosetta 2026-04-18)
- `pi_review_2026-04-26/` — PI brief + files index

### Sub-agent (Claude Code)
- **Location:** `~/.claude/agents/proteins.md`
- **Auto-fires:** for any session in `~/proteins-workspace/` or any session about BM5.5 / protein docking / Meiler-lab work.
- **Knows:** dataset (257 + 7 sources + 27 inputs + 6 protocols + 416,340 lock), three findings, repo + ACCRE + vault layouts, figure naming, voice rules.

---

## Repo map (GitHub, public)

| Repo | URL | Role |
|---|---|---|
| `Protein_Relax_Pipeline` | github.com/dreamlessx/Protein_Relax_Pipeline | BLUE pipeline + canonical analysis (`red_analysis/scripts/`) + DB schema/build + manuscript code. **Primary reference repo.** |
| `Protein_Ideal` | github.com/dreamlessx/Protein_Ideal | GREEN pipeline + dataset prep |
| `Protein_Data_Analysis` | github.com/dreamlessx/Protein_Data_Analysis | MolProbity validation toolchain (Phase 1 pilot, 20 proteins) |

All three are clean at upstream HEAD as of 2026-04-28. DB at 100% coverage: `rosetta_metrics` 416,340 + `prerosetta_metrics` 13,364 + `tm_scores` 104,765 + `rosetta_energy` 416,340 + `targets` 257 (full metadata + parent_pdb_id) + `qc_quarantine` 0. Loader: `Protein_Relax_Pipeline/db/scripts/build_db_supplements.py`. Energy extractor: `Protein_Relax_Pipeline/red_analysis/scripts/extract_rosetta_energy.py`. Live SQLite + 5 raw TSVs in the `db-2026-04-27a-supp` GitHub Release.

---

## Number reference (verified against locked DB)

```
257 BM5.5 complexes      = 162 rigid + 60 medium + 35 difficult (incl. 4 non-standard: BAAD, BOYV, BP57, CP57)
605 chains, 122,966 res  = post-strip across 257; 119 homo-multimers deduped to unique seqs; 41 his-tags removed
6,939 input structures   = 257 crystal + 1,285 AF ranked + 1,285 AF unrelaxed + 1,285 Boltz + 257 amber_crystal + 1,285 amber_af + 1,285 amber_boltz
416,340 Rosetta MP rows  = 5 × 77,100 (5-model sources) + 2 × 15,420 (1-model sources). NOT 2× — Blue (208,170 in Protein_Relax_Pipeline) + Green (208,170 in Protein_Ideal), DB unifies.
208,170 per pipeline     = 257 × 810 outputs/target = 257 × 27 × 6 × 5
13,364 pre-Rosetta MP    = the initials' MolProbity, post-Stage-C. 13,344 from combined TSV + 20 Stage C Blue crystal backfill. Now 6,682 Blue + 6,682 Green; symmetric.
12,065 pre-Rosetta TM    = initials vs crystal TM-score
92,700 post-Rosetta TM   = Rosetta outputs vs crystal TM-score
184,352 Rosetta energy   = energy values from Rosetta runs
1,720 filtered at lock   = 663 exact-duplicate + 27 legacy-source + 1,030 other dedup
6,820 (Phase 1 pilot)    = 20 proteins × 341 structures (Protein_Data_Analysis only, pre-full-benchmark)
```

---

## Three findings (manuscript-ready, source-of-truth)

1. **AMBER free lunch.** Cliff's d clashscore = -0.99 at TM Cliff's d = -0.01. 257/257 AF + 256/257 Boltz improved.
2. **Crystal worst MolProbity.** Crystal clashscore 13.85 vs relaxed predictions 0.68-2.82. Frame as idealization artifact.
3. **Rosetta protocol ranking.** dualspace_beta wins integrated MP (0.22, clashscore 0.68) at 0.019-0.025 TM cost. beta_nov16 dominates ref2015 on MP/clash/rama (40-42 of 42 triples), reverses on rotamer outliers (21/42).
