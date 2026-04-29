# BM5.5 Docking-Relaxation Benchmark

A locked, queryable evaluation of structure-prediction and relaxation pipelines on the BM5.5 protein-protein docking dataset. Contains the Blue pipeline implementation, the canonical analysis (`red_analysis/`), the SQLite schema and build code, and the locked snapshot 2026-04-27a published as a GitHub Release.

## What this benchmark answers

Existing structure-prediction benchmarks evaluate accuracy via TM-score, which is insensitive to local geometry. MolProbity reveals stereochemical defects that TM-score does not. The relaxation step exists to fix those defects, but no benchmark had systematically evaluated relaxation across a complete docking dataset with statistical pairing. This work fills that gap.

We ran AlphaFold 2.3.2 (5 ranked + 5 unrelaxed models) and Boltz-1 v0.4.1 (5 single-sequence models) on all 257 BM5.5 complexes. We then exposed every prediction (and the crystal reference) to a relaxation matrix: AMBER as a pre-conditioning step, followed by Rosetta 3.15 under six protocols (cartesian / dualspace / normal × beta_nov16 / ref2015), each replicated five times. Two independent pipelines (Blue, primary; Green, matched-parameters re-run) produce 416,340 Rosetta MolProbity rows, scored against TM-score and total energy at every cell.

## Three findings

**1. AMBER fixes local geometry without touching global fold.** AMBER relaxation reduces clashscore with Cliff's d = -0.99 across all paired comparisons (n = 257 per pipeline) at TM-score Cliff's d = -0.01. AMBER improves MolProbity for 257 of 257 AlphaFold targets and 256 of 257 Boltz targets. Defensible recommendation: AMBER is a near-complete clash-removal step at near-zero accuracy cost.

**2. Crystal structures carry the worst pre-Rosetta MolProbity.** Crystal clashscore averages 13.85; AlphaFold-relaxed averages 2.82; AMBER(Boltz) averages 1.60. The framing is idealization, not failure: predictions come from neural networks trained on energy-minimized targets, so their local geometry is closer to ideal than crystallographic structures whose stereochemistry reflects cryogenic packing constraints. Crystals remain the ground truth for global fold; predictions plus AMBER produce cleaner local geometry.

**3. dualspace_beta wins integrated MolProbity at small TM cost.** Across all 42 (pipeline, source, move-set) triples, beta_nov16 beats ref2015 on MP score, clashscore, and Rama-favored. dualspace_beta produces the lowest aggregate MolProbity score (0.222 on amber_boltz inputs) at a TM-score cost of 0.019 to 0.025 versus the least-perturbing protocols. cartesian_beta is the runner-up with comparable MP and tied TM-score retention. ref2015 retains a niche advantage on rotamer outliers (21 of 42 triples) but loses everywhere else.

## Dataset

Source: [Protein-Protein Docking Benchmark 5.5](https://zlab.wenglab.org/benchmark/) (Vreven 2015; Guest 2021).

| Quantity | Value |
|---|---|
| Total complexes | 257 |
| Rigid-body | 162 |
| Medium | 60 |
| Difficult | 35 |
| Total chains | 605 |
| Total residues | 122,966 |
| Non-standard zlab IDs | 4 (BAAD, BOYV, BP57, CP57) |

The four non-standard IDs map to canonical parent PDBs with explicit chain selections: BAAD = 3AAD_A:B (double bromodomain + ASF1), BOYV = 1OYV_B:I (subtilisin + tomato inhibitor), BP57 = 3P57_AB:P (MEF2A dimer + p300 TAZ2, AB:P selection), CP57 = 3P57_CD:P (same parent, CD:P selection). The DB stores `parent_pdb_id` for each non-standard target.

FASTAs are derived from crystal coordinates, not RCSB canonical sequences. Of 257 targets, 241 differ from RCSB; the net change is -3,956 residues from trimming unresolved termini, removing expression tags (41 targets had His-tags), and dropping construct extensions absent from electron density. Crystal-derived FASTAs improve prediction quality: Boltz mean confidence rises by 0.026 (72% of targets improved); AlphaFold mean pLDDT rises by 1.40 (73% of targets improved, 22 by more than 5 pLDDT).

DNA and RNA chains are excluded; AlphaFold-Multimer accepts protein only, Boltz-1 standard models are protein-focused, and BM5.5 evaluates protein-protein interfaces.

## Relaxation matrix

Each complex feeds 27 input structures into the Rosetta protocol grid:

| Source bucket | Models per target |
|---|---|
| `crystal` | 1 (cleaned crystal) |
| `amber_crystal` | 1 (AMBER-relaxed crystal) |
| `af_relaxed` | 5 (AlphaFold ranked, AlphaFold-internal AMBER) |
| `af_unrelaxed` | 5 (AlphaFold ranked, no relaxation) |
| `amber_af` | 5 (standalone AMBER on `af_unrelaxed`) |
| `boltz` | 5 (Boltz-1 single-sequence) |
| `amber_boltz` | 5 (standalone AMBER on `boltz`) |

Per target: 27 inputs × 6 protocols × 5 replicates = 810 Rosetta runs per pipeline. Per pipeline: 810 × 257 = 208,170. Two pipelines: 416,340 cells in the locked DB.

The six Rosetta 3.15 protocols cross three move sets with two scoring functions:

| Protocol | Move set | Scoring function |
|---|---|---|
| `cartesian_beta` | Cartesian | beta_nov16 |
| `cartesian_ref15` | Cartesian | REF2015 |
| `dualspace_beta` | Dualspace (bond geometry) | beta_nov16 |
| `dualspace_ref15` | Dualspace (bond geometry) | REF2015 |
| `normal_beta` | Torsion | beta_nov16 |
| `normal_ref15` | Torsion | REF2015 |

AMBER relaxation uses AlphaFold's built-in OpenMM with the AMBER ff14SB force field (energy tolerance 2.39 kcal/mol, position restraint 10.0 kcal/mol/Å²) on GPU. Standalone AMBER (the `amber_af` and `amber_boltz` sources) applies the same force field independently to AlphaFold-unrelaxed and Boltz-1 outputs.

## Locked snapshot 2026-04-27a

The DB is the canonical artifact. Single citation point for the manuscript.

| Table | Rows | Notes |
|---|---|---|
| `rosetta_metrics` | 416,340 | Locked. 663 exact-duplicate + 27 legacy-source rows filtered at ingest. |
| `prerosetta_metrics` | 13,364 | Pre-Rosetta MolProbity on the 27 input sources. 13,344 from `combined_molprobity.tsv` plus 20 Stage C Blue crystal rows backfilled from Green (PDBs verified byte-identical, MolProbity deterministic). |
| `tm_scores` | 104,765 | 12,065 pre-Rosetta TM-score against crystal reference, plus 92,700 post-Rosetta TM-score per Rosetta cell. |
| `rosetta_energy` | 416,340 | Total score and per-residue energy per cell. 1:1 with `rosetta_metrics`. Recovered from `relax.fasc` plus per-rep `score_*.sc` sidecars plus `POSE_ENERGIES_TABLE` parsing of PDB outputs as a final fallback. |
| `targets` | 257 | Difficulty (162R/60M/35D), category (zlab AA/AS/EI/ER/ES/OG/OR/OX), n_chains, n_residues, non_standard_flag, parent_pdb_id. |
| `qc_quarantine` | 0 | Clean. |
| `pipelines` / `sources` / `protocols` | 2 / 7 / 6 | Dimension tables. |

Build: `db/scripts/build_db.py` produces the `rosetta_metrics` lock from per-target TSVs. `db/scripts/build_db_supplements.py` is the idempotent companion that loads the supplemental tables, applies the schema migration for `rosetta_energy` and `targets.parent_pdb_id`, and refreshes the build manifest hash. Both run under the same snapshot ID; the manuscript cites `snapshot 2026-04-27a` once.

Release: [`db-2026-04-27a-supp`](https://github.com/dreamlessx/Protein_Relax_Pipeline/releases/tag/db-2026-04-27a-supp) on this repo. Six assets: `bm55.sqlite` (~225 MB, VACUUMed) plus the five raw TSVs (`combined_molprobity.tsv`, `combined_tmscore.tsv`, `combined_rosetta_molprobity.tsv`, `combined_rosetta_tmscore.tsv`, `combined_rosetta_energy.tsv`).

## Repository layout

```
Protein_Relax_Pipeline/
├── cleaned/                    257 cleaned crystal PDBs (homo-multimer dedup applied)
├── fasta/                      257 crystal-derived FASTAs
├── predictions/
│   ├── alphafold/              5 AMBER-relaxed ranked PDBs per target
│   ├── alphafold_unrelaxed/    5 unrelaxed ranked PDBs per target
│   └── boltz/                  5 single-sequence Boltz-1 PDBs per target
├── scripts/                    Pipeline scripts (Blue + shared)
├── db/
│   ├── sql/schema.sql          Schema (10 tables, 2 views)
│   ├── scripts/build_db.py     Locks rosetta_metrics from per-target TSVs
│   ├── scripts/build_db_supplements.py  Loads supplements, applies schema migration
│   └── data/                   bm55_difficulty.csv, bm55_chains.csv (loader inputs)
├── red_analysis/               Canonical analysis: figures, tables, scripts, target list, paper findings
├── DATA_INDEX.md               "I want X. Where do I get it?"
├── PROJECT_STATUS.md           Current state at snapshot lock
└── README.md                   This file
```

## Pipelines

Blue is primary. Green is an independent matched-parameters re-run that lives in [`Protein_Ideal`](https://github.com/dreamlessx/Protein_Ideal). The two pipelines use the same 257 targets, the same 27 input structures per target, the same 6 Rosetta protocols with identical flags, and the same 5-replicate count. Blue runs at `/data/p_csb_meiler/agarwm5/protein_pipeline/` on ACCRE with prefix `blue_`; Green runs at `/data/p_csb_meiler/agarwm5/protein_ideal_test/` with prefix `green_`. The DB unifies both under the same snapshot.

Reproducibility (Blue vs Green agreement): pre-Rosetta TM Pearson r = 0.997 (n = 1,128); pre-Rosetta RMSD r = 0.994; post-Rosetta TM r = 0.999 (n = 60); per-source clashscore r = 0.867 to 0.991; per-source MP score r = 0.941 to 0.984. The independent run reproduces Blue's three findings.

## Computational resources

All predictions and relaxations ran on ACCRE (Vanderbilt University HPC).

| Resource | Specification |
|---|---|
| AlphaFold | NVIDIA RTX A6000, partition `csb_gpu_acc`, 80 GB RAM |
| Boltz-1 | NVIDIA L40S 48 GB, partition `p_meiler_acc` |
| Rosetta | CPU, partition `batch` (`p_csb_meiler`) |
| AMBER (standalone) | GPU OpenMM, on AlphaFold partition |

All SLURM array scripts include `#SBATCH --exclude=cn1340` after a node-specific failure pattern produced 1,614 instant-failure jobs traced to that node.

## Quickstart

```bash
git clone git@github.com:dreamlessx/Protein_Relax_Pipeline.git
cd Protein_Relax_Pipeline

# Pull the locked DB + raw TSVs from the GitHub Release
gh release download db-2026-04-27a-supp --repo dreamlessx/Protein_Relax_Pipeline --dir release_assets/

# Or rebuild end-to-end from the raw TSVs in the release
python db/scripts/build_db.py \
  --raw-root release_assets/ \
  --out      bm55.sqlite \
  --snapshot 2026-04-27a

python db/scripts/build_db_supplements.py \
  --db              bm55.sqlite \
  --raw-root        release_assets/ \
  --difficulty-csv  db/data/bm55_difficulty.csv \
  --chains-csv      db/data/bm55_chains.csv \
  --snapshot-id     2026-04-27a
```

Query examples in `DATA_INDEX.md`. Canonical analysis scripts in `red_analysis/scripts/`. Paper findings summary in `red_analysis/PAPER_FINDINGS.md`.

## Companion repos

| Repo | Role |
|---|---|
| [`Protein_Ideal`](https://github.com/dreamlessx/Protein_Ideal) | Green pipeline (matched-parameters re-run). Independent verification of Blue. |
| [`Protein_Data_Analysis`](https://github.com/dreamlessx/Protein_Data_Analysis) | Phase 1 pilot (20 proteins, 6,820 structures). Established the validation methodology before scaling to BM5.5. Kept for the historical record. |

## References

1. Jumper, J., Evans, R., Pritzel, A. et al. Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583-589 (2021).
2. Wohlwend, J., Corso, G., Passaro, S. et al. Boltz-1: Democratizing Biomolecular Interaction Modeling. *bioRxiv* (2024).
3. Eastman, P., Swails, J., Chodera, J.D. et al. OpenMM 7: Rapid Development of High Performance Algorithms for Molecular Dynamics. *PLOS Comput. Biol.* 13(7): e1005659 (2017).
4. Conway, P., Tyka, M.D., DiMaio, F., Konerding, D.E., Baker, D. Relaxation of backbone bond geometry improves protein energy landscape modeling. *Protein Sci.* 23, 47-55 (2014).
5. Park, H., Bradley, P., Greisen, P. Jr et al. Simultaneous Optimization of Biomolecular Energy Functions on Features from Small Molecules and Macromolecules. *J. Chem. Theory Comput.* 12, 6201-6212 (2016) (beta_nov16).
6. Alford, R.F., Leaver-Fay, A., Jeliazkov, J.R. et al. The Rosetta All-Atom Energy Function for Macromolecular Modeling and Design. *J. Chem. Theory Comput.* 13, 3031-3048 (2017) (REF2015).
7. Williams, C.J., Headd, J.J., Moriarty, N.W. et al. MolProbity: More and better reference data for improved all-atom structure validation. *Protein Sci.* 27, 293-315 (2018).
8. Zhang, Y., Skolnick, J. TM-align: a protein structure alignment algorithm based on the TM-score. *Nucleic Acids Res.* 33, 2302-2309 (2005).
9. Vreven, T., Moal, I.H., Vangone, A. et al. Updates to the Integrated Protein-Protein Interaction Benchmarks. *J. Mol. Biol.* 427, 3031-3041 (2015) (BM5.0).
10. Guest, J.D., Vreven, T., Zhou, J. et al. An expanded benchmark for antibody-antigen docking and affinity prediction reveals insights into antibody recognition determinants. *J. Mol. Biol.* 433, 166983 (2021) (BM5.5).

## License

MIT.

---

*Snapshot 2026-04-27a, locked at 100.000% on 2026-04-27. Every fact table at 100% coverage; qc_quarantine clean. Last verified 2026-04-28.*
