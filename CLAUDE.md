# Protein Relax Pipeline: BM5.5 Prediction + Relaxation

## What This Is
End-to-end benchmarking pipeline: predict 257 protein-protein complexes with
AlphaFold 2.3.2 and Boltz-1 v0.4.1, then relax with AMBER and 6 Rosetta protocols
(5 replicates each). 416,340 total Rosetta MolProbity rows across both pipelines
(208,170 per pipeline).

## Pipeline Status (snapshot 2026-04-27a, 100.000% locked)
- AlphaFold predictions: 257/257 DONE
- Boltz-1 predictions: 257/257 DONE
- Standalone AMBER: 257/257 DONE
- Rosetta relaxation: 416,340 / 416,340 rows complete (208,170 per pipeline)
- 0 gap cells, 0 missing rows, 0 NaN
- 663 exact dups + 27 legacy-source rows filtered at ingest

## Tools Required
| Tool | Version | Purpose |
|------|---------|---------|
| AlphaFold | 2.3.2 | Structure prediction |
| Boltz-1 | 0.4.1 | Diffusion-based prediction |
| Rosetta | 3.15 | 6 relaxation protocols |
| AMBER | ff14SB | GPU-accelerated relaxation |
| Phenix/Reduce | Latest | Hydrogen placement |
| MolProbity | Latest | Validation scoring |

## ACCRE SLURM Configs
- AlphaFold: batch_gpu, A6000, 80GB RAM, 72h, array 1-257 (10 concurrent)
- Boltz: p_meiler_acc, L40S, 80GB RAM, 2 days, array 1-257 (10 concurrent)
- Rosetta: batch CPU, 8GB RAM, 7d, array 1-257 (50 concurrent), 27 inputs x 6 protocols x 5 reps = 810 runs per target
- All SLURM scripts include `#SBATCH --exclude=cn1340` (1,614 failures traced to that node)

## Key Scripts
- derive_all_from_crystal.py: extract crystal-derived FASTAs
- strip_crystals.py: remove homo-multimer duplicates
- af_consistency_rerun.slurm: AlphaFold batch prediction
- boltz_consistency_rerun.slurm: Boltz-1 batch prediction
- amber_relax.py + amber_relax.slurm: standalone AMBER relaxation
- amber_relax_crystal_v5.py: chain-split AMBER for crystals (1ACB, 1ATN fix)

## Data Layout
cleaned/          257 crystal PDBs (unique chains)
fasta/            257 crystal-derived sequences
predictions/
  alphafold/      AMBER-relaxed AF2 models
  alphafold_unrelaxed/
  boltz/          Boltz-1 predictions
db/               SQLite schema, build script, queries
red_analysis/     metrics, figures, scripts, tables, presentations

## Rosetta Protocols
6 variants: cartesian / dualspace / normal x beta_nov16 / ref2015.
5 replicates each for statistical robustness.
Input types per pipeline (27 total): 5 AF relaxed + 5 AF unrelaxed + 5 Boltz + 5 AMBER(AF) + 5 AMBER(Boltz) + 1 crystal + 1 AMBER(crystal).

## Available Slash Commands
- /validate: run MolProbity on completed structures
- /status: check ecosystem health and progress
- /submit: submit SLURM jobs
- /compare: generate AF2 vs Boltz1 comparison tables
- /report: generate publication-ready figures and tables
