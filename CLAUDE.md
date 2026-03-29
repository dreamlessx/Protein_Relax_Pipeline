# Protein Relax Pipeline — BM5.5 Prediction + Relaxation

## What This Is
End-to-end benchmarking pipeline: predict 257 protein-protein complexes with
AlphaFold 2.3.2 and Boltz-1 v0.4.1, then relax with AMBER and 6 Rosetta protocols
(5 replicates each). ~200,000 total relaxation runs.

## Pipeline Status
- AlphaFold predictions: 257/257 DONE
- Boltz-1 predictions: 257/257 DONE
- Standalone AMBER: 257/257 DONE
- Rosetta relaxation: IN PROGRESS (~14,600/~200,000 runs)

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
- AlphaFold: batch_gpu, A6000, 80GB RAM, 72h, array 1-236 (10 concurrent)
- Boltz: p_meiler_acc, L40S, 80GB RAM, 2 days, array 1-236 (10 concurrent)
- Rosetta: batch CPU, 8GB RAM, 72h, array 1-257 (50 concurrent)

## Key Scripts
- derive_all_from_crystal.py — Extract crystal-derived FASTAs
- strip_crystals.py — Remove homo-multimer duplicates
- af_consistency_rerun.slurm — AlphaFold batch prediction
- boltz_consistency_rerun.slurm — Boltz-1 batch prediction
- rosetta_relax.slurm — Rosetta relaxation (6 protocols x 5 replicates)

## Data Layout
cleaned/          — 257 crystal PDBs (unique chains)
fasta/            — 257 crystal-derived sequences
predictions/
  alphafold/      — AMBER-relaxed AF2 models
  alphafold_unrelaxed/
  boltz/          — Boltz-1 predictions

## Rosetta Protocols
6 variants: Cartesian/dualspace/normal x beta_nov16/ref2015
5 replicates each for statistical robustness
Input types: AF relaxed, AF unrelaxed, Boltz, AMBER variants, crystal

## Available Slash Commands
- /validate — Run MolProbity on completed structures
- /status — Check Rosetta relaxation progress
- /submit — Submit SLURM jobs for remaining runs
- /compare — Generate AF2 vs Boltz1 comparison tables
- /report — Generate publication-ready figures and tables
