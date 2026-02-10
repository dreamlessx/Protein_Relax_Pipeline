# Protein Structure Prediction & Relaxation Pipeline

A comprehensive benchmarking framework for evaluating AI-based protein structure prediction methods against the Protein-Protein Docking Benchmark 5.5 (BM5.5) dataset.

## Overview

This repository contains structure predictions and relaxation protocols for protein-protein complexes from the BM5.5 benchmark. We systematically compare AlphaFold 2.3.2 and Boltz-1 predictions against experimental crystal structures, with subsequent relaxation across multiple scoring functions.

**Key Features:**
- Full BM5.5 coverage (257 complexes)
- Both unrelaxed and AMBER-relaxed AlphaFold outputs
- 7 relaxation protocols (1 AMBER + 6 Rosetta)
- Automatic database fallback for antibody sequences

## Dataset

**Source:** [Protein-Protein Docking Benchmark 5.5](https://zlab.wenglab.org/benchmark/)

| Category | Count |
|----------|-------|
| Total complexes | 257 |
| Rigid-body | 162 |
| Medium difficulty | 60 |
| Difficult | 35 |

### Non-Standard BM5.5 Entries

BM5.5 includes 4 non-standard entries representing alternate chain combinations:

| ID | Parent PDB | Description |
|----|------------|-------------|
| BAAD | 3AAD | Alternate chains A:B (Double bromodomain + ASF1) |
| BOYV | 1OYV | Chains B:I (Subtilisin + Tomato inhibitor) |
| BP57 | 3P57 | Chains AB:P (MEF2A dimer + p300 TAZ2) |
| CP57 | 3P57 | Chains CD:P (MEF2A dimer + p300 TAZ2) |

## Structure Prediction Methods

### AlphaFold 2.3.2

| Parameter | Value |
|-----------|-------|
| Database | `full_dbs` (BFD + UniRef30 via HHblits) with auto-fallback to `reduced_dbs` |
| Model preset | Auto-detect (monomer/multimer) |
| Models per target | 5 ranked models |
| Relaxation | All 5 models AMBER-relaxed (`--models_to_relax=all`) |
| Output | Both relaxed and unrelaxed saved |
| Memory | 80 GB RAM |
| GPU | NVIDIA RTX A6000 |

### Boltz-1 v0.4.1

| Parameter | Value |
|-----------|-------|
| MSA | MSA server |
| Recycling steps | 10 |
| Sampling steps | 200 |
| Diffusion samples | 5 per target |
| Output | Unrelaxed (native Boltz) |

## Directory Structure

```
data/
├── {PDB_ID}/
│   ├── af_out/                    # AlphaFold AMBER-relaxed
│   │   ├── ranked_0.pdb
│   │   ├── ranked_1.pdb
│   │   ├── ranked_2.pdb
│   │   ├── ranked_3.pdb
│   │   └── ranked_4.pdb
│   ├── af_out_unrelaxed/          # AlphaFold unrelaxed
│   │   ├── ranked_0.pdb
│   │   ├── ranked_1.pdb
│   │   ├── ranked_2.pdb
│   │   ├── ranked_3.pdb
│   │   └── ranked_4.pdb
│   ├── boltz_out/                 # Boltz-1 predictions (unrelaxed)
│   │   ├── boltz_input_model_0.pdb
│   │   ├── boltz_input_model_1.pdb
│   │   ├── boltz_input_model_2.pdb
│   │   ├── boltz_input_model_3.pdb
│   │   └── boltz_input_model_4.pdb
│   ├── db_preset_used.txt         # Records full_dbs or reduced_dbs
│   ├── boltz_input.fasta          # Boltz-format FASTA
│   └── sequence.fasta             # Standard FASTA
└── ...

scripts/
├── prediction/
│   ├── af_full.slurm              # AlphaFold with full_dbs + fallback
│   ├── alphafold_array.slurm      # AlphaFold array job
│   ├── alphafold_single.slurm     # AlphaFold single target
│   ├── boltz_array.slurm          # Boltz-1 array job
│   └── boltz_single.slurm         # Boltz-1 single target
├── relaxation/                    # Rosetta relaxation scripts
├── analysis/                      # Analysis scripts
├── data_preparation/              # FASTA preparation
└── validation/                    # Validation scripts
```

## Relaxation Protocols

### AMBER Relaxation (AlphaFold Native)

AlphaFold's built-in AMBER relaxation using OpenMM:

| Parameter | Value |
|-----------|-------|
| Force field | AMBER ff14SB |
| Energy tolerance | 2.39 kcal/mol |
| Position restraint | 10.0 kcal/mol/Å² |
| Acceleration | GPU (`--use_gpu_relax`) |

### Rosetta Relaxation Protocols

Six Rosetta relaxation protocols with 5 replicates each:

| Protocol | Space | Scoring Function |
|----------|-------|------------------|
| `cartesian_beta` | Cartesian | beta_nov16 |
| `cartesian_ref15` | Cartesian | REF2015 |
| `dualspace_beta` | Dual (bond geometry) | beta_nov16 |
| `dualspace_ref15` | Dual (bond geometry) | REF2015 |
| `normal_beta` | Torsion | beta_nov16 |
| `normal_ref15` | Torsion | REF2015 |

### Summary: 7 Relaxation Types

| # | Type | Method | Applied To |
|---|------|--------|------------|
| 1 | AMBER (native) | AlphaFold/OpenMM | AlphaFold predictions |
| 2 | cartesian_beta | Rosetta | All predictions |
| 3 | cartesian_ref15 | Rosetta | All predictions |
| 4 | dualspace_beta | Rosetta | All predictions |
| 5 | dualspace_ref15 | Rosetta | All predictions |
| 6 | normal_beta | Rosetta | All predictions |
| 7 | normal_ref15 | Rosetta | All predictions |

## Computational Resources

All predictions generated on ACCRE (Vanderbilt University HPC).

| Resource | Specification |
|----------|---------------|
| Partition | batch_gpu (csb_gpu_acc) |
| GPU | NVIDIA RTX A6000 / L40S |
| Memory | 80 GB per job |
| Time limit | 72 hours |

## Technical Notes

### Database Fallback Strategy

The pipeline attempts `full_dbs` first for maximum accuracy, then automatically falls back to `reduced_dbs` if HHblits fails:

```
full_dbs (HHblits + BFD) → fails on antibodies → reduced_dbs (MMseqs2)
```

**Why fallback?** HHblits has a hard-coded limit (32763 residues) that antibody/immunoglobulin sequences exceed when matching titin-like proteins in BFD.

**Impact:** Minimal accuracy loss (~0.5-1% TM-score) since antibodies have extensive UniRef90 coverage.

### Boltz-1 FASTA Format

Boltz-1 requires specific header format:
```
>A|protein|chain A
MKTAYIAKQRQISFVKSH...
>B|protein|chain B
DIVLTQSPASLAVSLGQR...
```

Simple headers like `>A` fail with "Invalid record id" error.

### Storage Management

Intermediate files cleaned after prediction:
- `*.sto` (MSA) - 50-500 MB each
- `features.pkl` - 100-500 MB
- `result_model_*.pkl` - 200-800 MB each

Only PDB outputs and FASTA files retained.

## Run Log

| Date | Event | Details |
|------|-------|---------|
| 2026-02-07 | Initial batch | AF2 with reduced_dbs |
| 2026-02-08 | Cleanup | Removed 15 non-BM5.5 targets |
| 2026-02-08 | Added | 4 non-standard entries (BAAD, BOYV, BP57, CP57) |
| 2026-02-09 | Full re-run | All 257 targets, full_dbs + fallback, both relaxed/unrelaxed |

### Current Progress

| Method | Status | Details |
|--------|--------|---------|
| AlphaFold | Running | Job 8849933, 257 targets, full_dbs + fallback |
| Boltz-1 | Pending | After AlphaFold completion |
| Rosetta relax | Pending | After predictions complete |

## References

1. Jumper, J. et al. Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583–589 (2021).
2. Wohlwend, J. et al. Boltz-1: Democratizing Biomolecular Interaction Modeling. *bioRxiv* (2024).
3. Vreven, T. et al. Updates to the Integrated Protein–Protein Interaction Benchmarks. *J. Mol. Biol.* 427, 3031–3041 (2015).

## License

MIT License

---
*Last updated: 2026-02-09*
