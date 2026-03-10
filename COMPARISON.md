# Comparison Framework

**Last updated:** 2026-03-09

## Overview

This document describes the benchmarking framework for comparing protein structure prediction and relaxation methods against experimental crystal structures from the BM5.5 dataset.

## Input Sources (6 Types)

| # | Input Type | Source | Description |
|---|-----------|--------|-------------|
| 1 | `af_relaxed` | AlphaFold 2.3.2 | Built-in AMBER relaxation (control) |
| 2 | `af_unrelaxed` | AlphaFold 2.3.2 | Raw neural network output |
| 3 | `boltz` | Boltz-1 v0.4.1 | Diffusion model prediction |
| 4 | `amber_af` | Standalone AMBER | AF unrelaxed + standalone AMBER relaxation (test) |
| 5 | `amber_boltz` | Standalone AMBER | Boltz + standalone AMBER relaxation (test) |
| 6 | `crystal` | PDB/RCSB | Experimental X-ray/NMR structures |

## Relaxation Protocols (6 Rosetta + 1 AMBER)

### AMBER (AlphaFold Native)

| Parameter | Value |
|-----------|-------|
| Force field | AMBER ff14SB |
| Energy tolerance | 2.39 kcal/mol |
| Position restraint | 10.0 kcal/mol/A^2 |
| Role | Control (AF built-in) and test (standalone) |

### Rosetta 3.15 Protocols

| Protocol | Space | Scoring Function |
|----------|-------|------------------|
| `cartesian_beta` | Cartesian | beta_nov16 |
| `cartesian_ref15` | Cartesian | REF2015 |
| `dualspace_beta` | Dual (bond geometry) | beta_nov16 |
| `dualspace_ref15` | Dual (bond geometry) | REF2015 |
| `normal_beta` | Torsion | beta_nov16 |
| `normal_ref15` | Torsion | REF2015 |

Each Rosetta protocol uses 5 replicates for statistical robustness.

## Experimental Design

### Key Comparisons

1. **AF built-in AMBER vs standalone AMBER**: Same relaxation method, different implementations
   - AF built-in AMBER = control
   - Standalone AMBER = test
2. **AMBER vs Rosetta**: Different force fields and relaxation strategies
3. **Input quality**: How does input source (AF vs Boltz vs crystal) affect relaxation outcome?
4. **Rosetta protocol comparison**: Cartesian vs dualspace vs normal; beta_nov16 vs REF2015

### Runs per Target

| Component | Count |
|-----------|-------|
| Input types | 6 |
| Rosetta protocols | 6 |
| Replicates | 5 |
| **Total Rosetta runs** | **180** (6 x 6 x 5) |
| **Total per target** | **780** (6 inputs x 6 protocols x 5 reps, max) |

**Note:** Not all input types have 5 models (Boltz has 5 diffusion samples, AF has 5 ranked models), so actual run count may vary.

## FASTA Strategy: Crystal-Derived

All predictions use FASTAs derived directly from crystal PDB coordinates to ensure exact sequence-structure consistency.

### Why Crystal-Derived?

- RCSB canonical sequences differ from resolved crystal constructs
- UniProt full-length sequences include unresolved termini
- Expression tags (His-tags in 41 targets) are artifacts
- Ensures fair comparison: predicted structure covers exactly the same residues as the crystal reference

### Impact

| Metric | Value |
|--------|-------|
| Targets changed | 236/257 (92%) |
| Net residue change | -3,956 |
| Targets shortened | 227 |
| Targets lengthened | 14 |
| Additional chains gained | 13 targets |
| Mean change per target | -16.4 residues (-2.8%) |

### Prediction Quality Improvement

| Method | Mean Improvement | Targets Improved |
|--------|-----------------|------------------|
| Boltz confidence | +0.026 | 72% |
| AF pLDDT | +1.40 | 73% |
| AF targets with +5 pLDDT | 22 | -- |

## Quality Metrics (Planned)

### Structural Metrics
- RMSD to crystal (backbone + all-atom)
- TM-score
- GDT-TS
- lDDT

### Geometry Validation
- MolProbity (clashscore, Ramachandran, rotamer outliers)
- PoseBusters (bond lengths, angles, planarity)

### Energy Metrics
- Rosetta total score (REU)
- Per-residue energy profiles

## Current Progress

| Stage | Status |
|-------|--------|
| Crystal structures | 257/257 cleaned and stripped |
| FASTAs | 257/257 crystal-derived, verified consistent |
| AlphaFold predictions | 257/257 COMPLETE |
| AF built-in AMBER | 257/257 COMPLETE |
| AF unrelaxed | 257/257 COMPLETE |
| Boltz predictions | 257/257 COMPLETE |
| Standalone AMBER | 256/257 (1KTZ finishing, Job 9399617) |
| Rosetta relaxation | IN PROGRESS (~2,512/~200k runs, Jobs 9373165 + 9371774) |
| MolProbity | Pending |
| PoseBusters | Pending |

## Green Pipeline

This is our independent verification ("green") of the Blue pipeline's protocol. The green pipeline reproduces the exact same Rosetta flags and AMBER parameters to validate reproducibility.

- All SLURM jobs tagged with `green_` prefix
- Output files are `.pdb.gz` compressed
- Blue's Rosetta flags matched exactly
