# Comparison Framework

**Last updated:** 2026-04-27 (snapshot 2026-04-27a, pipeline locked at 100.000%)

## Overview

This document describes the BM5.5 benchmarking framework from Blue's perspective. Blue is the primary pipeline housed in this repo (`Protein_Relax_Pipeline`); Green is the matched-parameters re-run housed in companion repo `Protein_Ideal`. The two pipelines share dataset, FASTA strategy, Rosetta flags, and AMBER parameters. They differ on infrastructure (job naming, directory layout, Rosetta version 3.14 vs 3.15, AF AMBER on GPU vs CPU). The Protein_Ideal/COMPARISON.md mirror documents the same diff from Green's perspective.

This document describes the benchmarking framework for comparing protein structure prediction and relaxation methods against experimental crystal structures from the BM5.5 dataset.

## Input Sources (7 buckets)

| # | Input Type | Source | Description |
|---|-----------|--------|-------------|
| 1 | `af_relaxed` | AlphaFold 2.3.2 | Built-in AMBER relaxation (control), 5 ranked models per target |
| 2 | `af_unrelaxed` | AlphaFold 2.3.2 | Raw neural network output, 5 ranked models per target |
| 3 | `boltz` | Boltz-1 v0.4.1 | Diffusion model prediction, 5 samples per target |
| 4 | `amber_af` | Standalone AMBER | AF unrelaxed + standalone AMBER, 5 models per target |
| 5 | `amber_boltz` | Standalone AMBER | Boltz + standalone AMBER, 5 models per target |
| 6 | `crystal` | PDB/RCSB | Experimental X-ray/NMR structures, 1 per target |
| 7 | `amber_crystal` | Standalone AMBER | Standalone AMBER on crystal, 1 model per target |

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
| Input structures per target | 27 (5 AF relaxed + 5 AF unrelaxed + 5 Boltz + 5 AMBER(AF) + 5 AMBER(Boltz) + 1 crystal + 1 AMBER(crystal)) |
| Rosetta protocols | 6 |
| Replicates | 5 |
| **Per target per pipeline** | **810** (27 x 6 x 5) |
| **Per pipeline (257 targets)** | **208,170** |
| **Both pipelines combined** | **416,340** |

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
| Targets changed | 241/257 (94%) |
| Net residue change | -3,956 |
| Targets shortened | 227 |
| Targets lengthened | 14 |
| Additional chains gained | 13 targets |
| Mean change per target | -16.4 residues (-2.8%) |

### Prediction Quality Improvement

| Method | Mean Improvement | Targets Improved |
|--------|------------------|------------------|
| Boltz confidence | +0.026 | 72% |
| AF pLDDT | +1.40 | 73% |
| AF targets with +5 pLDDT | 22 | n/a |

## Quality Metrics (Computed)

### Structural Metrics
- RMSD to crystal (backbone + all-atom)
- TM-score
- GDT-TS
- lDDT

### Geometry Validation
- MolProbity (clashscore, Ramachandran, rotamer outliers)
- PoseBusters: deprecated per paper scope

### Energy Metrics
- Rosetta total score (REU)
- Per-residue energy profiles

## Pipeline Status (snapshot 2026-04-27a)

The pipeline is locked at 100.000%.

| Stage | Blue | Green |
|-------|------|-------|
| Crystal structures | 257/257 | 257/257 |
| FASTAs | 257/257 | 257/257 |
| AlphaFold predictions | 257/257 | 257/257 |
| AF built-in AMBER | 257/257 | 257/257 |
| AF unrelaxed | 257/257 | 257/257 |
| Boltz predictions | 257/257 | 257/257 |
| Standalone AMBER (AF + Boltz + crystal) | 257/257 | 257/257 |
| Rosetta MolProbity rows | 208,170 / 208,170 | 208,170 / 208,170 |
| Combined Rosetta MolProbity rows | 416,340 / 416,340 | (locked) |

## Green Pipeline

Green is the independent verification of the Blue pipeline's protocol. The green pipeline reproduces the exact same Rosetta flags and AMBER parameters to validate reproducibility.

- All SLURM jobs tagged with `green_` prefix
- Output files are `.pdb.gz` compressed
- Blue's Rosetta flags matched exactly

## Summary

Snapshot 2026-04-27a is the locked single citation point. Both Blue and Green hit 100.000% (208,170 / 208,170 each, 416,340 / 416,340 combined). The matched-parameter Green re-run reproduces Blue's three findings: (1) AMBER fixes local geometry without touching global fold (clashscore Cliff's d = -0.99 at TM Cliff's d = -0.01); (2) Crystal carries the worst pre-Rosetta MolProbity (idealization artifact); (3) `dualspace_beta` wins integrated MolProbity at small TM cost. Canonical figures, tables, and pptx live in `red_analysis/`. Phase 1 pilot (20 proteins) is in companion repo `Protein_Data_Analysis`.
