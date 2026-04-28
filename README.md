# Protein Structure Prediction & Relaxation Pipeline

A comprehensive benchmarking framework for evaluating AI-based protein structure prediction methods against the Protein-Protein Docking Benchmark 5.5 (BM5.5) dataset.

## Overview

This repository contains structure predictions and relaxation protocols for protein-protein complexes from the BM5.5 benchmark. We systematically compare AlphaFold 2.3.2 and Boltz-1 predictions against experimental crystal structures, with subsequent relaxation across multiple scoring functions.

**Key Features:**
- Full BM5.5 coverage (257 complexes)
- Crystal-derived FASTAs: all sequences extracted directly from crystal PDB coordinates
- Perfect consistency: Crystal == AF FASTA == Boltz FASTA for all 257 targets
- No duplicate chains in any FASTA or crystal PDB
- Both unrelaxed and AMBER-relaxed AlphaFold outputs
- 7 relaxation protocols (1 AMBER + 6 Rosetta)

## Dataset

**Source:** [Protein-Protein Docking Benchmark 5.5](https://zlab.wenglab.org/benchmark/)

| Category | Count |
|----------|-------|
| Total complexes | 257 |
| Rigid-body | 162 |
| Medium difficulty | 60 |
| Difficult | 35 |

**Total chains:** 605 across all targets
**Total residues:** 122,966

### Non-Standard BM5.5 Entries

BM5.5 includes 4 non-standard entries representing alternate chain combinations:

| ID | Parent PDB | Description |
|----|------------|-------------|
| BAAD | 3AAD | Alternate chains A:B (Double bromodomain + ASF1) |
| BOYV | 1OYV | Chains B:I (Subtilisin + Tomato inhibitor) |
| BP57 | 3P57 | Chains AB:P (MEF2A dimer + p300 TAZ2) |
| CP57 | 3P57 | Chains CD:P (MEF2A dimer + p300 TAZ2) |

### DNA/RNA Exclusion Policy

DNA/RNA chains present in crystal structures are excluded from predictions:
- AlphaFold Multimer only accepts protein sequences
- Boltz-1 standard models are protein-focused
- BM5.5 evaluates protein-protein interfaces only

## FASTA Derivation from Crystal Structures

All 257 FASTAs are derived directly from crystal PDB coordinates (ATOM records, CA atoms), ensuring exact sequence-structure consistency for benchmarking.

**Why crystal-derived?** RCSB canonical sequences often differ from what is resolved in the crystal:
- Terminal extensions not present in electron density
- Expression tags (His-tags found in 41 targets)
- Missing chains or extra entities in RCSB metadata
- Construct variants vs canonical UniProt sequences

**Impact of crystal derivation (241/257 targets changed):**
- Net change: -3,956 residues (crystal constructs shorter than RCSB canonical)
- 227 targets got shorter (trimmed unresolved termini)
- 14 targets got longer (gained chains missing from RCSB FASTA)
- 13 targets gained additional chains from crystal
- Mean change: -16.4 residues per target (-2.8%)

**Prediction quality improved:**
- Boltz: +0.026 mean confidence score (72% of targets improved)
- AF: +1.40 mean pLDDT (73% of targets improved)
- 22 AF targets improved by +5 pLDDT or more

### Crystal PDB Stripping

36 crystal PDBs were stripped of homo-multimer duplicate chains to match FASTAs. Only unique chains are retained (e.g., 1EXB: 8 chains reduced to 2).

### Consistency Verification

All 257 targets verified across 10 checks:
1. Crystal chain sequences == AF FASTA sequences (exact, bidirectional)
2. AF FASTA == Boltz FASTA (same sequences, same order, same count)
3. No duplicate sequences in any FASTA
4. No duplicate chains in any crystal PDB
5. Crystal chain count == FASTA entry count
6. Boltz headers have correct `|PROTEIN|` format
7. No empty sequences
8. No non-standard amino acids
9. No His-tags or expression artifacts
10. No extra non-protein chains

## Repository Structure

```
Protein_Relax_Pipeline/
├── cleaned/                    # Crystal structures (257 PDBs, stripped to unique chains)
├── fasta/                      # Crystal-derived FASTA sequences (257 files)
├── predictions/
│   ├── alphafold/              # AF2 AMBER-relaxed ranked PDBs
│   │   └── {PDB_ID}/ranked_0..4.pdb
│   ├── alphafold_unrelaxed/    # AF2 unrelaxed ranked PDBs
│   │   └── {PDB_ID}/ranked_0..4.pdb
│   └── boltz/                  # Boltz-1 predictions
│       └── {PDB_ID}/boltz_input_model_0..4.pdb
├── scripts/                    # SLURM and processing scripts
├── db/                         # SQLite schema, build script, queries
├── red_analysis/               # Metrics, figures, scripts, tables, presentations
├── AMBER_FIX_LOG.txt           # Log of AMBER-related FASTA fixes
├── consistency_fix_log.txt     # Log of chain composition fixes
├── histag_fix_log.txt          # Log of His-tag removals
└── test_subset/                # Test data (20 targets)
```

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
| Output format | PDB |
| GPU | NVIDIA L40S 48GB |

## Relaxation Protocols

### AMBER Relaxation (AlphaFold Native)

AlphaFold's built-in AMBER relaxation using OpenMM:

| Parameter | Value |
|-----------|-------|
| Force field | AMBER ff14SB |
| Energy tolerance | 2.39 kcal/mol |
| Position restraint | 10.0 kcal/mol/A^2 |
| Acceleration | GPU (`--use_gpu_relax`) |

### Rosetta 3.15 Relaxation

Six relaxation protocols with 5 replicates each:

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
| 2-7 | Rosetta protocols | Rosetta 3.15 | All 7 input source buckets (af_relaxed, af_unrelaxed, boltz, amber_af, amber_boltz, crystal, amber_crystal) |

Each target produces 810 Rosetta runs per pipeline: 27 inputs (5 AF relaxed + 5 AF unrelaxed + 5 Boltz + 5 AMBER(AF) + 5 AMBER(Boltz) + 1 crystal + 1 AMBER(crystal)) x 6 protocols x 5 replicates. 257 targets per pipeline = 208,170. Both pipelines combined = 416,340.

## Computational Resources

All predictions generated on ACCRE (Vanderbilt University HPC).

| Resource | Specification |
|----------|---------------|
| AF GPU | NVIDIA RTX A6000 (`csb_gpu_acc`) |
| Boltz GPU | NVIDIA L40S 48GB (`p_meiler_acc`) |
| Rosetta | CPU-only, `batch` partition (`p_csb_meiler`) |
| Memory | 80 GB per prediction job |

## Technical Notes

### Database Fallback Strategy

```
full_dbs (HHblits + BFD) -> fails on antibodies -> reduced_dbs (MMseqs2)
```

HHblits has a hard-coded limit (32763 residues) that some antibody/immunoglobulin sequences exceed. Minimal accuracy loss (~0.5-1% TM-score) with fallback.

### Boltz-1 FASTA Format

Boltz-1 requires specific header format:
```
>A|PROTEIN|
MKTAYIAKQRQISFVKSH...
>B|PROTEIN|
DIVLTQSPASLAVSLGQR...
```

### Storage Management

Intermediate files cleaned after prediction:
- `*.sto` (MSA) - 50-500 MB each
- `features.pkl` - 100-500 MB
- `result_model_*.pkl` - 200-800 MB each

Only PDB outputs and FASTA files retained.

## Processing Scripts

| Script | Purpose |
|--------|---------|
| `derive_all_from_crystal.py` | Derive all 257 FASTAs from crystal PDB coordinates |
| `strip_crystals.py` | Strip crystal PDBs to match FASTA (remove homo-multimer duplicates) |
| `af_consistency_rerun.slurm` | AlphaFold prediction for crystal-derived FASTAs |
| `boltz_consistency_rerun.slurm` | Boltz prediction for crystal-derived FASTAs |
| `amber_relax.py` + `amber_relax.slurm` | Standalone AMBER relaxation (AF + Boltz models) |
| `blue_fill_all.slurm` / `green_fill_all.slurm` | Rosetta relaxation: 27 inputs x 6 protocols x 5 reps |
| `blue_amber_crystal_fill.slurm` / `green_amber_crystal_fill.slurm` | AMBER-crystal Rosetta fills with v5 chain-split preprocessing |
| `blue_scan_gaps.sh` / `green_scan_gaps.sh` | Output completeness scanners |

## Pipeline Status (snapshot 2026-04-27a)

The pipeline is locked at 100.000% completion as of 2026-04-27.

| Method | Status | Details |
|--------|--------|---------|
| Input consistency | 257/257 | Crystal == AF FASTA == Boltz FASTA |
| AlphaFold (built-in AMBER) | 257/257 | COMPLETE |
| AlphaFold (unrelaxed) | 257/257 | COMPLETE |
| Boltz-1 | 257/257 | COMPLETE |
| Standalone AMBER (AF) | 257/257 | COMPLETE |
| Standalone AMBER (Boltz) | 257/257 | COMPLETE |
| Standalone AMBER (crystal) | 257/257 | COMPLETE (1ACB + 1ATN resolved via v5 chain-split preprocessing) |
| Rosetta MolProbity rows (Blue) | 208,170 / 208,170 | COMPLETE |
| Rosetta MolProbity rows (Green) | 208,170 / 208,170 | COMPLETE |
| Rosetta MolProbity rows (combined) | 416,340 / 416,340 | 100.000% locked |
| Gap cells | 0 | 0 NaN, 0 missing rows |
| Filtered at ingest | 663 exact dups + 27 legacy-source rows | (does not affect locked count) |

Pipeline definitions:
- Blue pipeline: primary. Job prefix: `blue_`
- Green pipeline: independent verification of Blue. Job prefix: `green_`
- Red analysis: metrics and figures. Job prefix: `red_`

Per-target arithmetic: 27 input structures x 6 protocols x 5 replicates = 810 Rosetta runs per target per pipeline.

Source buckets (7): af_relaxed, af_unrelaxed, amber_af, amber_boltz, amber_crystal, boltz, crystal.

Protocols (6): cartesian_beta, cartesian_ref15, dualspace_beta, dualspace_ref15, normal_beta, normal_ref15.

All SLURM array scripts include `#SBATCH --exclude=cn1340` after 1,614 instant-failure jobs were traced to that node.

## References

1. Jumper, J. et al. Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583-589 (2021).
2. Wohlwend, J. et al. Boltz-1: Democratizing Biomolecular Interaction Modeling. *bioRxiv* (2024).
3. Vreven, T. et al. (2015). Updates to the Integrated Protein-Protein Interaction Benchmarks. *J. Mol. Biol.* 427:3031-3041 (BM5.0).
4. Guest, J.D. et al. (2021). An expanded benchmark for antibody-antigen docking and affinity prediction reveals insights into antibody recognition determinants. *J. Mol. Biol.* 433:166983 (BM5.5).

## License

MIT License

---
*Last updated: 2026-04-27 (snapshot 2026-04-27a, pipeline locked at 100.000%)*
