# Protein Structure Prediction & Relaxation Pipeline

A comprehensive benchmarking framework for evaluating AI-based protein structure prediction methods against the Protein-Protein Docking Benchmark 5.5 (BM5.5) dataset.

## Overview

This repository contains structure predictions and relaxation protocols for protein-protein complexes from the BM5.5 benchmark. We systematically compare AlphaFold 2.3.2 and Boltz-1 predictions against experimental crystal structures, with subsequent relaxation across multiple scoring functions.

**Key Features:**
- Full BM5.5 coverage (257 complexes)
- Both unrelaxed and AMBER-relaxed AlphaFold outputs
- 7 relaxation protocols (1 AMBER + 6 Rosetta)
- Deduplicated FASTAs (unique chains only, consistent across AF and Boltz)
- PI-cleaned crystal structures for all 257 targets
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

### DNA/RNA Exclusion Policy

DNA/RNA chains present in crystal structures are excluded from predictions:
- AlphaFold Multimer only accepts protein sequences
- Boltz-1 standard models are protein-focused
- BM5.5 evaluates protein-protein interfaces only

Targets with DNA-mediated interfaces (e.g., 3P57) may have reduced DockQ scores due to missing nucleic acid mediators.

## Repository Structure

```
Protein_Relax_Pipeline/
├── cleaned/                    # PI-cleaned crystal structures (257 PDBs)
├── fasta/                      # Deduplicated FASTA sequences (257 files)
├── predictions/
│   ├── alphafold/              # AF2 AMBER-relaxed ranked PDBs (257 targets)
│   │   └── {PDB_ID}/ranked_0..4.pdb
│   ├── alphafold_unrelaxed/    # AF2 unrelaxed ranked PDBs (257 targets)
│   │   └── {PDB_ID}/ranked_0..4.pdb
│   └── boltz/                  # Boltz-1 predictions (122 correct, 135 re-running)
│       └── {PDB_ID}/boltz_input_model_0..4.pdb
├── scripts/                    # All SLURM scripts
├── AMBER_FIX_LOG.txt           # Log of FASTA modifications for AMBER fixes
└── test_subset/                # Test data (20 targets)
```

### Input Data Folders

| Folder | Files | Description |
|--------|-------|-------------|
| `cleaned/` | 257 PDBs | PI-cleaned crystal structures (Rosetta clean_pdb.py) |
| `fasta/` | 257 FASTAs | Deduplicated sequences from RCSB (unique chains only) |

### FASTA Deduplication

All FASTAs use **unique chains only** (deduplicated). This ensures AlphaFold and Boltz-1 predict the same complex:

- RCSB FASTAs list unique protein entities (e.g., "Chains A, D" = one entry for both copies)
- Both AF and Boltz receive identical sequences per target
- Homo-multimeric copies are NOT duplicated
- Boltz FASTAs use `>A|PROTEIN|` header format, regenerated from `sequence.fasta`

**Impact of deduplication:** 7 targets previously classified as "OOM" were only OOM due to duplicated chains. After dedup, all fit on L40S 48GB.

| Target | Before (duplicated) | After (deduplicated) |
|--------|--------------------|--------------------|
| 1DE4 | 9 chains, 3042 res | 3 chains, 1014 res |
| 1GXD | 4 chains, 1650 res | 2 chains, 825 res |
| 1K5D | 12 chains, 3212 res | 3 chains, 803 res |
| 1N2C | 8 chains, 3182 res | 3 chains, 1302 res |
| 1WDW | 12 chains, 3798 res | 2 chains, 633 res |
| 1ZM4 | 6 chains, 3147 res | 2 chains, 1049 res |
| 6EY6 | 16 chains, 3624 res | 2 chains, 453 res |

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
| GPU | NVIDIA L40S (primary), H100 80GB (for large targets) |

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
| Partition | batch_gpu (p_meiler_acc) |
| GPU | NVIDIA L40S 48GB / H100 80GB |
| Memory | 80-256 GB per job |
| Rosetta | CPU-only, batch partition |

### Rosetta Binary Path

```
/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/source/bin/relax.linuxgccrelease
```
Database: `/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/database`

## Technical Notes

### Database Fallback Strategy

The pipeline attempts `full_dbs` first for maximum accuracy, then automatically falls back to `reduced_dbs` if HHblits fails:

```
full_dbs (HHblits + BFD) → fails on antibodies → reduced_dbs (MMseqs2)
```

**Why fallback?** HHblits has a hard-coded limit (32763 residues) that antibody/immunoglobulin sequences exceed when matching titin-like proteins in BFD.

**Impact:** Minimal accuracy loss (~0.5-1% TM-score) since antibodies have extensive UniRef90 coverage.

### FASTA Verification

All input FASTAs verified against RCSB and crystal structures:
- 10 targets regenerated from RCSB to fix missing/incomplete chains
- All sequences match authoritative source (chain order may differ)
- DNA/RNA chains excluded per policy
- 135 Boltz FASTAs regenerated to remove homo-multimer chain duplication

### AMBER Failure Root Cause

7 original AMBER failures were caused by **non-standard amino acid codes** (X = unknown, Z = ambiguous Glu/Gln) in terminal positions of FASTA sequences.

**Fix:** Trimmed non-standard terminal residues from input FASTAs. Changes logged in `AMBER_FIX_LOG.txt`.

| Target | Chain | Residues Removed | Code | New Length |
|--------|-------|------------------|------|------------|
| 1ATN | 0 (N-term) | 1 `X` | Unknown AA | 372 |
| 1DFJ | 1 (N-term) | 1 `X` | Unknown AA | 456 |
| 1FC2 | 0 (C-term) | 3 `XXK` | Unknown + Lys | 55 |
| 1WEJ | 2 (N-term) | 1 `X` | Unknown AA | 104 |
| 2BTF | 0,1 (N-term) | 1 `X` each | Unknown AA | 374, 139 |
| 4CPA | 1 (N-term) | 2 `ZZ` | Ambiguous Glu/Gln | 36 |
| 5JMO | 2 (both terms) | 2 `X` | Unknown AA | 4 |

### Superseded PDBs

| Original | Superseded By | Resolution |
|----------|---------------|------------|
| 1A2K | 5BXQ | FASTA from 5BXQ, re-predicted |
| 3RVW | 5VPG | FASTA from 5VPG, re-predicted |

### Boltz-1 FASTA Format

Boltz-1 requires specific header format with deduplicated unique chains:
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

## Run Log

| Date | Event | Details |
|------|-------|---------|
| 2026-02-07 | Initial batch | AF2 with reduced_dbs |
| 2026-02-08 | Cleanup | Removed 15 non-BM5.5 targets |
| 2026-02-08 | Added | 4 non-standard entries (BAAD, BOYV, BP57, CP57) |
| 2026-02-09 | Full re-run | All 257 targets, full_dbs + fallback, both relaxed/unrelaxed |
| 2026-02-10 | Input verification | Cross-checked FASTAs, fixed 10 targets |
| 2026-02-10 | FASTA regeneration | Regenerated 10 FASTAs from RCSB, all match authoritative |
| 2026-02-11 | AF complete | 257/257, 7 AMBER failures with unrelaxed saved |
| 2026-02-15 | Boltz batch 1 | 233/257 (24 missing) |
| 2026-02-15 | Rosetta attempt 1 | Job 9011797 — FAILED (wrong binary path) |
| 2026-02-20 | AMBER root cause | X/Z non-standard residues, FASTAs trimmed |
| 2026-02-20 | AMBER re-run | 5/7 fixed, 2 pending (1FC2, 5JMO) |
| 2026-02-22 | Boltz fixes | 12 header/FASTA fixes resubmitted (Job 9195326) |
| 2026-02-22 | Rosetta resubmit | Job 9195328 — correct binary (Rosetta 3.15) |
| 2026-02-22 | GitHub backup | All predictions synced to repo |
| 2026-03-02 | Boltz H100 complete | 5 targets succeeded on H100 80GB |
| 2026-03-02 | 1FC2 + 5JMO fixed | All 7 AMBER failures now relaxed |
| 2026-03-02 | Rosetta progress | 70 done, 17 resumed (time limit), 120 pending |
| 2026-03-03 | FASTA deduplication | Found 135 Boltz FASTAs had duplicated homo-multimer chains |
| 2026-03-03 | OOM targets recovered | All 7 "OOM" targets now feasible after dedup (max 1302 res) |
| 2026-03-03 | Boltz re-run | Job 9304637: 135 targets with deduplicated FASTAs |
| 2026-03-03 | Full 257 restored | AF 257/257 complete, Boltz 122/257 correct + 135 re-running |
| 2026-03-03 | GitHub sync | Crystal structures, FASTAs, AF predictions all pushed |

### Current Progress

**Active benchmark: 257 targets** (full BM5.5)

| Method | Status | Details |
|--------|--------|---------|
| AlphaFold (relaxed) | 257/257 | All targets complete including former OOMs |
| AlphaFold (unrelaxed) | 257/257 | All targets complete |
| Boltz-1 | 122/257 correct | 135 re-running with deduplicated FASTAs (Job 9304637) |
| Rosetta relax | Running | Job 9195328 + 9292713 (needs restart after Boltz completes) |
| MolProbity | Pending | After Rosetta completion |
| PoseBusters | Pending | After Rosetta completion |

### Chain Deduplication Fix (2026-03-03)

**Root cause:** The original Boltz FASTA generation script expanded RCSB headers like "Chains A, D" into separate entries for each chain letter, creating duplicate sequences for homo-multimeric targets.

**Impact:** 135 of 257 targets had Boltz predicting larger complexes than intended (e.g., 1AHW: 6 chains instead of 3). Additionally, 7 targets classified as "OOM" were only OOM due to these inflated chain counts.

**Fix:** Regenerated all `boltz_input.fasta` from `sequence.fasta` (unique chains only). Boltz re-predictions submitted as Job 9304637.

## References

1. Jumper, J. et al. Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583-589 (2021).
2. Wohlwend, J. et al. Boltz-1: Democratizing Biomolecular Interaction Modeling. *bioRxiv* (2024).
3. Vreven, T. et al. Updates to the Integrated Protein-Protein Interaction Benchmarks. *J. Mol. Biol.* 427, 3031-3041 (2015).

## License

MIT License

---
*Last updated: 2026-03-03*
