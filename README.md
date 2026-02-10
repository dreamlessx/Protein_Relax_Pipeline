# Protein Structure Prediction & Relaxation Pipeline

A comprehensive benchmarking framework for evaluating AI-based protein structure prediction methods against the Protein-Protein Docking Benchmark 5.5 (BM5.5) dataset.

## Overview

This repository contains structure predictions and relaxation protocols for protein-protein complexes from the BM5.5 benchmark. We systematically compare AlphaFold 2.3.2 and Boltz-1 predictions against experimental crystal structures, with subsequent relaxation across multiple scoring functions.

## Dataset

**Source:** [Protein-Protein Docking Benchmark 5.5](https://zlab.wenglab.org/benchmark/)

| Category | Count |
|----------|-------|
| Total complexes | 257 |
| Rigid-body | 162 |
| Medium difficulty | 60 |
| Difficult | 35 |

### Non-Standard BM5.5 Entries

BM5.5 includes 4 non-standard entries representing alternate chain combinations from the same PDB structures:

| ID | Parent PDB | Description |
|----|------------|-------------|
| BAAD | 3AAD | Alternate chains A:B (Double bromodomain + ASF1) |
| BOYV | 1OYV | Chains B:I (Subtilisin + Tomato inhibitor) |
| BP57 | 3P57 | Chains AB:P (MEF2A dimer + p300 TAZ2) |
| CP57 | 3P57 | Chains CD:P (MEF2A dimer + p300 TAZ2) |

These are counted separately from their parent structures (3AAD, 1OYV, 3P57) as distinct docking cases.

## Structure Prediction Methods

### AlphaFold 2.3.2
- Database preset: `reduced_dbs` (UniRef30 via MMseqs2)
- Model preset: Auto-detect (monomer for 1 chain, multimer for 2+ chains)
- 5 ranked models per target
- **Both unrelaxed and AMBER-relaxed versions saved** (`--models_to_relax=all`)
- Template date: unrestricted
- Memory: 60 GB RAM, NVIDIA RTX A6000 GPU

### Boltz-1 v0.4.1
- MSA server for sequence alignments
- 10 recycling steps
- 200 sampling steps
- 5 diffusion samples per target
- Outputs are unrelaxed (native Boltz predictions)

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
│   ├── boltz_input.fasta
│   └── sequence.fasta
└── ...
```

## Relaxation Protocols

### AMBER Relaxation (AlphaFold Native)

AlphaFold's built-in AMBER relaxation using OpenMM:
- Force field: AMBER ff14SB
- Energy tolerance: 2.39 kcal/mol
- Position restraint stiffness: 10.0 kcal/mol/A^2
- GPU-accelerated (`--use_gpu_relax`)

### Rosetta Relaxation Protocols

Six Rosetta relaxation protocols with 5 replicates each:

| Protocol | Description |
|----------|-------------|
| `cartesian_beta` | Cartesian space, beta_nov16 scoring |
| `cartesian_ref15` | Cartesian space, REF2015 scoring |
| `dualspace_beta` | Dual space with bond geometry optimization, beta_nov16 |
| `dualspace_ref15` | Dual space with bond geometry optimization, REF2015 |
| `normal_beta` | Torsion space, beta_nov16 scoring |
| `normal_ref15` | Torsion space, REF2015 scoring |

### Summary of Relaxation Types

| # | Type | Method | Applied To |
|---|------|--------|------------|
| 1 | AMBER (native) | AlphaFold OpenMM | AlphaFold predictions |
| 2 | cartesian_beta | Rosetta | All predictions |
| 3 | cartesian_ref15 | Rosetta | All predictions |
| 4 | dualspace_beta | Rosetta | All predictions |
| 5 | dualspace_ref15 | Rosetta | All predictions |
| 6 | normal_beta | Rosetta | All predictions |
| 7 | normal_ref15 | Rosetta | All predictions |

## Computational Resources

All predictions generated on the ACCRE high-performance computing cluster at Vanderbilt University.

- **Partition:** batch_gpu (csb_gpu_acc)
- **GPU:** NVIDIA RTX A6000 / L40S
- **Memory:** 40-64 GB per job

## References

1. Jumper, J. et al. Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583–589 (2021).
2. Wohlwend, J. et al. Boltz-1: Democratizing Biomolecular Interaction Modeling. *bioRxiv* (2024).
3. Vreven, T. et al. Updates to the Integrated Protein–Protein Interaction Benchmarks. *J. Mol. Biol.* 427, 3031–3041 (2015).

## License

MIT License

## Technical Notes & Troubleshooting

### HHblits/BFD Memory Issues

**Problem:** AlphaFold crashes during MSA generation on antibody/immunoglobulin sequences with error:
```
RuntimeError: HHblits failed
WARNING: maximum number of residues 32763 exceeded
```

**Cause:** Immunoglobulin domains match titin-like proteins in the BFD metagenomic database, causing HHblits to exceed memory limits.

**Solution:** Use `--db_preset=reduced_dbs` with `--small_bfd_database_path` instead of full BFD. This bypasses HHblits and uses MMseqs2 for UniRef30 searches.

**Impact:** Minimal accuracy loss (~0.5-1% TM-score) for well-characterized protein families. Antibodies have extensive homologs in UniRef90, making BFD sequences redundant.

### Boltz-1 FASTA Format

**Requirement:** Boltz-1 expects headers in format `>CHAIN|protein|description`

**Example:**
```
>A|protein|chain A
MKTAYIAKQRQISFVKSH...
>B|protein|chain B
DIVLTQSPASLAVSLGQR...
```

**Note:** Simple headers like `>A` will fail with "Invalid record id" error. The `|protein|` designation is required.

### Storage Management

Large intermediate files are cleaned after successful prediction:
- `*.sto` (MSA alignments) - 50-500 MB each
- `features.pkl` - 100-500 MB
- `result_model_*.pkl` - 200-800 MB each

Only ranked PDB outputs and FASTA files are retained.

## Run Log

| Date | Event | Details |
|------|-------|---------|
| 2026-02-07 | Initial batch | AF2 with reduced_dbs preset |
| 2026-02-08 | Cleanup | Removed 15 non-BM5.5 targets from older benchmarks |
| 2026-02-08 | Added | 4 non-standard BM5.5 entries (BAAD, BOYV, BP57, CP57) |
| 2026-02-09 | Full re-run | All 257 targets with native AMBER relaxation, saving both relaxed and unrelaxed |

### Coverage Summary

| Metric | Count |
|--------|-------|
| Standard BM5.5 entries | 253 |
| Non-standard BM5.5 entries | 4 |
| **Total** | **257** |

### Current Progress

**AlphaFold:** Re-running all 257 targets (Job 8849208)
- Saving both unrelaxed and AMBER-relaxed versions
- Native AlphaFold AMBER relaxation (OpenMM ff14SB)

**Boltz-1:** Pending (after AlphaFold completion)

---
*Last updated: 2026-02-09*
