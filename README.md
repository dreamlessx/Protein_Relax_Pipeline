# Protein Structure Prediction & Relaxation Pipeline

A comprehensive benchmarking framework for evaluating AI-based protein structure prediction methods against the Protein-Protein Docking Benchmark 5.5 (BM5.5) dataset.

## Overview

This repository contains structure predictions and relaxation protocols for protein-protein complexes from the BM5.5 benchmark. We systematically compare AlphaFold 2.3.2 and Boltz-1 predictions against experimental crystal structures, with subsequent Rosetta relaxation across multiple scoring functions.

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
- 5 AMBER-relaxed ranked models per target (`--use_gpu_relax --models_to_relax=all`)
- Template date: unrestricted
- Memory: 60 GB RAM, NVIDIA RTX A6000 GPU

### Boltz-1 v0.4.1
- MSA server for sequence alignments
- 10 recycling steps
- 200 sampling steps  
- 5 diffusion samples per target

## Directory Structure

```
data/
├── {PDB_ID}/
│   ├── af_out/
│   │   ├── ranked_0.pdb
│   │   ├── ranked_1.pdb
│   │   ├── ranked_2.pdb
│   │   ├── ranked_3.pdb
│   │   └── ranked_4.pdb
│   ├── boltz_out/
│   │   ├── boltz_input_model_0.pdb
│   │   ├── boltz_input_model_1.pdb
│   │   ├── boltz_input_model_2.pdb
│   │   ├── boltz_input_model_3.pdb
│   │   └── boltz_input_model_4.pdb
│   ├── boltz_input.fasta
│   └── sequence.fasta
└── ...

test_subset/
├── {PDB_ID}/
│   ├── AF/                    # AlphaFold predictions
│   ├── Boltz/                 # Boltz-1 predictions
│   ├── relax/                 # Rosetta relaxation results
│   │   ├── AF/
│   │   └── Boltz/
│   ├── cartesian_beta/        # Relaxed experimental structure
│   ├── cartesian_ref15/
│   ├── dualspace_beta/
│   ├── dualspace_ref15/
│   ├── normal_beta/
│   └── normal_ref15/
└── ...
```

## Relaxation Protocols

Six Rosetta relaxation protocols with 5 replicates each:

| Protocol | Description |
|----------|-------------|
| `cartesian_beta` | Cartesian space, beta_nov16 scoring |
| `cartesian_ref15` | Cartesian space, REF2015 scoring |
| `dualspace_beta` | Dual space with bond geometry optimization, beta_nov16 |
| `dualspace_ref15` | Dual space with bond geometry optimization, REF2015 |
| `normal_beta` | Torsion space, beta_nov16 scoring |
| `normal_ref15` | Torsion space, REF2015 scoring |

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

**Requirement:** Boltz-1 expects headers in format `>CHAIN|PROTEIN|description`

**Example:**
```
>A|PROTEIN|1AK4 chain A
MKTAYIAKQRQISFVKSH...
>B|PROTEIN|1AK4 chain B
DIVLTQSPASLAVSLGQR...
```

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
| 2026-02-08 | Script update | Added `--models_to_relax=all` for full AMBER relaxation |

### Coverage Summary

| Metric | Count |
|--------|-------|
| Standard BM5.5 entries | 253 |
| Non-standard BM5.5 entries | 4 |
| **Total** | **257** |

### Current Progress

**Status:** 247/257 complete (96%)

**Remaining (10 targets):**

| Target | AlphaFold | Boltz-1 |
|--------|-----------|---------|
| 1WEJ | pending | done |
| 2BTF | pending | done |
| 4CPA | pending | done |
| 5JMO | pending | done |
| 4GXU | done | pending |
| 5HYS | done | pending |
| BAAD | done | pending |
| BOYV | done | pending |
| BP57 | done | pending |
| CP57 | done | pending |

### AMBER Relaxation Note

Predictions completed before 2026-02-08 have only `ranked_0.pdb` AMBER-relaxed (default AlphaFold behavior). Predictions after this date use `--models_to_relax=all` to relax all 5 ranked models.

---
*Last updated: 2026-02-08*
