# Protein Structure Prediction & Relaxation Pipeline

A comprehensive benchmarking framework for evaluating AI-based protein structure prediction methods against the Protein-Protein Docking Benchmark 5.5 (BM5.5) dataset.

## Overview

This repository contains structure predictions and relaxation protocols for protein-protein complexes from the BM5.5 benchmark. We systematically compare AlphaFold 2.3.2 and Boltz-1 predictions against experimental crystal structures, with subsequent Rosetta relaxation across multiple scoring functions.

## Dataset

**Source:** [Protein-Protein Docking Benchmark 5.5](https://zlab.umassmed.edu/benchmark/)

| Category | Count |
|----------|-------|
| Total complexes | 225 |
| Rigid-body | 162 |
| Medium difficulty | 60 |
| Difficult | 35 |

## Structure Prediction Methods

### AlphaFold 2.3.2
- Database preset: `reduced_dbs` (UniRef30 via MMseqs2)
- 5 ranked models per target
- GPU-accelerated relaxation enabled
- Template date: unrestricted

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
