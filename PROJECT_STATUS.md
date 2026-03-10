# Project Status

**Last updated:** 2026-03-09

## Dataset Summary

| Metric | Value |
|--------|-------|
| BM5.5 targets | 257 |
| Total chains | 605 |
| Total residues | 122,966 |
| Rigid-body | 162 |
| Medium difficulty | 60 |
| Difficult | 35 |

## Input Consistency

All 257 FASTAs are crystal-derived and verified identical across all prediction methods.

| Check | Status |
|-------|--------|
| Crystal == AF FASTA | 257/257 |
| AF FASTA == Boltz FASTA | 257/257 |
| No duplicate chains | 257/257 |
| No His-tags | 257/257 |

## Prediction Status

| Method | Progress | Status |
|--------|----------|--------|
| AlphaFold 2.3.2 | 257/257 | COMPLETE |
| AF built-in AMBER | 257/257 | COMPLETE (all have ranked_0..4.pdb) |
| AF unrelaxed | 257/257 | COMPLETE |
| Boltz-1 v0.4.1 | 257/257 | COMPLETE |
| Standalone AMBER (AF) | 256/257 | 1KTZ finishing (Job 9399617) |
| Standalone AMBER (Boltz) | 256/257 | 1KTZ finishing (Job 9399617) |

## Rosetta Relaxation Status

**Status:** IN PROGRESS

| Detail | Value |
|--------|-------|
| SLURM Jobs | 9373165, 9371774 |
| Targets started | ~50/257 |
| Runs complete | ~2,512 / ~200,000 |
| Input types | 6 (af_relaxed, af_unrelaxed, boltz, amber_af, amber_boltz, crystal) |
| Protocols | 6 (cartesian_beta, cartesian_ref15, dualspace_beta, dualspace_ref15, normal_beta, normal_ref15) |
| Replicates | 5 per protocol |
| Runs per target | 780 (max: 6 inputs x 6 protocols x 5 reps) |

## Pipeline Design

- **Green pipeline** = our independent verification of Blue's protocol
- Standalone AMBER = test condition; AF built-in AMBER = control condition
- AF built-in AMBER + AF unrelaxed come from the SAME AlphaFold run
- Rosetta output: .pdb.gz compressed, matching Blue's flags exactly

## Key Fixes Applied

| # | Fix | Impact |
|---|-----|--------|
| 1 | Crystal-derived FASTAs | 236 targets changed from UniProt full-length; net -3,956 residues |
| 2 | Crystal chain deduplication | 36 PDBs stripped of homo-multimer duplicate chains |
| 3 | FASTA deduplication | 135 Boltz FASTAs had duplicated homo-multimer chains |
| 4 | His-tag removal | 41 targets had expression artifacts removed |
| 5 | 1KTZ template bug | AF workaround using max_template_date=1900 |
| 6 | Rosetta AMBER naming collision | MODEL_LABELS array fix to avoid overwriting |

## Remaining Work

- [ ] 1KTZ standalone AMBER finishing (Job 9399617)
- [ ] Rosetta relaxation: ~197k runs remaining
- [ ] MolProbity validation (after Rosetta)
- [ ] PoseBusters validation (after Rosetta)
- [ ] Final analysis and comparison

## Timeline

| Date | Event |
|------|-------|
| 2026-02-07 | Initial AF batch (reduced_dbs) |
| 2026-02-09 | Full AF re-run (full_dbs + fallback) |
| 2026-02-11 | AF complete, 7 AMBER failures identified |
| 2026-02-15 | Boltz batch 1 (233/257) |
| 2026-02-20 | AMBER root cause: non-standard residues |
| 2026-02-22 | Rosetta started (Job 9195328) |
| 2026-03-02 | Boltz complete (257/257) |
| 2026-03-03 | FASTA dedup (135 targets) + consistency fix |
| 2026-03-04 | Crystal-derived FASTAs (all 257) |
| 2026-03-04 | Crystal stripping (36 PDBs) |
| 2026-03-05 | AF + Boltz re-predictions with crystal FASTAs |
| 2026-03-08 | All predictions complete (AF, Boltz, built-in AMBER) |
| 2026-03-09 | Rosetta relaxation in progress (~50/257 targets) |
