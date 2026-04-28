# Blue / Green Shared Context

**Last updated 2026-04-27, pipeline locked at 100.000% (snapshot 2026-04-27a)**

## Identities

- **Blue** = primary pipeline (Claude). Job prefix: `blue_`
- **Green** = independent verification of Blue's protocol. Job prefix: `green_`
- **Red** = analysis and metrics. Job prefix: `red_`
- Blue and Green follow identical Rosetta flags, identical inputs, identical output structure.

## Attribution

- ALL GitHub commits: `dreamlessx <259208371+dreamlessx@users.noreply.github.com>`
- NEVER add Co-Authored-By to any commit.
- GitHub repo: `dreamlessx/Protein_Relax_Pipeline`
- This applies to Blue, Green, AND Red, no exceptions.

## Locked-State Counts (snapshot 2026-04-27a)

| Pipeline | Rosetta MolProbity rows |
|----------|-------------------------|
| Blue | 208,170 / 208,170 |
| Green | 208,170 / 208,170 |
| Combined | 416,340 / 416,340 |

Per target per pipeline: 27 input structures x 6 protocols x 5 replicates = 810 runs.

7 source buckets: af_relaxed, af_unrelaxed, amber_af, amber_boltz, amber_crystal, boltz, crystal.

6 protocols: cartesian_beta, cartesian_ref15, dualspace_beta, dualspace_ref15, normal_beta, normal_ref15.

## Output Directory Structure (Blue)

```
af_work/{TARGET}/rosetta_out/
  af_relaxed_ranked_{0..4}/         # AF built-in AMBER
  af_unrelaxed_ranked_{0..4}/       # AF unrelaxed
  boltz_boltz_input_model_{0..4}/   # Boltz
  amber_af_ranked_{0..4}/           # Standalone AMBER(AF), MODEL_LABELS
  amber_boltz_model_{0..4}/         # Standalone AMBER(Boltz), MODEL_LABELS
  crystal_{TARGET}/                 # Crystal structure
  amber_crystal_{TARGET}/           # Standalone AMBER on crystal (v5 chain-split)
  each has: {protocol}/*_r{1..5}.pdb.gz
```

## Output Directory Structure (Green)

```
benchmarking/data/{TARGET}/rosetta_out/
  af_relaxed_ranked_{0..4}/
  af_unrelaxed_ranked_{0..4}/
  boltz_boltz_input_model_{0..4}/
  amber_af_ranked_{0..4}/
  amber_boltz_model_{0..4}/
  crystal_{TARGET}/
  amber_crystal_{TARGET}/
  each has: {protocol}/*_r{1..5}.pdb.gz
```

## Input Paths (Blue)

- AF relaxed: `af_work/{TARGET}/af_out/ranked_*.pdb`
- AF unrelaxed: `af_work/{TARGET}/af_out_unrelaxed/ranked_*.pdb`
- Boltz: `af_work/{TARGET}/boltz_out/boltz_results_boltz_input/predictions/boltz_input/boltz_input_model_*.pdb`
- AMBER AF: `af_work/{TARGET}/amber_out/af_unrelaxed_ranked_*/relaxed.pdb`
- AMBER Boltz: `af_work/{TARGET}/amber_out/boltz_model_*/relaxed.pdb`
- Crystal: `protein_ideal_test/benchmarking/cleaned/{TARGET}.pdb`

## Input Paths (Green)

- AF relaxed: `benchmarking/data/{TARGET}/AF/ranked_*.pdb`
- AF unrelaxed: `benchmarking/data/{TARGET}/af_out_unrelaxed/ranked_*.pdb`
- Boltz: `benchmarking/data/{TARGET}/Boltz/boltz_input_model_*.pdb`
- AMBER AF: `benchmarking/data/{TARGET}/amber_out/af_unrelaxed_${i}/relaxed.pdb`
- AMBER Boltz: `benchmarking/data/{TARGET}/amber_out/boltz_model_${i}/relaxed.pdb`
- Crystal: `benchmarking/cleaned/{TARGET}.pdb`

## Rosetta Flags (Must Match Exactly)

```
Common: -ignore_zero_occupancy false -nstruct 1 -no_nstruct_label -out:pdb_gz
        -flip_HNQ -fa_max_dis 9.0 -optimization:default_max_cycles 200
        -out:levels all:warning

cartesian_beta:   -relax:cartesian -beta_nov16 -score:weights beta_nov16_cart
cartesian_ref15:  -relax:cartesian -score:weights ref2015_cart
dualspace_beta:   -relax:dualspace -beta_nov16 -score:weights beta_nov16_cart -nonideal -relax:minimize_bond_angles -relax:minimize_bond_lengths
dualspace_ref15:  -relax:dualspace -score:weights ref2015_cart -nonideal -relax:minimize_bond_angles -relax:minimize_bond_lengths
normal_beta:      -beta_nov16 -score:weights beta_nov16
normal_ref15:     -score:weights ref2015
```

All SLURM scripts include `#SBATCH --exclude=cn1340` (1,614 instant-failure jobs were traced to that node).

## Key Script Registry

- Blue Rosetta fills: `scripts/blue_fill_all.slurm`
- Green Rosetta fills: `scripts/green_fill_all.slurm`
- Blue mopups: `scripts/blue_mopup.slurm`
- Green mopups: `scripts/green_mopup.slurm`, `scripts/green_mopup_v2.slurm`
- Blue AMBER-crystal: `scripts/blue_amber_crystal_fill.slurm`, `scripts/blue_rosetta_amber_crystal_v4.slurm`
- Green AMBER-crystal: `scripts/green_amber_crystal_fill.slurm`, `scripts/green_rosetta_amber_crystal_v4.slurm`
- Gap scanners: `scripts/blue_scan_gaps.sh`, `scripts/green_scan_gaps.sh`
- AlphaFold: `scripts/af_consistency_rerun.slurm`
- Boltz: `scripts/boltz_consistency_rerun.slurm`
- Standalone AMBER: `scripts/amber_relax.py`, `scripts/amber_relax.slurm`

## Cluster Notes

- SLURM account: `p_csb_meiler`, partition: `batch`
- Rosetta binary: `/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/source/bin/relax.linuxgccrelease`
- Bad node: cn1340 (excluded in all SLURM scripts)
