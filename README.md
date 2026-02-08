# Protein Relaxation Benchmarking Pipeline

A systematic pipeline for benchmarking AI-predicted protein-protein complex structures refined with Rosetta relaxation protocols. Compares AlphaFold 2.3.2 and Boltz-1 predictions across six relaxation protocols, validated with MolProbity.

---

## Table of Contents

1. [Overview](#overview)
2. [Dataset](#dataset)
3. [Prerequisites](#prerequisites)
4. [Pipeline Stages](#pipeline-stages)
   - [Stage 1: Data Preparation](#stage-1-data-preparation)
   - [Stage 2: Structure Prediction](#stage-2-structure-prediction)
   - [Stage 3: Rosetta Relaxation](#stage-3-rosetta-relaxation)
   - [Stage 4: Structural Validation](#stage-4-structural-validation)
   - [Stage 5: Metric Collection](#stage-5-metric-collection)
5. [Directory Structure](#directory-structure)
6. [Reproducing](#reproducing)
7. [Acknowledgments](#acknowledgments)

---

## Overview

This pipeline predicts protein-protein complex structures using two AI methods and refines them with six Rosetta relaxation protocols to determine which combinations yield the most physically realistic models.

- **2 prediction methods:** AlphaFold 2.3.2 (multimer) and Boltz-1 v0.4.1
- **6 relaxation protocols:** 3 modes (Cartesian, Dualspace, Normal) x 2 scoring functions (beta_nov16, ref2015)
- **5 replicates** per protocol per structure
- **MolProbity validation** for structural quality assessment

## Dataset

**234 protein-protein complexes** from the [Protein-Protein Docking Benchmark 5.5](https://zlab.umassmed.edu/benchmark/) (Weng Lab, UMass):

| Category | Count | Description |
|----------|-------|-------------|
| AA | 48 | Antibody-Antigen |
| AS | 9 | Antigen-Single Domain Ab |
| EI | 9 | Enzyme-Inhibitor |
| ER | 12 | Enzyme-Regulatory |
| ES | 6 | Enzyme-Substrate |
| OG | 1 | Others (G-Protein) |
| OR | 6 | Others (Receptor) |
| OX | 18 | Others (Miscellaneous) |

The BM5.5 provides structures in a split format per complex:

| File Suffix | Meaning |
|-------------|---------|
| `_r_b.pdb` | Receptor, bound conformation (from complex crystal) |
| `_r_u.pdb` | Receptor, unbound conformation (separate crystal) |
| `_l_b.pdb` | Ligand, bound conformation |
| `_l_u.pdb` | Ligand, unbound conformation |

For this pipeline, bound conformations are merged into complete complexes.

---

## Prerequisites

### Software

| Software | Version | Purpose | Installation |
|----------|---------|---------|-------------|
| Python | 3.10+ | Data processing | System or conda |
| AlphaFold | 2.3.2 | Structure prediction | [GitHub](https://github.com/deepmind/alphafold) |
| Boltz-1 | 0.4.1 | Structure prediction | [GitHub](https://github.com/jwohlwend/boltz) |
| Rosetta | 3.14 | Relaxation protocols | [RosettaCommons](https://www.rosettacommons.org/) |
| Phenix | Latest | MolProbity validation | [Phenix](https://www.phenix-online.org/) |
| Reduce | Latest | Hydrogen addition | Included with Phenix or [SBGrid](https://sbgrid.org/) |
| PyMOL | 2.x+ | RMSD/energy analysis | [PyMOL](https://pymol.org/) |
| requests | Latest | FASTA downloads | `pip install requests` |

### ACCRE-Specific Paths (Vanderbilt)

```bash
# AlphaFold 2.3.2
AF2_MINICONDA=/sb/apps/alphafold232/miniconda3
AF2_REPO=/sb/apps/alphafold232/alphafold
AF2_DATADIR=/sb/apps/alphafold-data.230

# Boltz-1 v0.4.1
BOLTZ_MINICONDA=/sb/apps/boltz1-v0.4.1/miniconda3

# Rosetta 3.14
RELAX=/dors/meilerlab/apps/rosetta/rosetta-3.14/main/source/bin/relax.linuxgccrelease
ROSETTA_DB=/dors/meilerlab/apps/rosetta/rosetta-3.14/main/database
CLEAN_SCRIPT=/dors/meilerlab/apps/rosetta/rosetta-3.14/main/tools/protein_tools/scripts/clean_pdb.py

# MolProbity (SBGrid)
REDUCE=/programs/x86_64-linux/system/sbgrid_bin/reduce
PHENIX=/programs/x86_64-linux/system/sbgrid_bin/phenix.molprobity
```

### Compute Resources

| Job Type | GPU | Memory | Walltime |
|----------|-----|--------|----------|
| AlphaFold (single) | NVIDIA RTX A6000 | 100 GB | 72 h |
| AlphaFold (array) | NVIDIA RTX A6000 | 6 GB | 48 h |
| Boltz-1 | NVIDIA L40S | 256 GB | 2 days |
| Rosetta Relax | CPU only | 4 GB | 48 h |

---

## Pipeline Stages

### Stage 1: Data Preparation

#### 1a. Download and Merge BM5.5 Structures

Download from [https://zlab.umassmed.edu/benchmark/](https://zlab.umassmed.edu/benchmark/) and merge bound partner structures into complete complexes:

```bash
# Merge bound receptor + ligand into one PDB per complex
for pdb_id in $(ls structures/*_r_b.pdb | sed 's/.*\///' | sed 's/_r_b.pdb//'); do
    cat "structures/${pdb_id}_r_b.pdb" "structures/${pdb_id}_l_b.pdb" | \
        grep -E '^(ATOM|HETATM|TER)' > "merged/${pdb_id}.pdb"
    echo "END" >> "merged/${pdb_id}.pdb"
done
```

Or download full complex PDBs directly from RCSB:
```bash
for pdb_id in $(cat pdb_list.txt); do
    wget -q "https://files.rcsb.org/download/${pdb_id}.pdb" -O "merged/${pdb_id}.pdb"
    sleep 0.2
done
```

**Expected output:** ~230-266 PDB files in `merged/`

#### 1b. Clean PDB Structures

Standardizes atom naming, removes HETATM records, and merges cleaned chains into a single file per complex.

```bash
bash scripts/data_preparation/clean_pdbs.sh <merged_dir> <cleaned_dir> <rosetta_clean_script>
```

**What the script does:**
1. Extracts unique chain IDs from ATOM records in each PDB
2. Runs Rosetta's `clean_pdb.py` on each chain individually
3. Merges cleaned chains into one file with TER records between chains
4. Renumbers atom serials sequentially from 1

**Expected output:** One `{PDBID}.pdb` per complex in `cleaned/`, same count as `merged/`

**Verify:**
```bash
ls cleaned/*.pdb | wc -l          # Same count as merged/
tail -3 cleaned/1AK4.pdb          # Should end with TER, then END
head -1 cleaned/1AK4.pdb          # Atom serial should start at 1
```

#### 1c. Download FASTA Sequences

Downloads sequences from RCSB with fallback to PDBe. Handles obsolete entries by looking up replacements.

```bash
python3 scripts/data_preparation/download_fastas.py merged/ fasta/
```

**What the script does:**
1. Extracts 4-character PDB IDs from filenames
2. Tries three download endpoints: RCSB primary, RCSB fallback, PDBe
3. For obsolete entries, queries RCSB API for replacement PDB IDs
4. Retries with exponential backoff; 0.15s delay between requests
5. Logs all results to `fasta_download_log.csv`

**Expected output:** One `.fasta` per PDB in `fasta/`, plus `fasta_download_log.csv`

**Verify:**
```bash
grep -c 'ok' fasta/fasta_download_log.csv      # ~260-266 successes
grep 'fail' fasta/fasta_download_log.csv        # 0-3 failures
head -2 fasta/1AK4.fasta                        # Should start with >
```

#### 1d. Organize FASTAs into Per-PDB Directories

```bash
python3 scripts/data_preparation/organize_fastas.py fasta/ data/ --rename sequence.fasta
```

Creates `data/{PDBID}/sequence.fasta` for each downloaded FASTA.

**Verify:**
```bash
ls data/ | wc -l                                # Same count as FASTAs
ls data/1AK4/sequence.fasta                     # Should exist
```

#### 1e. Convert to Boltz-1 Input Format

Boltz-1 requires headers in the format `>{CHAIN}|PROTEIN|`. This script parses chain IDs from RCSB headers and rewrites them.

```bash
python3 scripts/data_preparation/prepare_boltz_fastas.py data/
```

**What the script does:**
1. Reads each `data/*/sequence.fasta`
2. Parses chain IDs from headers (looks for `Chains A, B` pattern)
3. Rewrites as `>A|PROTEIN|`, `>B|PROTEIN|`, etc.
4. Falls back to chain `A` if header format is unrecognized

**Expected output:** `boltz_input.fasta` in each `data/{PDBID}/` directory

**Verify:**
```bash
head -4 data/1AK4/boltz_input.fasta             # >A|PROTEIN| format
# Validate all Boltz FASTAs:
for f in data/*/boltz_input.fasta; do
    if ! grep -q '|PROTEIN|' "$f"; then echo "BAD: $f"; fi
done
```

---

### Stage 2: Structure Prediction

#### AlphaFold 2.3.2

**Array job** (all targets):
```bash
N=$(ls -d data/*/ | wc -l)
sbatch --array=1-${N}%5 scripts/prediction/alphafold_array.slurm
```

**Single target:**
```bash
sbatch --chdir=data/1AK4 scripts/prediction/alphafold_single.slurm
```

**What the scripts do:**
- Auto-detect monomer vs multimer from sequence count (`>` lines)
- Run AlphaFold with `full_dbs` preset and GPU relaxation
- Monomer: uses `pdb70` database
- Multimer: uses `pdb_seqres` + `uniprot`, 1 prediction per model (5 total)
- Completion guard: skips targets with >= 5 ranked models already

**Key AlphaFold flags:**
```
--db_preset=full_dbs
--use_gpu_relax
--max_template_date=9999-12-31
--num_multimer_predictions_per_model=1  (multimer only)
```

**AlphaFold databases used:**
- `uniref90/uniref90.fasta`
- `mgnify/mgy_clusters_2022_05.fa`
- `uniref30/UniRef30_2021_03`
- `bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt`
- `pdb70/pdb70` (monomer) or `pdb_seqres/pdb_seqres.txt` + `uniprot/uniprot.fasta` (multimer)
- `pdb_mmcif/mmcif_files/` + `pdb_mmcif/obsolete.dat`

**Expected output per target:**
```
data/{PDBID}/af_out/sequence/
├── ranked_0.pdb through ranked_4.pdb   # 5 ranked models
├── ranking_debug.json                   # pLDDT/pTM scores
├── timings.json                         # Runtime
├── features.pkl                         # Input features (~large)
├── msas/                                # MSA alignments (~140 MB)
├── result_model_*.pkl                   # Model outputs (~50 MB each)
├── unrelaxed_model_*.pdb                # Pre-relaxation
└── relaxed_model_*.pdb                  # AMBER-relaxed
```

**Verify:**
```bash
# Count completed targets
for d in data/*/af_out/sequence/; do
    [ $(ls "$d"/ranked_*.pdb 2>/dev/null | wc -l) -ge 5 ] && echo "OK: $(basename $(dirname $(dirname $d)))"
done | wc -l
```

#### Boltz-1 v0.4.1

**Array job** (all targets):
```bash
N=$(ls -d data/*/ | wc -l)
sbatch --array=1-${N} scripts/prediction/boltz_array.slurm
```

**Single target:**
```bash
sbatch --chdir=data/1AK4 scripts/prediction/boltz_single.slurm boltz_input.fasta
```

**What the scripts do:**
- Use MSA server for fast multiple sequence alignment (`--use_msa_server`)
- Generate 5 diffusion samples per target
- Validate FASTA headers match `>X|PROTEIN|` format before running
- Skip guard checks for existing models (override with `FORCE=1`)

**Key Boltz parameters:**
```
--diffusion_samples 5
--recycling_steps 10
--sampling_steps 200
--output_format pdb
--use_msa_server
```

**Expected output per target:**
```
data/{PDBID}/boltz_out_dir/boltz_results_boltz_input/
├── predictions/boltz_input/
│   ├── boltz_input_model_0.pdb through model_4.pdb
│   ├── confidence_*.json
│   └── plddt_*.npz
├── msa/
└── processed/
```

**Verify:**
```bash
for d in data/*/boltz_out_dir/; do
    count=$(find "$d" -name 'boltz_input_model_*.pdb' 2>/dev/null | wc -l)
    [ "$count" -ge 5 ] && echo "OK"
done | wc -l
```

---

### Stage 3: Rosetta Relaxation

#### Organize Predictions

Before relaxation, copy predictions into the `test/` directory structure:

```bash
for pdb_id in $(ls data/); do
    mkdir -p test/$pdb_id/AF test/$pdb_id/Boltz

    # Crystal structure
    [ -f cleaned/$pdb_id.pdb ] && cp cleaned/$pdb_id.pdb test/$pdb_id/$pdb_id.pdb

    # AF predictions
    cp data/$pdb_id/af_out/sequence/ranked_*.pdb test/$pdb_id/AF/ 2>/dev/null

    # Boltz predictions
    find data/$pdb_id/boltz_out_dir -name 'boltz_input_model_*.pdb' \
        -exec cp {} test/$pdb_id/Boltz/ \; 2>/dev/null
done
```

#### Run Relaxation

```bash
N=$(ls -d test/*/ | wc -l)
sbatch --array=1-${N} scripts/relaxation/relax_predictions.slurm
```

**What the script does:**
1. For each PDB directory, finds the crystal structure
2. Runs 6 relaxation protocols x 5 replicates = 30 relaxed structures
3. Each replicate uses `-nstruct 1` with a different suffix (`_r1` through `_r5`)

#### The 6 Protocols

| Protocol | Mode | Scoring Function | Rosetta Flags |
|----------|------|-----------------|---------------|
| `cartesian_beta` | Cartesian | beta_nov16 | `-relax:cartesian -beta_nov16 -score:weights beta_nov16_cart` |
| `cartesian_ref15` | Cartesian | ref2015 | `-relax:cartesian -score:weights ref2015_cart` |
| `dualspace_beta` | Dualspace | beta_nov16 | `-relax:dualspace -beta_nov16 -score:weights beta_nov16_cart -nonideal -relax:minimize_bond_angles -relax:minimize_bond_lengths` |
| `dualspace_ref15` | Dualspace | ref2015 | `-relax:dualspace -score:weights ref2015_cart -nonideal -relax:minimize_bond_angles -relax:minimize_bond_lengths` |
| `normal_beta` | Normal FastRelax | beta_nov16 | `-beta_nov16 -score:weights beta_nov16` |
| `normal_ref15` | Normal FastRelax | ref2015 | `-score:weights ref2015` |

**Common Rosetta flags (all protocols):**
```
-ignore_zero_occupancy false    # Process all atoms
-nstruct 1                      # One structure per run
-no_nstruct_label               # Clean output naming
-out:pdb_gz                     # Gzip output (saves disk)
-flip_HNQ                       # Optimize His/Asn/Gln hydrogens
-fa_max_dis 9.0                 # Interaction distance cutoff
-optimization:default_max_cycles 200   # Convergence iterations
```

**Expected output per PDB (crystal relaxation):**
```
test/{PDBID}/
├── cartesian_beta/
│   ├── {PDBID}_r1.pdb.gz ... _r5.pdb.gz    # 5 replicates
│   ├── log/{PDBID}_cart_beta_r1.log ... _r5.log
│   └── relax.fasc                            # Rosetta scores
├── cartesian_ref15/
├── dualspace_beta/
├── dualspace_ref15/
├── normal_beta/
└── normal_ref15/
```

**For AI prediction relaxation (extended):**
```
test/{PDBID}/relax/
├── AF/ranked_0/{protocol}/ranked_0_r{1-5}.pdb.gz
├── AF/ranked_1/...
├── ...
├── Boltz/boltz_input_model_0/{protocol}/...
└── ...
```

**Verify:**
```bash
# Crystal relaxation: 30 files per PDB (6 protocols x 5 replicates)
find test/1AK4/ -maxdepth 2 -name '*.pdb.gz' | wc -l

# AI prediction relaxation: 300 files per PDB (10 models x 6 x 5)
find test/1AK4/relax/ -name '*.pdb.gz' | wc -l

# Check a log for errors
tail -5 test/1AK4/cartesian_beta/log/1AK4_cart_beta_r1.log
```

---

### Stage 4: Structural Validation

```bash
bash scripts/validation/run_molprobity.sh [base_dir] [output_dir]
```

**What the script does:**
1. For each PDB, validates all structures: crystal, raw predictions, and relaxed structures
2. Decompresses `.pdb.gz` files as needed
3. Adds hydrogens with `reduce -FLIP -Quiet` (falls back to `-Quiet` only, then raw input)
4. Runs `phenix.molprobity` with 600-second timeout per structure
5. Parses output for 16 metrics
6. Writes per-protein CSV files (supports resume — skips already-validated structures)

**Structures validated per PDB:**
- 1 crystal structure
- 10 raw predictions (5 AF + 5 Boltz)
- 30 crystal relaxations (6 protocols x 5 replicates)
- 300 prediction relaxations (10 models x 6 protocols x 5 replicates)
- **341 structures per PDB**

**Metrics collected:**

| Metric | Description | Good Values |
|--------|-------------|-------------|
| Clashscore | Steric clashes per 1000 atoms | < 10 |
| Ramachandran Favored | % in favored backbone regions | > 95% |
| Ramachandran Outliers | % in outlier backbone regions | < 0.5% |
| Rotamer Outliers | % side-chain rotamer outliers | < 2% |
| C-beta Deviations | Backbone geometry deviations | < 5 |
| Bond RMSZ | Bond length deviation from ideal | < 1.5 |
| Angle RMSZ | Bond angle deviation from ideal | < 1.5 |
| MolProbity Score | Combined quality metric | < 2.0 |

**Verify:**
```bash
wc -l molprobity_results/csvs/1AK4_validation.csv     # 342 (1 header + 341)
grep -c 'FAILED' molprobity_results/csvs/*.csv          # Should be 0

# Merge all CSVs:
head -1 molprobity_results/csvs/$(ls molprobity_results/csvs/ | head -1) > all_validation.csv
tail -n +2 -q molprobity_results/csvs/*_validation.csv >> all_validation.csv
```

---

### Stage 5: Metric Collection

```python
# Inside PyMOL:
run scripts/analysis/collect_metrics.py
collect_metrics ref=crystal, sel=polymer and name CA, do_fit=1, out=metrics.tsv
```

**What the script does:**
1. Aligns each loaded object to the reference structure using C-alpha atoms
2. Calculates RMSD with `cmd.align()` (no outlier rejection, `cycles=0`)
3. Extracts Rosetta `total_score` from PDB REMARK lines or external scorefile
4. Outputs TSV with: `object`, `rmsd_to_{ref}`, `pairs`, `energy_total_score`

**With Rosetta scorefile:**
```python
collect_metrics ref=crystal, scorefile=test/1AK4/cartesian_beta/relax.fasc, out=metrics.tsv
```

---

## Directory Structure

```
Protein_Relax_Pipeline/
├── README.md
├── LICENSE
├── scripts/
│   ├── data_preparation/
│   │   ├── download_fastas.py         # Download FASTA from RCSB/PDBe
│   │   ├── organize_fastas.py         # Organize into per-PDB directories
│   │   ├── prepare_boltz_fastas.py    # Convert to Boltz-1 input format
│   │   └── clean_pdbs.sh              # Clean PDBs with Rosetta
│   ├── prediction/
│   │   ├── alphafold_array.slurm      # Array job: AlphaFold 2.3.2
│   │   ├── alphafold_single.slurm     # Single target: AlphaFold 2.3.2
│   │   ├── boltz_array.slurm          # Array job: Boltz-1
│   │   └── boltz_single.slurm         # Single target: Boltz-1
│   ├── relaxation/
│   │   └── relax_predictions.slurm    # 6 protocols x 5 replicates
│   ├── validation/
│   │   └── run_molprobity.sh          # MolProbity validation
│   └── analysis/
│       └── collect_metrics.py         # PyMOL RMSD/energy collection
├── data/                              # Per-PDB prediction data (147 targets)
│   └── {PDB_ID}/
│       ├── sequence.fasta             # RCSB format
│       ├── boltz_input.fasta          # Boltz-1 format
│       ├── af_out/ranked_{0-4}.pdb    # AlphaFold predictions
│       └── boltz_out/boltz_input_model_{0-4}.pdb  # Boltz-1 predictions
└── test_subset/                       # Complete 20-PDB benchmark set
    └── {PDB_ID}/
        ├── {PDB_ID}.pdb               # Experimental crystal structure
        ├── AF/ranked_{0-4}.pdb        # AlphaFold predictions
        ├── Boltz/boltz_input_model_{0-4}.pdb  # Boltz-1 predictions
        ├── {protocol}/                # Crystal relaxations (5 replicates each)
        │   ├── {PDB}_r{1-5}.pdb.gz
        │   └── log/
        └── relax/
            ├── AF/{model}/{protocol}/ # Relaxed AF predictions
            └── Boltz/{model}/{protocol}/  # Relaxed Boltz predictions
```

---

## Reproducing

1. Clone this repository
2. Install required software (see [Prerequisites](#prerequisites))
3. Update SLURM account, partition, and path settings in scripts for your cluster
4. Follow pipeline stages in order
5. Use verification checkpoints after each step to confirm outputs

For a complete step-by-step reproducibility walkthrough with script audits, see the companion repository: [Protein_Ideal](https://github.com/dreamlessx/Protein_Ideal)

## Acknowledgments

This work is conducted in the [Meiler Lab](https://www.meilerlab.org/) at Vanderbilt University, supported by ARPA-H.

## Citation

*Manuscript in preparation.*

## License

See [LICENSE](LICENSE) for details.
