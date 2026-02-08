#!/bin/bash --norc
#
# Run AlphaFold 2.3.2 on a single target directory.
# Auto-detects FASTA (prefers sequence.fasta, falls back to boltz_input.fasta)
# and sets model_preset to monomer/multimer based on sequence count.
#
# Usage:
#   sbatch --chdir=<target_dir> alphafold_single.slurm [fasta_basename]
#
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=batch_gpu
#SBATCH --constraint=csbtmp
#SBATCH --mail-user=mudit.agarwal@vanderbilt.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --gres=gpu:nvidia_rtx_a6000:1
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --job-name=af232
#SBATCH --output=%x.out

set -euo pipefail

# ACCRE paths
AF2_MINICONDA=/sb/apps/alphafold232/miniconda3
AF2_REPO=/sb/apps/alphafold232/alphafold
AF2_DATADIR=/sb/apps/alphafold-data.230

echo "Node: $SLURM_JOB_NODELIST"
nvidia-smi || true

CALCDIR=.
cd "$CALCDIR"
echo "CWD  : $(pwd)"

# Pick the FASTA
FASTA_BASENAME="${1:-}"
if [[ -z "$FASTA_BASENAME" ]]; then
  if   [[ -s sequence.fasta ]];      then FASTA_BASENAME="sequence.fasta"
  elif [[ -s boltz_input.fasta ]];   then FASTA_BASENAME="boltz_input.fasta"
  else
    echo "No sequence.fasta or boltz_input.fasta in $(pwd); provide a basename as arg." >&2
    exit 2
  fi
else
  [[ "$FASTA_BASENAME" == */* ]] && { echo "Pass a basename, not a path: $FASTA_BASENAME" >&2; exit 2; }
  [[ -s "$FASTA_BASENAME" ]] || { echo "FASTA not found here: $FASTA_BASENAME" >&2; exit 2; }
fi
echo "FASTA: $FASTA_BASENAME"

# Decide model preset
SEQCOUNT=$(grep -c '^>' "$FASTA_BASENAME" || true)
if (( SEQCOUNT > 1 )); then
  MODEL_PRESET=multimer
else
  MODEL_PRESET=monomer
fi
echo "MODEL_PRESET: $MODEL_PRESET (seqs: $SEQCOUNT)"

OUTDIR="$CALCDIR/af_out"
mkdir -p "$OUTDIR"

# Activate env
source "$AF2_MINICONDA/bin/activate" af232
export LD_LIBRARY_PATH="$AF2_MINICONDA/envs/af232/lib:${LD_LIBRARY_PATH:-}"

# Build args
ARGS=(
  --fasta_paths="$FASTA_BASENAME"
  --output_dir="$OUTDIR"
  --data_dir="$AF2_DATADIR"
  --max_template_date=9999-12-31
  --model_preset="$MODEL_PRESET"
  --db_preset=full_dbs
  --use_gpu_relax
  --uniref90_database_path="$AF2_DATADIR/uniref90/uniref90.fasta"
  --mgnify_database_path="$AF2_DATADIR/mgnify/mgy_clusters_2022_05.fa"
  --uniref30_database_path="$AF2_DATADIR/uniref30/UniRef30_2021_03"
  --bfd_database_path="$AF2_DATADIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
  --template_mmcif_dir="$AF2_DATADIR/pdb_mmcif/mmcif_files"
  --obsolete_pdbs_path="$AF2_DATADIR/pdb_mmcif/obsolete.dat"
)

# Preset-specific DB flags
if [[ "$MODEL_PRESET" == "monomer" ]]; then
  ARGS+=( --pdb70_database_path="$AF2_DATADIR/pdb70/pdb70" )
else
  ARGS+=( --pdb_seqres_database_path="$AF2_DATADIR/pdb_seqres/pdb_seqres.txt" )
  ARGS+=( --uniprot_database_path="$AF2_DATADIR/uniprot/uniprot.fasta" )
  ARGS+=( --num_multimer_predictions_per_model=1 )
fi

set -x
python "$AF2_REPO/run_alphafold.py" "${ARGS[@]}"
rc=$?
set +x
echo "AlphaFold exit code: $rc"
exit $rc
