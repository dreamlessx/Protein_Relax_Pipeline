#!/bin/bash --norc
#
# Array job for AlphaFold 2.3.2 structure prediction.
# Runs both monomer and multimer presets for each target.
# Skips targets that already have >= 5 ranked models.
#
# Submit with:  sbatch --array=1-N%K alphafold_array.slurm
#   N = number of targets in AFSET_DIR
#   K = desired concurrency cap
#
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=batch_gpu
#SBATCH --constraint=csbtmp
#SBATCH --mail-user=mudit.agarwal@vanderbilt.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:nvidia_rtx_a6000:1
#SBATCH --mem=6G
#SBATCH --time=48:00:00
#SBATCH --job-name=af-array
#SBATCH --output=logs/af_array_%A_%a.out
#SBATCH --error=logs/af_array_%A_%a.err

set -euo pipefail

# ======================= USER PATHS =======================
AFSET_DIR="${AFSET_DIR:-$PWD/afset}"

# ACCRE AlphaFold 2.3.2 install
AF2_MINICONDA="${AF2_MINICONDA:-/sb/apps/alphafold232/miniconda3}"
AF2_ENV_NAME="${AF2_ENV_NAME:-af232}"
AF2_DATADIR="${AF2_DATADIR:-/sb/apps/alphafold-data.230}"
AF2_REPO="${AF2_REPO:-/sb/apps/alphafold232/alphafold}"

# Optional: Singularity image (leave empty to use conda env)
AF_SIF="${AF_SIF:-}"
AF_SIF_RUN="${AF_SIF_RUN:-/app/alphafold/run_alphafold.py}"

# AF settings
MAX_TEMPLATE_DATE="${MAX_TEMPLATE_DATE:-9999-12-31}"
USE_GPU_RELAX="${USE_GPU_RELAX:-true}"
DB_PRESET="${DB_PRESET:-full_dbs}"

mkdir -p logs

# ======================= PICK TARGET BY ARRAY INDEX =======================
mapfile -d '' -t TARGETS < <(find "$AFSET_DIR" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)
NTOTAL=${#TARGETS[@]}
: "${SLURM_ARRAY_TASK_ID:?Submit with --array=1-$NTOTAL%K}"
IDX=$((SLURM_ARRAY_TASK_ID - 1))
(( IDX >= 0 && IDX < NTOTAL )) || { echo "Index out of range ($SLURM_ARRAY_TASK_ID/$NTOTAL)"; exit 1; }

TDIR="${TARGETS[$IDX]}"
TNAME="$(basename "$TDIR")"
SEQDIR="$TDIR/sequence"
OUTROOT="$TDIR/af_out"
MONO_OUT="$OUTROOT/monomer"
MULTI_OUT="$OUTROOT/multimer"
mkdir -p "$MONO_OUT" "$MULTI_OUT"

echo "Target: $TNAME"
nvidia-smi || true
echo "Running on node(s): $SLURM_JOB_NODELIST"

# ======================= FASTA SELECTION =======================
shopt -s nullglob
FA=("$SEQDIR"/*.fa "$SEQDIR"/*.fasta "$SEQDIR"/*.faa "$SEQDIR"/*.fas "$SEQDIR"/*.fna)
shopt -u nullglob
[[ ${#FA[@]} -gt 0 ]] || { echo "[$TNAME] No FASTA found in $SEQDIR"; exit 2; }
FASTA="${FA[0]}"
FA_BASENAME="$(basename "${FASTA%.*}")"

# ======================= COMPLETION CHECKS =======================
count_ranked () {
  local preset_out="$1"
  local sub="$preset_out/$FA_BASENAME"
  [[ -d "$sub" ]] || { echo 0; return; }
  ls -1 "$sub"/ranked_*.pdb 2>/dev/null | wc -l
}

NEED_MONO=1; NEED_MULTI=1
[[ $(count_ranked "$MONO_OUT")  -ge 5 ]]  && NEED_MONO=0
[[ $(count_ranked "$MULTI_OUT") -ge 5 ]]  && NEED_MULTI=0

if [[ $NEED_MONO -eq 0 && $NEED_MULTI -eq 0 ]]; then
  echo "[$TNAME] Already complete: monomer & multimer have >= 5 ranked models. Skipping."
  exit 0
fi

# ======================= ENVIRONMENT =======================
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0

# ======================= RUNNER =======================
run_af () {
  local preset="$1" outdir="$2"
  echo "[$TNAME][$preset] FASTA=$FASTA OUT=$outdir DB=$AF2_DATADIR"

  local db_flags=(
    --uniref90_database_path="$AF2_DATADIR/uniref90/uniref90.fasta"
    --mgnify_database_path="$AF2_DATADIR/mgnify/mgy_clusters_2022_05.fa"
    --uniref30_database_path="$AF2_DATADIR/uniref30/UniRef30_2021_03"
    --bfd_database_path="$AF2_DATADIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    --pdb70_database_path="$AF2_DATADIR/pdb70/pdb70"
    --template_mmcif_dir="$AF2_DATADIR/pdb_mmcif/mmcif_files"
    --obsolete_pdbs_path="$AF2_DATADIR/pdb_mmcif/obsolete.dat"
  )

  if [[ -n "$AF_SIF" ]]; then
    singularity exec --nv \
      -B "$AF2_DATADIR":"$AF2_DATADIR" -B "$TDIR":"$TDIR" -B "$outdir":"$outdir" \
      "$AF_SIF" python "$AF_SIF_RUN" \
        --fasta_paths="$FASTA" \
        --max_template_date="$MAX_TEMPLATE_DATE" \
        --data_dir="$AF2_DATADIR" \
        --output_dir="$outdir" \
        --use_gpu_relax \
        --model_preset="$preset" \
        --db_preset="$DB_PRESET" \
        --run_relax \
        "${db_flags[@]}"
  else
    if [[ -d "$AF2_MINICONDA" ]]; then
      source "$AF2_MINICONDA/bin/activate" "$AF2_ENV_NAME"
      export LD_LIBRARY_PATH="$AF2_MINICONDA/envs/$AF2_ENV_NAME/lib:$LD_LIBRARY_PATH"
    fi

    python "$AF2_REPO/run_alphafold.py" \
        --fasta_paths="$FASTA" \
        --max_template_date="$MAX_TEMPLATE_DATE" \
        --data_dir="$AF2_DATADIR" \
        --output_dir="$outdir" \
        --use_gpu_relax \
        --model_preset="$preset" \
        --db_preset="$DB_PRESET" \
        --run_relax \
        "${db_flags[@]}"
  fi
}

# ======================= EXECUTE =======================
if [[ $NEED_MONO  -eq 1 ]]; then run_af monomer  "$MONO_OUT";  fi
if [[ $NEED_MULTI -eq 1 ]]; then run_af multimer "$MULTI_OUT"; fi

echo "[$TNAME] Done. monomer_count=$(count_ranked "$MONO_OUT") multimer_count=$(count_ranked "$MULTI_OUT")"
