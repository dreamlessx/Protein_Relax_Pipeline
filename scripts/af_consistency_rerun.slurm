#!/bin/bash --norc
#SBATCH --constraint=csbtmp
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=batch_gpu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --gres=gpu:nvidia_rtx_a6000:1
#SBATCH --mem=80G
#SBATCH --time=72:00:00
#SBATCH --job-name=af_consist
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/af_consistency_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/af_consistency_%A_%a.err

# AlphaFold 2.3.2 re-prediction for 236 targets with crystal-derived FASTAs
# All FASTAs derived from crystal PDB coordinates for exact consistency
# Submit with: sbatch --array=1-236%10 af_consistency_rerun.slurm

set -uo pipefail

DIRLIST=/data/p_csb_meiler/agarwm5/protein_pipeline/af_crystal_derived_dirlist.txt
AF2_MINICONDA=/sb/apps/alphafold232/miniconda3
AF2_DATADIR=/csbtmp/alphafold-data.230
AF2_REPO=/sb/apps/alphafold232/alphafold

CALCDIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$DIRLIST")

if [[ -z "$CALCDIR" ]]; then
    echo "ERROR: Could not read directory for task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

TARGET=$(basename "$CALCDIR")
cd "$CALCDIR" || exit 1

echo "Task $SLURM_ARRAY_TASK_ID: $TARGET (consistency fix re-run)"
echo "Node: $SLURM_JOB_NODELIST"

# Archive old predictions
if [[ -d af_out ]]; then
    mv af_out af_out_old_$(date +%Y%m%d) 2>/dev/null || true
fi
if [[ -d af_out_unrelaxed ]]; then
    mv af_out_unrelaxed af_out_unrelaxed_old_$(date +%Y%m%d) 2>/dev/null || true
fi

FASTA=sequence.fasta
[[ -f "$FASTA" ]] || { echo "ERROR: No sequence.fasta"; exit 2; }

SEQCOUNT=$(grep -c '^>' "$FASTA")
echo "FASTA: $FASTA ($SEQCOUNT chains)"

source $AF2_MINICONDA/bin/activate af232
export LD_LIBRARY_PATH=$AF2_MINICONDA/envs/af232/lib:${LD_LIBRARY_PATH:-}

run_af() {
    local DB_MODE="$1"
    echo "=== Running AF with $DB_MODE ==="

    local COMMON_ARGS=(
        --fasta_paths=$FASTA
        --max_template_date=9999-12-31
        --data_dir=$AF2_DATADIR
        --output_dir=$CALCDIR/af_out
        --use_gpu_relax
        --models_to_relax=all
        --uniref90_database_path=$AF2_DATADIR/uniref90/uniref90.fasta
        --mgnify_database_path=$AF2_DATADIR/mgnify/mgy_clusters_2022_05.fa
        --template_mmcif_dir=$AF2_DATADIR/pdb_mmcif/mmcif_files
        --obsolete_pdbs_path=$AF2_DATADIR/pdb_mmcif/obsolete.dat
    )

    local DB_ARGS=()
    if [[ "$DB_MODE" == "full_dbs" ]]; then
        DB_ARGS=(
            --uniref30_database_path=$AF2_DATADIR/uniref30/UniRef30_2021_03
            --bfd_database_path=$AF2_DATADIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
        )
    else
        DB_ARGS=(
            --db_preset=reduced_dbs
            --small_bfd_database_path=$AF2_DATADIR/small_bfd/bfd-first_non_consensus_sequences.fasta
        )
    fi

    if (( SEQCOUNT > 1 )); then
        python $AF2_REPO/run_alphafold.py \
            "${COMMON_ARGS[@]}" \
            "${DB_ARGS[@]}" \
            --model_preset=multimer \
            --pdb_seqres_database_path=$AF2_DATADIR/pdb_seqres/pdb_seqres.txt \
            --uniprot_database_path=$AF2_DATADIR/uniprot/uniprot.fasta \
            --num_multimer_predictions_per_model=1
    else
        python $AF2_REPO/run_alphafold.py \
            "${COMMON_ARGS[@]}" \
            "${DB_ARGS[@]}" \
            --pdb70_database_path=$AF2_DATADIR/pdb70/pdb70
    fi
}

# Try full_dbs first, fallback to reduced_dbs
if run_af "full_dbs"; then
    echo "full_dbs succeeded for $TARGET"
else
    echo "full_dbs FAILED for $TARGET"
    unrelaxed_count=$(find "$CALCDIR/af_out" -name "unrelaxed_model_*.pdb" 2>/dev/null | wc -l)
    if [[ "$unrelaxed_count" -gt 0 ]]; then
        echo "WARNING: Found $unrelaxed_count unrelaxed models — AMBER likely failed."
        echo "Keeping unrelaxed output."
    else
        echo "No unrelaxed models — MSA failure. Trying reduced_dbs..."
        rm -rf "$CALCDIR/af_out" 2>/dev/null
        if run_af "reduced_dbs"; then
            echo "reduced_dbs succeeded for $TARGET"
        else
            echo "ERROR: Both full_dbs and reduced_dbs failed for $TARGET"
            exit 1
        fi
    fi
fi

# Extract ranked PDBs and organize output
SEQDIR=$(find $CALCDIR/af_out -maxdepth 1 -type d ! -name af_out | head -1)
if [[ -n "$SEQDIR" ]]; then
    # Save unrelaxed models
    mkdir -p $CALCDIR/af_out_unrelaxed
    for i in 1 2 3 4 5; do
        src=$(ls $SEQDIR/unrelaxed_model_${i}_*.pdb 2>/dev/null | head -1)
        [[ -f "$src" ]] && cp "$src" "$CALCDIR/af_out_unrelaxed/ranked_$((i-1)).pdb"
    done
    # Move ranked (AMBER-relaxed) PDBs up
    for f in $SEQDIR/ranked_*.pdb; do
        [[ -f "$f" ]] && mv "$f" $CALCDIR/af_out/
    done
fi

# Cleanup intermediates
for subdir in $CALCDIR/af_out/*/; do
    rm -rf "$subdir/msas" 2>/dev/null
    rm -f "$subdir"/features.pkl 2>/dev/null
    rm -f "$subdir"/result_model_*.pkl 2>/dev/null
    rm -f "$subdir"/relaxed_model_*.pdb 2>/dev/null
    rm -f "$subdir"/timings.json 2>/dev/null
done

echo "Done: $TARGET (relaxed + unrelaxed)"
ls $CALCDIR/af_out/ranked_*.pdb 2>/dev/null | wc -l
ls $CALCDIR/af_out_unrelaxed/ranked_*.pdb 2>/dev/null | wc -l
