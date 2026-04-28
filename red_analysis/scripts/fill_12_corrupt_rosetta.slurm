#!/bin/bash
#SBATCH --job-name=fill12_ros
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/fill12_ros_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/fill12_ros_%A_%a.err
#SBATCH --exclude=cn1340
#SBATCH --array=0-11

set -eu

RELAX=/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/source/bin/relax.linuxgccrelease
DB=/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/database
BLUE=/data/p_csb_meiler/agarwm5/af_work
GREEN=/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/data

# 12 missing rows: target pipeline src_type protocol rep input_pdb_basename
declare -a ENTRIES=(
  "1AZS blue boltz_boltz_input_model_0 dualspace_beta  2 boltz_input_model_0 boltz_out/boltz_results_boltz_input/predictions/boltz_input/boltz_input_model_0.pdb"
  "1BKD blue boltz_boltz_input_model_1 dualspace_ref15 2 boltz_input_model_1 boltz_out/boltz_results_boltz_input/predictions/boltz_input/boltz_input_model_1.pdb"
  "1E6E blue boltz_boltz_input_model_3 dualspace_beta  5 boltz_input_model_3 boltz_out/boltz_results_boltz_input/predictions/boltz_input/boltz_input_model_3.pdb"
  "1E6J blue boltz_boltz_input_model_2 normal_ref15    4 boltz_input_model_2 boltz_out/boltz_results_boltz_input/predictions/boltz_input/boltz_input_model_2.pdb"
  "1GLA blue boltz_boltz_input_model_1 cartesian_beta  4 boltz_input_model_1 boltz_out/boltz_results_boltz_input/predictions/boltz_input/boltz_input_model_1.pdb"
  "1GP2 blue boltz_boltz_input_model_0 cartesian_ref15 2 boltz_input_model_0 boltz_out/boltz_results_boltz_input/predictions/boltz_input/boltz_input_model_0.pdb"
  "1GXD blue boltz_boltz_input_model_1 dualspace_ref15 4 boltz_input_model_1 boltz_out/boltz_results_boltz_input/predictions/boltz_input/boltz_input_model_1.pdb"
  "1H1V blue boltz_boltz_input_model_3 cartesian_ref15 3 boltz_input_model_3 boltz_out/boltz_results_boltz_input/predictions/boltz_input/boltz_input_model_3.pdb"
  "1HE8 blue af_unrelaxed_ranked_1 normal_ref15 5 ranked_1 af_out_unrelaxed/ranked_1.pdb"
  "2Z0E blue boltz_boltz_input_model_4 dualspace_ref15 2 boltz_input_model_4 boltz_out/boltz_results_boltz_input/predictions/boltz_input/boltz_input_model_4.pdb"
  "4DN4 green boltz_boltz_input_model_4 normal_beta 4 boltz_input_model_4 boltz_out/boltz_results_boltz_input/predictions/boltz_input/boltz_input_model_4.pdb"
  "6AL0 blue af_unrelaxed_ranked_0 cartesian_beta 4 ranked_0 af_out_unrelaxed/ranked_0.pdb"
)

E=${ENTRIES[$SLURM_ARRAY_TASK_ID]}
read T P SRC PROTO REP PREFIX REL_INPUT <<< "$E"

if [ "$P" = "blue" ]; then ROOT=$BLUE; else ROOT=$GREEN; fi
INPUT=$ROOT/$T/$REL_INPUT
OUTDIR=$ROOT/$T/rosetta_out/$SRC/$PROTO
OUTPATH=$OUTDIR/${PREFIX}_r${REP}.pdb.gz

echo "Task $SLURM_ARRAY_TASK_ID: $T $P $SRC $PROTO r$REP"
echo "Input: $INPUT"
echo "Output: $OUTPATH"
echo "Host: $(hostname); Start: $(date)"

test -f "$INPUT" || { echo "MISSING INPUT PDB"; exit 1; }
mkdir -p "$OUTDIR"

# Remove corrupt file if present
if [ -f "$OUTPATH" ]; then
  size=$(stat -c %s "$OUTPATH")
  if ! gunzip -t "$OUTPATH" 2>/dev/null; then
    echo "Removing corrupt $OUTPATH (size=$size)"
    rm -f "$OUTPATH"
  fi
fi

declare -A CFG
CFG[cartesian_beta]="-relax:cartesian -beta_nov16 -score:weights beta_nov16_cart"
CFG[cartesian_ref15]="-relax:cartesian -score:weights ref2015_cart"
CFG[dualspace_beta]="-relax:dualspace -beta_nov16 -score:weights beta_nov16_cart -nonideal -relax:minimize_bond_angles -relax:minimize_bond_lengths"
CFG[dualspace_ref15]="-relax:dualspace -score:weights ref2015_cart -nonideal -relax:minimize_bond_angles -relax:minimize_bond_lengths"
CFG[normal_beta]="-beta_nov16 -score:weights beta_nov16"
CFG[normal_ref15]="-score:weights ref2015"

cd "$OUTDIR"

$RELAX \
  -in:file:s "$INPUT" \
  -database "$DB" \
  ${CFG[$PROTO]} \
  -ignore_zero_occupancy false \
  -nstruct 1 -no_nstruct_label -out:pdb_gz -flip_HNQ \
  -fa_max_dis 9.0 -optimization:default_max_cycles 200 \
  -out:levels all:warning \
  -suffix "_fill${REP}" 2>&1 | tail -5

BASENAME=$(basename "$INPUT" .pdb)
CANDIDATE="${BASENAME}_fill${REP}.pdb.gz"
if [ -f "$CANDIDATE" ]; then
  mv "$CANDIDATE" "${PREFIX}_r${REP}.pdb.gz"
  echo "wrote ${PREFIX}_r${REP}.pdb.gz"
  gunzip -t "${PREFIX}_r${REP}.pdb.gz" && echo "gzip integrity OK"
else
  echo "RELAX DID NOT WRITE OUTPUT"
fi
echo "End: $(date)"
