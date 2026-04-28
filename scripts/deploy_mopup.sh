#!/bin/bash
# Deploy and submit mop-up jobs for both Blue and Green pipelines
# Run from local machine: bash scripts/deploy_mopup.sh
#
# What this does:
#   1. Copies scan + SLURM scripts to ACCRE
#   2. Generates worklists (scans for all missing runs)
#   3. Patches array size into SLURM scripts
#   4. Submits both mop-up arrays
#
# Safe to run while existing jobs are active - all scripts have skip-logic

set -euo pipefail

ACCRE="accre"
BLUE_SCRIPTS="/data/p_csb_meiler/agarwm5/protein_pipeline/scripts"
BLUE_BASE="/data/p_csb_meiler/agarwm5/protein_pipeline"
GREEN_BASE="/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking"

echo "=== Deploying mop-up scripts to ACCRE ==="

# Copy scripts
scp scripts/blue_scan_gaps.sh scripts/blue_mopup.slurm "${ACCRE}:${BLUE_SCRIPTS}/"
scp scripts/green_scan_gaps.sh scripts/green_mopup.slurm "${ACCRE}:${GREEN_BASE}/"

echo "=== Generating Blue worklist (this takes ~5 min) ==="
ssh "$ACCRE" "cd ${BLUE_SCRIPTS} && bash blue_scan_gaps.sh > ${BLUE_BASE}/blue_worklist.tsv && wc -l ${BLUE_BASE}/blue_worklist.tsv"

echo "=== Generating Green worklist (this takes ~5 min) ==="
ssh "$ACCRE" "cd ${GREEN_BASE} && bash green_scan_gaps.sh > ${GREEN_BASE}/green_worklist.tsv && wc -l ${GREEN_BASE}/green_worklist.tsv"

echo "=== Patching array sizes and submitting ==="
ssh "$ACCRE" bash <<'REMOTE'
set -euo pipefail

BLUE_BASE="/data/p_csb_meiler/agarwm5/protein_pipeline"
GREEN_BASE="/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking"

# Blue
BLUE_N=$(wc -l < "${BLUE_BASE}/blue_worklist.tsv")
if [[ "$BLUE_N" -gt 0 ]]; then
  sed -i "s/^# SBATCH --array=1-NLINES%100/#SBATCH --array=1-${BLUE_N}%100/" "${BLUE_BASE}/scripts/blue_mopup.slurm"
  echo "Blue: $BLUE_N work items"
  cd "${BLUE_BASE}" && sbatch scripts/blue_mopup.slurm
else
  echo "Blue: nothing to do!"
fi

# Green
GREEN_N=$(wc -l < "${GREEN_BASE}/green_worklist.tsv")
if [[ "$GREEN_N" -gt 0 ]]; then
  sed -i "s/^# SBATCH --array=1-NLINES%100/#SBATCH --array=1-${GREEN_N}%100/" "${GREEN_BASE}/green_mopup.slurm"
  echo "Green: $GREEN_N work items"
  cd "${GREEN_BASE}" && sbatch green_mopup.slurm
else
  echo "Green: nothing to do!"
fi
REMOTE

echo "=== Done. Check with: ssh accre squeue -u \$(whoami) ==="
