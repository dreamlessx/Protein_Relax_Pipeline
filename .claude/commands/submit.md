Submit SLURM jobs for remaining Rosetta relaxation runs.
1. SSH to ACCRE
2. Identify which targets/protocols still need runs
3. Generate SLURM array job script for remaining work
4. Use: --account=p_meiler_acc, --partition=batch, --mem=8G, --time=72:00:00
5. Submit with sbatch, return job ID
6. Confirm with squeue
