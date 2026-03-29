Check Rosetta relaxation progress on ACCRE.
1. SSH to ACCRE, count completed .pdb files in relaxation output dirs
2. Compare against expected total (~200,000 = 257 targets x 6 protocols x 5 replicates x ~26 input types)
3. Report: completed/total per protocol, per input type
4. Check squeue for running/pending jobs
5. Identify any failed jobs (check slurm output files for errors)
6. Estimate time to completion based on current throughput
