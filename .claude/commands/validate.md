Run MolProbity validation on completed structures.
1. SSH to Curium or ACCRE
2. Activate protein conda env
3. Run run_validation_parallel.py with 12 workers on completed relaxation outputs
4. Run molprobity_extended.py for RMSZ scores
5. Run posebusters.py for geometry checks
6. Parse results, report per-protocol averages
7. Flag any structures with MolProbity score > 2.0 or clashscore > 20
