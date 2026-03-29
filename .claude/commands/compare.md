Generate AF2 vs Boltz1 comparison analysis.
1. Load validation results from molprobity_full.csv
2. Group by prediction method (AF2 relaxed, AF2 unrelaxed, Boltz, crystal)
3. Compute mean +/- std for: Ramachandran favored, clashscore, MolProbity score, rotamer outliers
4. Statistical tests: paired t-test or Wilcoxon between methods
5. Generate comparison table (markdown + CSV)
6. Generate bar charts with error bars
7. Save to validation_results/comparison_report/
