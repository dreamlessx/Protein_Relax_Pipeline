#!/usr/bin/env python3
"""
Fill the 12 exact missing MolProbity rows identified in the audit.
Each row: (target, pipeline, src_type, protocol, rep, pdb_path).
Appends one line per target TSV using the schema of compute_rosetta_molprobity.py.
"""
import os, sys, gzip, tempfile, warnings
warnings.filterwarnings('ignore')

# Environment setup
conda_prefix = os.environ.get('CONDA_PREFIX', '')
conda_bin = os.path.join(conda_prefix, 'bin')
if conda_bin not in os.environ.get('PATH', ''):
    os.environ['PATH'] = conda_bin + ':' + os.environ.get('PATH', '')

from mmtbx.validation.molprobity import molprobity

BLUE_BASE = "/data/p_csb_meiler/agarwm5/af_work"
GREEN_BASE = "/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/data"
METRICS = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"

# Each entry: (target, pipeline, src_type, protocol, rep, pdb_gz_basename_pattern)
MISSING = [
    ("1AZS", "blue",  "boltz_boltz_input_model_0", "dualspace_beta",  2, "boltz_input_model_0_r2.pdb.gz"),
    ("1BKD", "blue",  "boltz_boltz_input_model_1", "dualspace_ref15", 2, "boltz_input_model_1_r2.pdb.gz"),
    ("1E6E", "blue",  "boltz_boltz_input_model_3", "dualspace_beta",  5, "boltz_input_model_3_r5.pdb.gz"),
    ("1E6J", "blue",  "boltz_boltz_input_model_2", "normal_ref15",    4, "boltz_input_model_2_r4.pdb.gz"),
    ("1GLA", "blue",  "boltz_boltz_input_model_1", "cartesian_beta",  4, "boltz_input_model_1_r4.pdb.gz"),
    ("1GP2", "blue",  "boltz_boltz_input_model_0", "cartesian_ref15", 2, "boltz_input_model_0_r2.pdb.gz"),
    ("1GXD", "blue",  "boltz_boltz_input_model_1", "dualspace_ref15", 4, "boltz_input_model_1_r4.pdb.gz"),
    ("1H1V", "blue",  "boltz_boltz_input_model_3", "cartesian_ref15", 3, "boltz_input_model_3_r3.pdb.gz"),
    ("1HE8", "blue",  "af_unrelaxed_ranked_1",     "normal_ref15",    5, "ranked_1_r5.pdb.gz"),
    ("2Z0E", "blue",  "boltz_boltz_input_model_4", "dualspace_ref15", 2, "boltz_input_model_4_r2.pdb.gz"),
    ("4DN4", "green", "boltz_boltz_input_model_4", "normal_beta",     4, "boltz_input_model_4_r4.pdb.gz"),
    ("6AL0", "blue",  "af_unrelaxed_ranked_0",     "cartesian_beta",  4, "ranked_0_r4.pdb.gz"),
]

def score_pdb_gz(pdb_gz_path):
    if not os.path.exists(pdb_gz_path):
        return None
    with tempfile.NamedTemporaryFile(mode='wb', suffix='.pdb', delete=False) as t:
        with gzip.open(pdb_gz_path, 'rb') as g:
            t.write(g.read())
        tmp = t.name
    try:
        m = molprobity.molprobity(pdb_file=tmp, keep_hydrogens=False,
                                   outliers_only=False, quiet=True)
        rama = m.ramalyze
        rota = m.rotalyze
        cb = m.cbetadev
        clash_str = f"{m.clashes.get_clashscore():.2f}" if m.clashes else "nan"
        mp_str = f"{m.molprobity_score():.2f}" if m.molprobity_score() is not None else "nan"
        rb = m.restraints.bonds if m.restraints else None
        ra = m.restraints.angles if m.restraints else None
        rms_b = f"{rb.mean:.4f}" if rb else "nan"
        rms_a = f"{ra.mean:.4f}" if ra else "nan"
        return (clash_str, f"{rama.percent_outliers:.2f}", f"{rama.percent_favored:.2f}",
                f"{rota.percent_outliers:.2f}", mp_str, cb.get_outlier_count(), rms_b, rms_a)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def resolve_pdb_path(target, pipeline, src_type, protocol, fname):
    base = BLUE_BASE if pipeline == "blue" else GREEN_BASE
    return os.path.join(base, target, "rosetta_out", src_type, protocol, fname)


def main():
    for target, pipeline, src_type, protocol, rep, fname in MISSING:
        pdb = resolve_pdb_path(target, pipeline, src_type, protocol, fname)
        if not os.path.exists(pdb):
            print(f"MISSING PDB: {pdb}")
            continue
        out = os.path.join(METRICS, f"rosetta_molprobity_{pipeline}_{target}.tsv")
        if not os.path.exists(out):
            print(f"MISSING TSV: {out}")
            continue
        print(f"scoring {target} {pipeline} {src_type} {protocol} r{rep}...")
        res = score_pdb_gz(pdb)
        if res is None:
            print("  score failed")
            continue
        clash, ro, rf, rro, mp, cb, rb, ra = res
        line = f"{target}\t{pipeline}\t{src_type}\t{protocol}\t{rep}\t{clash}\t{ro}\t{rf}\t{rro}\t{mp}\t{cb}\t{rb}\t{ra}\n"
        with open(out, 'a') as f:
            f.write(line)
        print(f"  appended: {line.strip()}")


if __name__ == "__main__":
    main()
