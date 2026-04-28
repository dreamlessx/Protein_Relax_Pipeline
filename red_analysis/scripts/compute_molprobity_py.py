#!/usr/bin/env python3
"""
compute_molprobity_py.py - Red Analysis Pipeline

Computes MolProbity validation metrics using cctbx Python API.
Designed for SLURM array parallelization (one target per task).

Metrics extracted:
  - Clashscore (steric clashes per 1000 atoms)
  - Ramachandran outliers (%)
  - Ramachandran favored (%)
  - Rotamer outliers (%)
  - MolProbity score (composite)
  - C-beta deviations
  - RMS bonds / RMS angles (from geometry restraints)

Uses individual validators (ramalyze, rotalyze, cbetadev, clashscore)
instead of the composite molprobity.molprobity() which requires cablam data
we don't have on this cluster.

MolProbity composite score formula (Chen et al. 2010):
  MP = 0.426*ln(1+clashscore) + 0.33*ln(1+max(0,rota_pct-1)) + 0.25*ln(1+max(0,100-rama_fav_pct-2))
  (plus 0.5 if any cbeta outliers, but we omit that correction)

COMMENT: MolProbity is the key complement to TM-score. TM-score tells us
"how close to the right answer?" MolProbity tells us "how physically
reasonable is this structure?" A structure can have perfect TM-score but
terrible MolProbity (e.g., steric clashes, bad rotamers). The hypothesis
is that AMBER relaxation fixes MolProbity metrics without changing global
fold (TM-score) - which is exactly what we saw in Phase 1.

If that's true, it would mean:
  - AMBER is doing what it should: optimizing local geometry
  - TM-score is insensitive to what AMBER fixes
  - MolProbity should show clear improvement after AMBER

This is actually the INTERESTING finding. If AMBER doesn't help MolProbity
either, then what's the point of AMBER relaxation?
"""

import os
import sys
import math
import tempfile
import warnings
warnings.filterwarnings('ignore')

# === Environment setup ===
# Set monomer library path (ccp4_mon_lib in conda, underscore not hyphen)
conda_prefix = os.environ.get('CONDA_PREFIX', '')
if not conda_prefix:
    # Fallback: derive from python executable path
    conda_prefix = os.path.dirname(os.path.dirname(sys.executable))
mon_lib_path = os.path.join(conda_prefix, 'share', 'ccp4_mon_lib')
if os.path.exists(mon_lib_path):
    os.environ['MMTBX_CCP4_MONOMER_LIB'] = mon_lib_path

# === Monkey-patch 1: SS disulfide link alias ===
# CCP4 monomer library uses 'disulf' key; cctbx expects 'SS'
import mmtbx.monomer_library.server as _mls
_orig_server_init = _mls.server.__init__
def _patched_server_init(self, *args, **kwargs):
    _orig_server_init(self, *args, **kwargs)
    if "SS" not in self.link_link_id_dict and "disulf" in self.link_link_id_dict:
        self.link_link_id_dict["SS"] = self.link_link_id_dict["disulf"]
_mls.server.__init__ = _patched_server_init

# === Monkey-patch 2: probe/reduce module detection ===
# cctbx checks libtbx.env.has_module("probe") but we installed standalone
import libtbx.env_config
_orig_has_module = libtbx.env_config.environment.has_module
def _patched_has_module(self, name):
    if name in ("probe", "reduce"):
        return True
    return _orig_has_module(self, name)
libtbx.env_config.environment.has_module = _patched_has_module

import iotbx.pdb
import mmtbx.model
from mmtbx.validation import ramalyze, rotalyze, cbetadev, clashscore

BLUE_BASE = "/data/p_csb_meiler/agarwm5/af_work"
GREEN_BASE = "/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/data"
CRYSTAL_DIR = "/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/cleaned"
TARGET_LIST = "/data/p_csb_meiler/agarwm5/red_analysis/target_list.txt"
OUTDIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"


def compute_mp_score(clash_val, rota_pct, rama_fav_pct):
    """Compute MolProbity composite score (Chen et al. 2010).

    Formula: 0.426*ln(1+clashscore) + 0.33*ln(1+max(0,rota%-1)) + 0.25*ln(1+max(0,100-ramafav%-2))
    """
    try:
        score = (0.426 * math.log(1 + float(clash_val))
                 + 0.33 * math.log(1 + max(0, float(rota_pct) - 1))
                 + 0.25 * math.log(1 + max(0, 100 - float(rama_fav_pct) - 2)))
        return score
    except (ValueError, TypeError):
        return float('nan')


def run_molprobity(pdb_path, target, pipeline, source, model_idx, outfile):
    """Run MolProbity individual validators on a single PDB and append results."""
    if not os.path.exists(pdb_path):
        print(f"  WARN: Missing {pdb_path}", file=sys.stderr)
        return

    tmppath = None
    try:
        # Strip hydrogen atoms - AF-relaxed and AMBER structures contain H atoms
        # that cause "unknown nonbonded energy type symbols" errors in cctbx.
        # MolProbity uses reduce to place its own hydrogens for clashscore anyway.
        pdb_inp = iotbx.pdb.input(file_name=pdb_path)
        hierarchy = pdb_inp.construct_hierarchy()
        h_sel = hierarchy.atom_selection_cache().selection("not (element H or element D)")
        hierarchy_stripped = hierarchy.select(h_sel)

        # Write stripped PDB to temp file and re-read (cctbx needs clean input)
        with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as tmp:
            tmp.write(hierarchy_stripped.as_pdb_string())
            tmppath = tmp.name

        pdb_inp = iotbx.pdb.input(file_name=tmppath)
        model = mmtbx.model.manager(model_input=pdb_inp)
        model.process(make_restraints=True)
        hierarchy = model.get_hierarchy()

        # Ramachandran analysis
        rama = ramalyze.ramalyze(pdb_hierarchy=hierarchy, outliers_only=False)
        rama_out = rama.percent_outliers
        rama_fav = rama.percent_favored

        # Rotamer analysis
        rota = rotalyze.rotalyze(pdb_hierarchy=hierarchy, outliers_only=False)
        rota_out = rota.percent_outliers

        # C-beta deviations
        cbeta = cbetadev.cbetadev(pdb_hierarchy=hierarchy, outliers_only=True)
        cbeta_count = cbeta.get_outlier_count()

        # Clashscore (requires probe binary)
        try:
            clash = clashscore.clashscore(pdb_hierarchy=hierarchy)
            clash_val = clash.get_clashscore()
        except Exception as e:
            print(f"  WARN: Clashscore failed for {pdb_path}: {e}", file=sys.stderr)
            clash_val = "NA"

        # MolProbity composite score
        mp_score = compute_mp_score(clash_val, rota_out, rama_fav)

        # Geometry from restraints (RMS bonds and angles)
        try:
            geom = model.get_restraints_manager().geometry
            sites_cart = model.get_sites_cart()
            energies = geom.energies_sites(sites_cart=sites_cart, compute_gradients=False)
            rms_b = f"{energies.bond_deviations()[2]:.4f}"
            rms_a = f"{energies.angle_deviations()[2]:.4f}"
        except Exception:
            rms_b = "NA"
            rms_a = "NA"

        # Format clashscore
        clash_str = f"{clash_val:.2f}" if isinstance(clash_val, (int, float)) else clash_val
        mp_str = f"{mp_score:.2f}" if not math.isnan(mp_score) else "NA"

        with open(outfile, 'a') as f:
            f.write(f"{target}\t{pipeline}\t{source}\t{model_idx}\t"
                    f"{clash_str}\t{rama_out:.2f}\t{rama_fav:.2f}\t"
                    f"{rota_out:.2f}\t{mp_str}\t{cbeta_count}\t{rms_b}\t{rms_a}\n")

    except Exception as e:
        print(f"  ERROR on {pdb_path}: {e}", file=sys.stderr)
    finally:
        if tmppath and os.path.exists(tmppath):
            os.unlink(tmppath)


def main():
    pipeline = sys.argv[1] if len(sys.argv) > 1 else 'blue'

    # Get target from SLURM or CLI
    slurm_id = os.environ.get('SLURM_ARRAY_TASK_ID')
    if slurm_id:
        with open(TARGET_LIST) as f:
            targets = f.read().strip().split('\n')
        target = targets[int(slurm_id) - 1]
    else:
        target = sys.argv[2] if len(sys.argv) > 2 else '1A2K'

    outfile = os.path.join(OUTDIR, f"molprobity_{pipeline}_{target}.tsv")

    # Write header if new file
    if not os.path.exists(outfile):
        with open(outfile, 'w') as f:
            f.write("target\tpipeline\tsource\tmodel_idx\tclashscore\t"
                    "rama_outliers\trama_favored\trota_outliers\t"
                    "molprobity_score\tcbeta_outliers\trms_bonds\trms_angles\n")

    print(f"=== Red MolProbity: {target} ({pipeline}) ===")

    # Crystal structure baseline
    crystal_pdb = os.path.join(CRYSTAL_DIR, f"{target}.pdb")
    run_molprobity(crystal_pdb, target, pipeline, "crystal", "0", outfile)

    if pipeline == 'blue':
        base = os.path.join(BLUE_BASE, target)
        for i in range(5):
            run_molprobity(f"{base}/af_out/ranked_{i}.pdb", target, pipeline, "af_relaxed", str(i), outfile)
        for i in range(5):
            pdb = f"{base}/af_out_unrelaxed/ranked_{i}.pdb"
            if not os.path.exists(pdb):
                pdb = f"{base}/af_out/unrelaxed_model_{i}_ptm.pdb"
            run_molprobity(pdb, target, pipeline, "af_unrelaxed", str(i), outfile)
        for i in range(5):
            run_molprobity(f"{base}/boltz_out/boltz_results_boltz_input/predictions/boltz_input/boltz_input_model_{i}.pdb",
                           target, pipeline, "boltz", str(i), outfile)
        for i in range(5):
            run_molprobity(f"{base}/amber_out/af_unrelaxed_ranked_{i}/relaxed.pdb",
                           target, pipeline, "amber_af", str(i), outfile)
        for i in range(5):
            run_molprobity(f"{base}/amber_out/boltz_model_{i}/relaxed.pdb",
                           target, pipeline, "amber_boltz", str(i), outfile)

    elif pipeline == 'green':
        base = os.path.join(GREEN_BASE, target)
        for i in range(5):
            run_molprobity(f"{base}/AF/ranked_{i}.pdb", target, pipeline, "af_relaxed", str(i), outfile)
        for i in range(5):
            run_molprobity(f"{base}/af_out_unrelaxed/ranked_{i}.pdb",
                           target, pipeline, "af_unrelaxed", str(i), outfile)
        for i in range(5):
            run_molprobity(f"{base}/Boltz/boltz_input_model_{i}.pdb",
                           target, pipeline, "boltz", str(i), outfile)
        for i in range(5):
            run_molprobity(f"{base}/amber_out/af_unrelaxed_{i}/relaxed.pdb",
                           target, pipeline, "amber_af", str(i), outfile)
        for i in range(5):
            run_molprobity(f"{base}/amber_out/boltz_model_{i}/relaxed.pdb",
                           target, pipeline, "amber_boltz", str(i), outfile)

    print(f"=== Done: {outfile} ===")


if __name__ == '__main__':
    main()
