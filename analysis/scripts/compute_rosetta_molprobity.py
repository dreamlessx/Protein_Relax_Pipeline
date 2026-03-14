#!/usr/bin/env python3
"""
compute_rosetta_molprobity.py — Red Analysis Pipeline

Computes MolProbity metrics for Rosetta-relaxed structures.
Handles .pdb.gz files (decompresses on the fly).

COMMENT: This is Phase 4 — the most important analysis for the paper.
The key question: does Rosetta improve MolProbity metrics?
We already know:
  - Rosetta degrades TM-score by ~0.02 (Phase 2)
  - AMBER improves MolProbity without changing TM-score (Phase 3)

If Rosetta ALSO improves MolProbity, then Rosetta is doing something
meaningful: trading a small amount of global accuracy for better local
geometry. If Rosetta doesn't improve MolProbity either, then what's
the point?

Expected: Rosetta should dramatically improve clashscore and rotamers
(that's what its energy function optimizes). The interesting question
is which protocol improves MolProbity the MOST while degrading TM-score
the LEAST.
"""

import os
import sys
import math
import gzip
import glob
import tempfile
import warnings
warnings.filterwarnings('ignore')

# Environment setup
conda_prefix = os.environ.get('CONDA_PREFIX', '')
if not conda_prefix:
    conda_prefix = os.path.dirname(os.path.dirname(sys.executable))
mon_lib_path = os.path.join(conda_prefix, 'share', 'ccp4_mon_lib')
if os.path.exists(mon_lib_path):
    os.environ['MMTBX_CCP4_MONOMER_LIB'] = mon_lib_path

# Monkey-patches
import mmtbx.monomer_library.server as _mls
_orig_server_init = _mls.server.__init__
def _patched_server_init(self, *args, **kwargs):
    _orig_server_init(self, *args, **kwargs)
    if "SS" not in self.link_link_id_dict and "disulf" in self.link_link_id_dict:
        self.link_link_id_dict["SS"] = self.link_link_id_dict["disulf"]
_mls.server.__init__ = _patched_server_init

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
TARGET_LIST = "/data/p_csb_meiler/agarwm5/red_analysis/target_list.txt"
OUTDIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"

PROTOCOLS = ["cartesian_beta", "cartesian_ref15", "dualspace_beta",
             "dualspace_ref15", "normal_beta", "normal_ref15"]


def compute_mp_score(clash_val, rota_pct, rama_fav_pct):
    """MolProbity composite score (Chen et al. 2010)."""
    try:
        return (0.426 * math.log(1 + float(clash_val))
                + 0.33 * math.log(1 + max(0, float(rota_pct) - 1))
                + 0.25 * math.log(1 + max(0, 100 - float(rama_fav_pct) - 2)))
    except (ValueError, TypeError):
        return float('nan')


def run_molprobity_gz(pdbgz_path, target, pipeline, src_type, protocol, rep, outfile):
    """Run MolProbity on a .pdb.gz file."""
    if not os.path.exists(pdbgz_path):
        return

    tmppath = None
    tmpgz = None
    try:
        # Decompress .pdb.gz to temp file
        with gzip.open(pdbgz_path, 'rt') as gz:
            pdb_content = gz.read()

        # Strip H atoms and write to temp file
        with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as tmp:
            tmpgz = tmp.name
            tmp.write(pdb_content)

        pdb_inp = iotbx.pdb.input(file_name=tmpgz)
        hierarchy = pdb_inp.construct_hierarchy()
        h_sel = hierarchy.atom_selection_cache().selection("not (element H or element D)")
        hierarchy_stripped = hierarchy.select(h_sel)

        with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as tmp:
            tmp.write(hierarchy_stripped.as_pdb_string())
            tmppath = tmp.name

        pdb_inp2 = iotbx.pdb.input(file_name=tmppath)
        model = mmtbx.model.manager(model_input=pdb_inp2)
        model.process(make_restraints=True)
        h = model.get_hierarchy()

        # Validators
        rama = ramalyze.ramalyze(pdb_hierarchy=h, outliers_only=False)
        rota = rotalyze.rotalyze(pdb_hierarchy=h, outliers_only=False)
        cb = cbetadev.cbetadev(pdb_hierarchy=h, outliers_only=True)

        try:
            clash = clashscore.clashscore(pdb_hierarchy=h)
            clash_val = clash.get_clashscore()
        except:
            clash_val = "NA"

        mp_score = compute_mp_score(clash_val, rota.percent_outliers, rama.percent_favored)

        try:
            geom = model.get_restraints_manager().geometry
            sites_cart = model.get_sites_cart()
            energies = geom.energies_sites(sites_cart=sites_cart, compute_gradients=False)
            rms_b = f"{energies.bond_deviations()[2]:.4f}"
            rms_a = f"{energies.angle_deviations()[2]:.4f}"
        except:
            rms_b = "NA"
            rms_a = "NA"

        clash_str = f"{clash_val:.2f}" if isinstance(clash_val, (int, float)) else clash_val
        mp_str = f"{mp_score:.2f}" if not math.isnan(mp_score) else "NA"

        with open(outfile, 'a') as f:
            f.write(f"{target}\t{pipeline}\t{src_type}\t{protocol}\t{rep}\t"
                    f"{clash_str}\t{rama.percent_outliers:.2f}\t{rama.percent_favored:.2f}\t"
                    f"{rota.percent_outliers:.2f}\t{mp_str}\t{cb.get_outlier_count()}\t"
                    f"{rms_b}\t{rms_a}\n")

    except Exception as e:
        print(f"  ERROR on {pdbgz_path}: {e}", file=sys.stderr)
    finally:
        for p in [tmppath, tmpgz]:
            if p and os.path.exists(p):
                os.unlink(p)


def main():
    pipeline = sys.argv[1] if len(sys.argv) > 1 else 'blue'

    slurm_id = os.environ.get('SLURM_ARRAY_TASK_ID')
    if slurm_id:
        with open(TARGET_LIST) as f:
            targets = f.read().strip().split('\n')
        target = targets[int(slurm_id) - 1]
    else:
        target = sys.argv[2] if len(sys.argv) > 2 else '1A2K'

    outfile = os.path.join(OUTDIR, f"rosetta_molprobity_{pipeline}_{target}.tsv")

    if not os.path.exists(outfile):
        with open(outfile, 'w') as f:
            f.write("target\tpipeline\tsrc_type\tprotocol\trep\tclashscore\t"
                    "rama_outliers\trama_favored\trota_outliers\t"
                    "molprobity_score\tcbeta_outliers\trms_bonds\trms_angles\n")

    if pipeline == 'blue':
        rosetta_dir = os.path.join(BLUE_BASE, target, 'rosetta_out')
    else:
        rosetta_dir = os.path.join(GREEN_BASE, target, 'rosetta_out')

    if not os.path.isdir(rosetta_dir):
        print(f"WARN: No rosetta_out for {target} ({pipeline})")
        return

    print(f"=== Red Rosetta MolProbity: {target} ({pipeline}) ===")

    count = 0
    for src_dir in sorted(glob.glob(os.path.join(rosetta_dir, '*/'))):
        src_name = os.path.basename(src_dir.rstrip('/'))

        for protocol in PROTOCOLS:
            protocol_dir = os.path.join(src_dir, protocol)
            if not os.path.isdir(protocol_dir):
                continue

            for pdbgz in sorted(glob.glob(os.path.join(protocol_dir, '*_r[1-5].pdb.gz'))):
                # Extract rep from filename
                import re
                m = re.search(r'_r(\d+)\.pdb\.gz$', pdbgz)
                if not m:
                    continue
                rep = m.group(1)

                run_molprobity_gz(pdbgz, target, pipeline, src_name, protocol, rep, outfile)
                count += 1

    print(f"=== Done: {count} structures scored -> {outfile} ===")


if __name__ == '__main__':
    main()
