#!/usr/bin/env python3
"""
backfill_green_molprobity.py - Red Analysis Pipeline

Compute missing Green pre-Rosetta MolProbity for crystal + af_unrelaxed targets.
Uses same mmtbx Python API approach as compute_molprobity_py.py.

Appends results to combined_molprobity.tsv.
"""

import os
import sys
import math
import tempfile
import gzip
import numpy as np
import pandas as pd

# === cctbx environment setup (must come before cctbx imports) ===
conda_prefix = os.environ.get('CONDA_PREFIX', '')
if not conda_prefix:
    conda_prefix = os.path.dirname(os.path.dirname(sys.executable))
mon_lib_path = os.path.join(conda_prefix, 'share', 'ccp4_mon_lib')
if os.path.exists(mon_lib_path):
    os.environ['MMTBX_CCP4_MONOMER_LIB'] = mon_lib_path

# === Monkey-patch 1: SS disulfide link alias ===
import mmtbx.monomer_library.server as _mls
_orig_server_init = _mls.server.__init__
def _patched_server_init(self, *args, **kwargs):
    _orig_server_init(self, *args, **kwargs)
    if "SS" not in self.link_link_id_dict and "disulf" in self.link_link_id_dict:
        self.link_link_id_dict["SS"] = self.link_link_id_dict["disulf"]
_mls.server.__init__ = _patched_server_init

# === Monkey-patch 2: probe/reduce module detection ===
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

GREEN_BASE = "/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/data"
CRYSTAL_DIR = "/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/cleaned"
METRICS_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"


def compute_mp_score(clash_val, rota_pct, rama_fav_pct):
    """MolProbity composite score (Chen et al. 2010)."""
    try:
        return (0.426 * math.log(1 + float(clash_val))
                + 0.33 * math.log(1 + max(0, float(rota_pct) - 1))
                + 0.25 * math.log(1 + max(0, 100 - float(rama_fav_pct) - 2)))
    except (ValueError, TypeError):
        return float('nan')


def run_molprobity(pdb_path):
    """Run MolProbity validators on a PDB file, return metrics dict."""
    tmppath = None
    try:
        # Strip hydrogen atoms
        pdb_inp = iotbx.pdb.input(file_name=pdb_path)
        hierarchy = pdb_inp.construct_hierarchy()
        h_sel = hierarchy.atom_selection_cache().selection("not (element H or element D)")
        hierarchy_stripped = hierarchy.select(h_sel)

        # Write stripped PDB to temp file and re-read
        with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as tmp:
            tmp.write(hierarchy_stripped.as_pdb_string())
            tmppath = tmp.name

        pdb_inp = iotbx.pdb.input(file_name=tmppath)
        model = mmtbx.model.manager(model_input=pdb_inp)
        model.process(make_restraints=True)
        hierarchy = model.get_hierarchy()

        # Ramachandran
        rama = ramalyze.ramalyze(pdb_hierarchy=hierarchy, outliers_only=False)
        rama_out = rama.percent_outliers
        rama_fav = rama.percent_favored

        # Rotamer
        rota = rotalyze.rotalyze(pdb_hierarchy=hierarchy, outliers_only=False)
        rota_out = rota.percent_outliers

        # C-beta
        cbeta = cbetadev.cbetadev(pdb_hierarchy=hierarchy, outliers_only=True)
        cbeta_count = cbeta.get_outlier_count()

        # Clashscore
        try:
            clash = clashscore.clashscore(pdb_hierarchy=hierarchy)
            clash_val = clash.get_clashscore()
        except Exception:
            clash_val = float('nan')

        # Geometry (bonds + angles)
        geom = model.get_restraints_manager()
        if geom is not None:
            energies = geom.geometry.energies_sites(
                sites_cart=model.get_sites_cart(), compute_gradients=False)
            rms_bonds = energies.bond_deviations()[2]
            rms_angles = energies.angle_deviations()[2]
        else:
            rms_bonds = float('nan')
            rms_angles = float('nan')

        mp_score = compute_mp_score(clash_val, rota_out, rama_fav)

        return {
            'clashscore': round(clash_val, 2) if not np.isnan(clash_val) else float('nan'),
            'rama_outliers': round(rama_out, 2),
            'rama_favored': round(rama_fav, 2),
            'rota_outliers': round(rota_out, 2),
            'molprobity_score': round(mp_score, 2) if not np.isnan(mp_score) else float('nan'),
            'cbeta_outliers': cbeta_count,
            'rms_bonds': round(rms_bonds, 4) if not np.isnan(rms_bonds) else float('nan'),
            'rms_angles': round(rms_angles, 4) if not np.isnan(rms_angles) else float('nan'),
        }
    except Exception as e:
        return {'error': str(e)}
    finally:
        if tmppath and os.path.exists(tmppath):
            os.unlink(tmppath)


def main():
    mp = pd.read_csv(os.path.join(METRICS_DIR, "combined_molprobity.tsv"), sep='\t')
    green_mp = mp[mp['pipeline'] == 'green']
    all_targets = sorted(mp['target'].unique())

    green_crystal = set(green_mp[green_mp['source'] == 'crystal']['target'].unique())
    green_af = set(green_mp[green_mp['source'] == 'af_unrelaxed']['target'].unique())
    missing_crystal = sorted(set(all_targets) - green_crystal)
    missing_af = sorted(set(all_targets) - green_af)

    print(f"Missing Green crystal: {len(missing_crystal)}")
    print(f"Missing Green af_unrelaxed: {len(missing_af)}")

    new_rows = []

    # Crystal backfill
    print("\n=== Crystal backfill ===")
    for target in missing_crystal:
        pdb = os.path.join(CRYSTAL_DIR, f"{target}.pdb")
        if not os.path.exists(pdb):
            print(f"  {target}: NO PDB")
            continue

        metrics = run_molprobity(pdb)
        if 'error' in metrics:
            print(f"  {target}: FAILED - {metrics['error'][:80]}")
            continue

        row = {'target': target, 'pipeline': 'green', 'source': 'crystal', 'model_idx': 0}
        row.update(metrics)
        new_rows.append(row)
        print(f"  {target}: OK  clash={metrics.get('clashscore', 'NA')}")

    # AF unrelaxed backfill
    print(f"\n=== AF unrelaxed backfill ({len(missing_af)} targets) ===")
    done = 0
    for target in missing_af:
        af_dir = os.path.join(GREEN_BASE, target, 'af_out_unrelaxed')
        if not os.path.isdir(af_dir):
            continue

        for model_idx in range(5):
            pdb = os.path.join(af_dir, f"ranked_{model_idx}.pdb")
            if not os.path.exists(pdb):
                continue

            metrics = run_molprobity(pdb)
            if 'error' in metrics:
                continue

            row = {'target': target, 'pipeline': 'green', 'source': 'af_unrelaxed',
                   'model_idx': model_idx}
            row.update(metrics)
            new_rows.append(row)

        done += 1
        if done % 10 == 0:
            print(f"  {done}/{len(missing_af)} targets, {len(new_rows)} rows...")

    print(f"\n=== Done: {len(new_rows)} new rows ===")

    if new_rows:
        new_df = pd.DataFrame(new_rows)
        cols = list(mp.columns)
        for c in cols:
            if c not in new_df.columns:
                new_df[c] = np.nan
        new_df = new_df[cols]

        combined = pd.concat([mp, new_df], ignore_index=True)
        outpath = os.path.join(METRICS_DIR, "combined_molprobity.tsv")
        combined.to_csv(outpath, sep='\t', index=False)
        print(f"Updated {outpath}: {len(combined)} total rows")

        green = combined[combined['pipeline'] == 'green']
        for src in ['crystal', 'af_unrelaxed', 'af_relaxed', 'boltz', 'amber_af', 'amber_boltz']:
            n = green[green['source'] == src]['target'].nunique()
            print(f"  Green {src}: {n} targets")
    else:
        print("No new rows to add.")


if __name__ == '__main__':
    main()
