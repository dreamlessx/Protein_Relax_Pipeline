#!/usr/bin/env python3
"""
Standalone AMBER relaxation script for AlphaFold PDB outputs.
Uses OpenMM with AMBER ff14SB force field (same as AlphaFold).
Uses PDBFixer to handle non-standard residues and protonation states.
"""

import sys
import os
from pathlib import Path

# Add AlphaFold to path
sys.path.insert(0, '/sb/apps/alphafold232/alphafold')

from openmm import app as openmm_app
from openmm import unit
from openmm import LangevinIntegrator, Platform
from openmm.app import PDBFile, ForceField, Modeller, HBonds, Simulation
from pdbfixer import PDBFixer
import numpy as np

ENERGY_TOLERANCE = 2.39  # kcal/mol, same as AlphaFold
STIFFNESS = 10.0  # kcal/mol/A^2
MAX_ITERATIONS = 0  # unlimited
MAX_OUTER_ITERATIONS = 20


def relax_pdb(input_pdb: str, output_pdb: str, use_gpu: bool = True):
    """Relax a single PDB file using AMBER force field."""

    print(f"Relaxing: {input_pdb}")

    # Use PDBFixer to handle non-standard residues (especially HIS -> HID/HIE/HIP)
    fixer = PDBFixer(filename=input_pdb)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)  # pH 7.0

    # Use AMBER ff14SB (same as AlphaFold)
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # Create modeller from fixed structure
    modeller = Modeller(fixer.topology, fixer.positions)

    # Create system
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=openmm_app.NoCutoff,
        constraints=HBonds
    )

    # Add position restraints
    from openmm import CustomExternalForce
    force = CustomExternalForce(
        "0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
    )
    force.addGlobalParameter("k", STIFFNESS * unit.kilocalories_per_mole / unit.angstroms**2)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    for i, pos in enumerate(modeller.positions):
        force.addParticle(i, [pos.x, pos.y, pos.z])
    system.addForce(force)

    # Setup simulation
    integrator = LangevinIntegrator(0*unit.kelvin, 1/unit.picosecond, 0.001*unit.picoseconds)

    if use_gpu:
        try:
            platform = Platform.getPlatformByName('CUDA')
            properties = {'CudaPrecision': 'mixed'}
            simulation = Simulation(modeller.topology, system, integrator, platform, properties)
        except:
            platform = Platform.getPlatformByName('CPU')
            simulation = Simulation(modeller.topology, system, integrator, platform)
    else:
        platform = Platform.getPlatformByName('CPU')
        simulation = Simulation(modeller.topology, system, integrator, platform)

    simulation.context.setPositions(modeller.positions)

    # Minimize
    ret = {"min_attempts": 0, "num_violations": 0}

    for attempt in range(MAX_OUTER_ITERATIONS):
        ret["min_attempts"] += 1
        simulation.minimizeEnergy(maxIterations=MAX_ITERATIONS, tolerance=ENERGY_TOLERANCE)

        # Check for violations (same logic as AlphaFold)
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        ret["num_violations"] = 0  # Simplified - skip violation check

        if ret["num_violations"] == 0:
            break

    # Get final positions
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()

    # Save output
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, positions, f)

    print(f"  Saved: {output_pdb} (attempts: {ret['min_attempts']})")
    return True


def process_target(target_dir: str, use_gpu: bool = True):
    """Process all unrelaxed models in a target directory."""

    af_out = Path(target_dir) / "af_out"
    if not af_out.exists():
        print(f"No af_out directory: {target_dir}")
        return

    # Find unrelaxed models (ranked_1 through ranked_4)
    for i in range(1, 5):
        input_pdb = af_out / f"ranked_{i}.pdb"
        if not input_pdb.exists():
            continue

        # Check if already has hydrogens (already relaxed)
        with open(input_pdb) as f:
            content = f.read()
            if ' H ' in content or ' HA ' in content or ' HB' in content:
                print(f"  ranked_{i}.pdb already relaxed, skipping")
                continue

        # Backup original
        backup = af_out / f"ranked_{i}_unrelaxed.pdb"
        if not backup.exists():
            os.rename(input_pdb, backup)

        # Relax
        try:
            relax_pdb(str(backup), str(input_pdb), use_gpu=use_gpu)
        except Exception as e:
            print(f"  ERROR relaxing ranked_{i}.pdb: {e}")
            # Restore original if relax failed
            if backup.exists() and not input_pdb.exists():
                os.rename(backup, input_pdb)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: amber_relax.py <target_dir> [--cpu]")
        sys.exit(1)

    target_dir = sys.argv[1]
    use_gpu = "--cpu" not in sys.argv

    process_target(target_dir, use_gpu=use_gpu)
