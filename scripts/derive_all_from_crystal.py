#!/usr/bin/env python3
"""
Derive ALL 257 FASTAs from crystal PDB coordinates.
This ensures crystal == FASTA for every target.
"""

import os
import shutil
from collections import defaultdict

AF_ROOT = "/data/p_csb_meiler/agarwm5/af_work"
CRYSTAL_ROOT = "/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/cleaned"

def three_to_one(resname):
    mapping = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
        'MSE': 'M', 'SEC': 'U', 'PYL': 'O'
    }
    return mapping.get(resname, 'X')

def parse_pdb_chains(pdb_path):
    chains = defaultdict(list)
    prev_resnum = {}
    for line in open(pdb_path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            chain = line[21]
            resnum = int(line[22:26])
            resname = line[17:20].strip()
            if chain not in prev_resnum or resnum != prev_resnum[chain]:
                chains[chain].append(resname)
                prev_resnum[chain] = resnum
    return chains

def parse_fasta(fasta_path):
    entries = []
    header = None
    seq = []
    for line in open(fasta_path):
        line = line.strip()
        if line.startswith('>'):
            if header is not None:
                entries.append((header, ''.join(seq)))
            header = line
            seq = []
        elif line:
            seq.append(line)
    if header is not None:
        entries.append((header, ''.join(seq)))
    return entries

def write_fasta(entries, path):
    with open(path, 'w') as f:
        for header, seq in entries:
            f.write(f"{header}\n")
            for i in range(0, len(seq), 80):
                f.write(f"{seq[i:i+80]}\n")

def write_boltz_fasta(entries, path):
    chain_letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    with open(path, 'w') as f:
        for i, (_, seq) in enumerate(entries):
            if i < len(chain_letters):
                f.write(f">{chain_letters[i]}|PROTEIN|\n")
            else:
                f.write(f">chain{i}|PROTEIN|\n")
            for j in range(0, len(seq), 80):
                f.write(f"{seq[j:j+80]}\n")

def get_unique_chains(crystal_seqs):
    """Get unique protein chain sequences (deduplicate homo-multimers)."""
    unique = {}
    for ch in sorted(crystal_seqs.keys()):
        seq = crystal_seqs[ch]
        is_dup = False
        for prev_ch, prev_seq in unique.items():
            # Check if sequences are essentially the same
            if seq == prev_seq:
                is_dup = True
                break
            if seq in prev_seq or prev_seq in seq:
                is_dup = True
                break
            # Check high identity (>95% in overlapping region)
            shorter = min(len(seq), len(prev_seq))
            if shorter >= 20:
                matches = sum(1 for a, b in zip(seq[:shorter], prev_seq[:shorter]) if a == b)
                if matches / shorter > 0.95:
                    is_dup = True
                    break
        if not is_dup:
            unique[ch] = seq
    return unique

targets = sorted([f.replace('.pdb', '') for f in os.listdir(CRYSTAL_ROOT) if f.endswith('.pdb')])
print(f"Processing {len(targets)} targets\n")

changed = []
unchanged = []
errors = []

for target in targets:
    crystal_path = os.path.join(CRYSTAL_ROOT, f"{target}.pdb")
    fasta_path = os.path.join(AF_ROOT, target, "sequence.fasta")
    boltz_path = os.path.join(AF_ROOT, target, "boltz_input.fasta")
    
    if not os.path.exists(fasta_path):
        errors.append(f"{target}: no sequence.fasta")
        continue
    
    # Parse crystal
    crystal_chains = parse_pdb_chains(crystal_path)
    crystal_seqs = {}
    for ch in sorted(crystal_chains.keys()):
        crystal_seqs[ch] = ''.join(three_to_one(r) for r in crystal_chains[ch])
    
    # Get unique chains
    unique = get_unique_chains(crystal_seqs)
    
    # Read current FASTA
    old_entries = parse_fasta(fasta_path)
    old_seqs = sorted([s for _, s in old_entries])
    
    # Create new crystal-derived entries
    new_entries = []
    for ch, seq in sorted(unique.items()):
        header = f">{target}|Chain {ch}|Crystal-derived|"
        new_entries.append((header, seq))
    
    new_seqs = sorted([s for _, s in new_entries])
    
    # Check if anything changed
    if old_seqs == new_seqs:
        unchanged.append(target)
    else:
        # Backup
        backup = fasta_path + ".pre_crystal_all"
        if not os.path.exists(backup):
            shutil.copy2(fasta_path, backup)
        
        write_fasta(new_entries, fasta_path)
        
        if os.path.exists(boltz_path):
            boltz_backup = boltz_path + ".pre_crystal_all"
            if not os.path.exists(boltz_backup):
                shutil.copy2(boltz_path, boltz_backup)
        write_boltz_fasta(new_entries, boltz_path)
        
        changed.append(target)

print(f"Changed: {len(changed)}")
print(f"Unchanged: {len(unchanged)}")
print(f"Errors: {len(errors)}")

if errors:
    print(f"\nErrors:")
    for e in errors:
        print(f"  {e}")

# Write new rerun dirlist for ALL changed targets
rerun_path = "/data/p_csb_meiler/agarwm5/protein_pipeline/af_crystal_derived_dirlist.txt"
with open(rerun_path, 'w') as f:
    for t in changed:
        f.write(f"{AF_ROOT}/{t}\n")
print(f"\nAF rerun dirlist: {rerun_path} ({len(changed)} targets)")

boltz_rerun_path = "/data/p_csb_meiler/agarwm5/protein_pipeline/boltz_crystal_derived_dirlist.txt"
with open(boltz_rerun_path, 'w') as f:
    for t in changed:
        f.write(f"{AF_ROOT}/{t}\n")
print(f"Boltz rerun dirlist: {boltz_rerun_path} ({len(changed)} targets)")

# Show what changed
if changed:
    print(f"\nChanged targets:")
    for t in changed:
        print(f"  {t}")
