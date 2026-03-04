#!/usr/bin/env python3
"""
Strip crystal PDBs to only keep chains present in FASTA.
Remove homo-multimer duplicate chains from crystal structures.
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

targets = sorted([f.replace('.pdb', '') for f in os.listdir(CRYSTAL_ROOT) if f.endswith('.pdb')])

stripped = 0
unchanged = 0

for target in targets:
    crystal_path = os.path.join(CRYSTAL_ROOT, f"{target}.pdb")
    fasta_path = os.path.join(AF_ROOT, target, "sequence.fasta")
    
    # Get FASTA sequences
    fasta_entries = parse_fasta(fasta_path)
    fasta_seqs = set(s for _, s in fasta_entries)
    
    # Get crystal chain sequences
    crystal_chains = parse_pdb_chains(crystal_path)
    crystal_seqs = {}
    for ch in sorted(crystal_chains.keys()):
        crystal_seqs[ch] = ''.join(three_to_one(r) for r in crystal_chains[ch])
    
    # Find which crystal chains to KEEP (exact match to a FASTA sequence)
    keep_chains = set()
    used_fasta = set()
    for ch, cseq in sorted(crystal_seqs.items()):
        if cseq in fasta_seqs and cseq not in used_fasta:
            keep_chains.add(ch)
            used_fasta.add(cseq)
    
    # If all chains kept, nothing to do
    if len(keep_chains) == len(crystal_seqs):
        unchanged += 1
        continue
    
    # Strip PDB to only keep matching chains
    backup = crystal_path + ".pre_strip"
    if not os.path.exists(backup):
        shutil.copy2(crystal_path, backup)
    
    with open(crystal_path) as f:
        lines = f.readlines()
    
    with open(crystal_path, 'w') as f:
        for line in lines:
            if line.startswith(("ATOM", "HETATM", "TER", "ANISOU")):
                chain = line[21] if len(line) > 21 else ''
                if chain in keep_chains:
                    f.write(line)
            else:
                f.write(line)
    
    removed = sorted(set(crystal_seqs.keys()) - keep_chains)
    stripped += 1
    print(f"{target}: kept {sorted(keep_chains)}, removed {removed}")

print(f"\nStripped: {stripped}")
print(f"Unchanged: {unchanged}")
print(f"Total: {stripped + unchanged}")
