#!/usr/bin/env python3
"""
Convert RCSB-format FASTA files to Boltz-1 input format.

Boltz-1 requires headers in the format: >{CHAIN}|PROTEIN|
This script parses chain IDs from standard RCSB FASTA headers and
rewrites each sequence with the Boltz-compatible header.

Usage:
    python3 prepare_boltz_fastas.py <data_dir>

Where <data_dir> contains subdirectories, each with a sequence.fasta file.
Outputs boltz_input.fasta in each subdirectory.
"""
import re
import sys
from pathlib import Path

def parse_chains(header: str):
    m = re.search(r'Chain[s]?\s+([A-Za-z0-9, ]+)', header, re.IGNORECASE)
    if not m:
        parts = header.split('|')
        if len(parts) > 1:
            t = parts[1].strip()
            if t and len(t) >= 1:
                return [t[0].upper()]
        return ['A']
    raw = m.group(1)
    chains = [c.strip().upper() for c in re.split(r'[,\s]+', raw) if c.strip()]
    return chains or ['A']

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 prepare_boltz_fastas.py <data_dir>")
        sys.exit(1)

    root = Path(sys.argv[1]).expanduser().resolve()

    for fasta in sorted(root.glob("*/sequence.fasta")):
        lines = fasta.read_text().splitlines()
        out = []
        seq = []
        chains = []
        for line in lines:
            if line.startswith('>'):
                if seq and chains:
                    seqstr = ''.join(seq)
                    for ch in chains:
                        out.append(f">{ch}|PROTEIN|")
                        out.append(seqstr)
                    seq = []
                chains = parse_chains(line)
            else:
                seq.append(line.strip())

        if seq and chains:
            seqstr = ''.join(seq)
            for ch in chains:
                out.append(f">{ch}|PROTEIN|")
                out.append(seqstr)

        dest = fasta.parent / "boltz_input.fasta"
        dest.write_text('\n'.join(out) + '\n')
        print(f"[OK] {dest}")

if __name__ == "__main__":
    main()
