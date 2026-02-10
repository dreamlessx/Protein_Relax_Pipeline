#!/usr/bin/env python3
import re
from pathlib import Path

root = Path("/home/agarwm5/benchmarking/data").expanduser().resolve()

def parse_chains(header: str):
    m = re.search(r'Chain[s]?\s+([A-Za-z0-9, ]+)', header, re.IGNORECASE)
    if not m:
        # fallback: try to extract something sane after the first '|'
        parts = header.split('|')
        if len(parts) > 1:
            t = parts[1].strip()
            if t and len(t) >= 1:
                return [t[0].upper()]
        return ['A']
    raw = m.group(1)
    chains = [c.strip().upper() for c in re.split(r'[,\s]+', raw) if c.strip()]
    return chains or ['A']

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
