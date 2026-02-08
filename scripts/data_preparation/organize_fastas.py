#!/usr/bin/env python3
"""
Organize FASTA files into per-PDB subdirectories.

Usage:
    python3 organize_fastas.py <fasta_dir> <data_dir> [--move] [--rename sequence.fasta] [--overwrite]
"""
import argparse
from pathlib import Path
import shutil

def main():
    ap = argparse.ArgumentParser(description="Put each FASTA into its own folder under a global data dir.")
    ap.add_argument("fasta_dir", help="Directory containing *.fasta files")
    ap.add_argument("data_dir", help="Output 'data' directory (created if needed)")
    ap.add_argument("--move", action="store_true", help="Move files instead of copying (default: copy)")
    ap.add_argument("--rename", default=None, help="Rename inside each folder (e.g., sequence.fasta)")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing files")
    args = ap.parse_args()

    src = Path(args.fasta_dir).expanduser().resolve()
    out_root = Path(args.data_dir).expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    fastas = sorted(src.glob("*.fasta"))
    if not fastas:
        print(f"No .fasta files found in {src}")
        return

    ok = skip = 0
    for fp in fastas:
        pdb_id = fp.stem.upper()
        dest_dir = out_root / pdb_id
        dest_dir.mkdir(parents=True, exist_ok=True)

        dest_name = args.rename if args.rename else f"{pdb_id}.fasta"
        dest = dest_dir / dest_name

        if dest.exists() and not args.overwrite:
            print(f"[SKIP] {pdb_id} (exists)")
            skip += 1
            continue

        if args.move:
            shutil.move(str(fp), str(dest))
        else:
            shutil.copy2(str(fp), str(dest))

        print(f"[OK] {pdb_id} -> {dest}")
        ok += 1

    print(f"\nDone. Created/updated {ok} entries; skipped {skip}.")

if __name__ == "__main__":
    main()
