#!/usr/bin/env bash
# Clean PDB files using Rosetta's clean_pdb.py.
# Writes ONE combined cleaned file per complex.
#
# Usage:
#   bash clean_pdbs.sh <merged_dir> <cleaned_dir> <rosetta_clean_script>
#
# Example:
#   bash clean_pdbs.sh ./merged ./cleaned /path/to/rosetta/tools/protein_tools/scripts/clean_pdb.py

set -euo pipefail

MERGED="${1:?Usage: clean_pdbs.sh <merged_dir> <cleaned_dir> <rosetta_clean_script>}"
CLEANED="${2:?Usage: clean_pdbs.sh <merged_dir> <cleaned_dir> <rosetta_clean_script>}"
CLEAN_SCRIPT="${3:?Usage: clean_pdbs.sh <merged_dir> <cleaned_dir> <rosetta_clean_script>}"

mkdir -p "$CLEANED"
shopt -s nullglob

for pdb in "$MERGED"/*.pdb; do
  base=$(basename "$pdb")
  stem="${base%.pdb}"

  # Get unique chain IDs from ATOM records (col 22). Map blank -> '_' for Rosetta.
  mapfile -t CHAINS < <(
    awk '/^ATOM/{c=substr($0,22,1); if(c==" ") c="_"; print toupper(c)}' "$pdb" | sort -u
  )

  if (( ${#CHAINS[@]} == 0 )); then
    echo "[SKIP] $base (no ATOM records)"
    continue
  fi

  echo "[CLEAN] $base  chains: ${CHAINS[*]}"

  tmpdir="$(mktemp -d "$CLEANED/.tmp_${stem}_XXXX")"

  for ch in "${CHAINS[@]}"; do
    ( cd "$tmpdir" && python "$CLEAN_SCRIPT" "$pdb" "$ch" > /dev/null )
  done

  out="$CLEANED/${stem}.pdb"
  python - "$tmpdir" "$out" "${CHAINS[@]}" <<'PY'
import sys, os, glob
tmp = sys.argv[1]
out = sys.argv[2]
chains = sys.argv[3:]

serial = 1
with open(out, 'w') as fo:
    for ch in chains:
        files = sorted(glob.glob(os.path.join(tmp, f"*_{ch}.pdb")))
        if not files:
            continue
        with open(files[0]) as f:
            for line in f:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                line = f"{line[:6]}{serial:5d}{line[11:]}"
                fo.write(line.rstrip("\n") + "\n")
                serial += 1
        fo.write("TER\n")
    fo.write("END\n")
PY

  rm -rf "$tmpdir"
  echo "  -> $out"
done

echo "Done. Combined cleaned PDBs in: $CLEANED"
