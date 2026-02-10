# /dors/meilerlab/home/agarwm5/benchmarking/scripts/clean_all_pdbs.sh
#!/usr/bin/env bash
# Clean every PDB in merged/ using Rosetta clean_pdb.py,
# but write ONE combined cleaned file per complex into cleaned/.
set -euo pipefail

ROOT="/dors/meilerlab/home/agarwm5/benchmarking"
MERGED="$ROOT/merged"
CLEANED="$ROOT/cleaned"
CLEAN_SCRIPT="/dors/meilerlab/apps/rosetta/rosetta-3.14/main/tools/protein_tools/scripts/clean_pdb.py"

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

  # Temporary dir for per-chain Rosetta outputs
  tmpdir="$(mktemp -d "$CLEANED/.tmp_${stem}_XXXX")"

  # 1) Run Rosetta per chain (this script requires exactly one chain per call)
  for ch in "${CHAINS[@]}"; do
    ( cd "$tmpdir" && python "$CLEAN_SCRIPT" "$pdb" "$ch" > /dev/null )
  done

  # 2) Merge cleaned chains -> single file; renumber atom serials; add TER between chains
  out="$CLEANED/${stem}.pdb"
  python - "$tmpdir" "$out" "${CHAINS[@]}" <<'PY'
import sys, os, glob
tmp = sys.argv[1]
out = sys.argv[2]
chains = sys.argv[3:]

serial = 1
with open(out, 'w') as fo:
    for ch in chains:
        # Rosetta names like <stem>_<CHAIN>.pdb (CHAIN may be '_')
        files = sorted(glob.glob(os.path.join(tmp, f"*_{ch}.pdb")))
        if not files:
            continue
        with open(files[0]) as f:
            for line in f:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                # renumber atom serials (cols 7–11)
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

