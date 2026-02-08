#!/usr/bin/env python3
"""
Download FASTA sequences from RCSB/PDBe for a directory of PDB files.

Tries multiple endpoints (RCSB primary, RCSB fallback, PDBe) and handles
obsolete entries by looking up replacement IDs.

Usage:
    python3 download_fastas.py <merged_dir> [fasta_out_dir] [--overwrite]
"""
import sys, re, time, csv
from pathlib import Path
import requests

UA = "Mozilla/5.0 (pdbs_to_fastas.py)"
HEADERS = {"User-Agent": UA, "Accept": "text/plain,*/*;q=0.8"}

RCSB_PRIMARY   = "https://www.rcsb.org/fasta/entry/{pid}"
RCSB_FALLBACK  = "https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList={pid}&compressionType=uncompressed"
PDBE_FALLBACK  = "https://www.ebi.ac.uk/pdbe/entry/pdb/{pid}/fasta?download=1"
RCSB_ENTRY_API = "https://data.rcsb.org/rest/v1/core/entry/{pid}"

PDB_ID_RE = re.compile(r"([0-9][A-Za-z0-9]{3})")

def looks_like_fasta(txt: str) -> bool:
    return bool(txt) and txt.lstrip().startswith(">")

def http_get(url, timeout=30):
    return requests.get(url, headers=HEADERS, timeout=timeout, allow_redirects=True)

def resolve_replacement(pid: str) -> str | None:
    """If an entry is obsolete, try to find its replacement ID."""
    try:
        r = http_get(RCSB_ENTRY_API.format(pid=pid))
        if not r.ok:
            return None
        data = r.json()
        rep = None
        for key in ("rcsb_accession_info", "rcsb_entry_container_identifiers", "rcsb_entry_info"):
            if key in data and isinstance(data[key], dict):
                rep = data[key].get("replaced_by") or data[key].get("replaced_entry_id")
                if rep:
                    break
        if isinstance(rep, list) and rep:
            rep = rep[0]
        return str(rep).upper() if rep else None
    except Exception:
        return None

def fetch_fasta(pid: str, retries=2, backoff=1.5) -> str:
    pid = pid.upper()
    endpoints = [
        RCSB_PRIMARY.format(pid=pid),
        RCSB_FALLBACK.format(pid=pid),
        PDBE_FALLBACK.format(pid=pid),
    ]
    for url in endpoints:
        for i in range(retries + 1):
            try:
                r = http_get(url)
                if r.ok and looks_like_fasta(r.text):
                    return r.text
            except requests.RequestException:
                pass
            time.sleep(backoff ** i)
    rep = resolve_replacement(pid)
    if rep and rep != pid:
        for url in [
            RCSB_PRIMARY.format(pid=rep),
            RCSB_FALLBACK.format(pid=rep),
            PDBE_FALLBACK.format(pid=rep),
        ]:
            try:
                r = http_get(url)
                if r.ok and looks_like_fasta(r.text):
                    return r.text
            except requests.RequestException:
                pass
    raise RuntimeError(f"no FASTA found for {pid}")

def extract_pdb_id(name: str) -> str | None:
    m = PDB_ID_RE.search(name.upper())
    return m.group(1) if m else None

def main():
    if len(sys.argv) < 2 or len(sys.argv) > 4:
        print("Usage: python3 download_fastas.py <merged_dir> [fasta_out_dir] [--overwrite]")
        sys.exit(1)

    merged_dir = Path(sys.argv[1]).expanduser().resolve()
    if len(sys.argv) >= 3 and not sys.argv[2].startswith("--"):
        fasta_dir = Path(sys.argv[2]).expanduser().resolve()
        overwrite = (len(sys.argv) == 4 and sys.argv[3] == "--overwrite")
    else:
        fasta_dir = merged_dir.parent / "fasta"
        overwrite = (len(sys.argv) >= 3 and sys.argv[-1] == "--overwrite")

    fasta_dir.mkdir(parents=True, exist_ok=True)

    pdb_files = sorted([p for p in merged_dir.glob("*") if p.is_file() and p.suffix.lower() in (".pdb", ".cif", ".ent")])
    if not pdb_files:
        print(f"No PDB/mmCIF files found in {merged_dir}")
        sys.exit(0)

    log_path = fasta_dir / "fasta_download_log.csv"
    with log_path.open("w", newline="") as logf:
        w = csv.writer(logf)
        w.writerow(["pdb_id", "source_file", "status", "note"])

        ok = fail = skip = 0
        for fp in pdb_files:
            pid = extract_pdb_id(fp.stem)
            if not pid:
                w.writerow(["", str(fp.name), "skip", "no PDB id found"])
                print(f"[SKIP] {fp.name} (no PDB id found)")
                skip += 1
                continue

            out = fasta_dir / f"{pid}.fasta"
            if out.exists() and not overwrite:
                w.writerow([pid, fp.name, "skip", "exists"])
                print(f"[SKIP] {pid} (exists)")
                skip += 1
                continue

            try:
                fasta = fetch_fasta(pid)
                out.write_text(fasta)
                w.writerow([pid, fp.name, "ok", ""])
                print(f"[OK]   {pid} -> {out.name}")
                ok += 1
            except Exception as e:
                w.writerow([pid, fp.name, "fail", str(e)])
                print(f"[FAIL] {pid}: {e}")
                fail += 1
            time.sleep(0.15)

    print(f"\nDone. OK={ok}, FAIL={fail}, SKIP={skip}")
    print(f"Log: {log_path}")

if __name__ == "__main__":
    main()
