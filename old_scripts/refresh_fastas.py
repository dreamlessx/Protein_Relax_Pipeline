#!/usr/bin/env python3
import os, time, shutil, sys
from pathlib import Path
import requests

UA = "Mozilla/5.0 (FASTA refresher)"
HEADERS = {"User-Agent": UA, "Accept": "text/plain,*/*;q=0.8"}

RCSB_PRIMARY   = "https://www.rcsb.org/fasta/entry/{pid}"  # current
RCSB_FALLBACK  = "https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList={pid}&compressionType=uncompressed"
PDBE_FALLBACK  = "https://www.ebi.ac.uk/pdbe/entry/pdb/{pid}/fasta?download=1"
RCSB_ENTRY_API = "https://data.rcsb.org/rest/v1/core/entry/{pid}"

def looks_like_fasta(txt: str) -> bool:
    return txt and txt.lstrip().startswith(">")

def http_get(url, timeout=30):
    return requests.get(url, headers=HEADERS, timeout=timeout, allow_redirects=True)

def resolve_replacement(pid: str) -> str | None:
    try:
        r = http_get(RCSB_ENTRY_API.format(pid=pid))
        if not r.ok:
            return None
        data = r.json()
        # Common place for replacement info
        rep = None
        for key in ("rcsb_accession_info", "rcsb_entry_container_identifiers", "rcsb_entry_info"):
            if key in data and isinstance(data[key], dict):
                rep = data[key].get("replaced_by") or data[key].get("replaced_entry_id")
                if rep:
                    break
        # Some APIs return a list
        if isinstance(rep, list) and rep:
            rep = rep[0]
        return str(rep).upper() if rep else None
    except Exception:
        return None

def fetch_fasta_for(pid: str, retries=2, backoff=1.6) -> str:
    pid = pid.upper()

    endpoints = [
        RCSB_PRIMARY.format(pid=pid),
        RCSB_FALLBACK.format(pid=pid),
        PDBE_FALLBACK.format(pid=pid),
    ]

    # Try direct endpoints
    for url in endpoints:
        for i in range(retries + 1):
            try:
                r = http_get(url)
                if r.ok and looks_like_fasta(r.text):
                    return r.text
            except requests.RequestException:
                pass
            time.sleep(backoff ** i)

    # Try resolving obsoleted -> replaced_by and re-fetch
    new_id = resolve_replacement(pid)
    if new_id and new_id != pid:
        for url in [
            RCSB_PRIMARY.format(pid=new_id),
            RCSB_FALLBACK.format(pid=new_id),
            PDBE_FALLBACK.format(pid=new_id),
        ]:
            try:
                r = http_get(url)
                if r.ok and looks_like_fasta(r.text):
                    return r.text
            except requests.RequestException:
                pass

    raise RuntimeError(f"no FASTA found for {pid}")

def refresh_folder(folder: Path):
    # collect IDs from existing .fasta names
    ids = [fp.stem for fp in folder.glob("*.fasta")]
    if not ids:
        print(f"[SKIP] {folder} (no .fasta files found)")
        return

    # wipe folder
    for item in folder.iterdir():
        if item.is_dir():
            shutil.rmtree(item)
        else:
            item.unlink(missing_ok=True)

    # re-download
    ok, fail = 0, 0
    for pid in ids:
        try:
            fasta = fetch_fasta_for(pid)
            (folder / f"{pid}.fasta").write_text(fasta)
            ok += 1
            print(f"[OK]  {folder.name}/{pid}.fasta")
        except Exception as e:
            fail += 1
            print(f"[FAIL] {folder.name}/{pid}: {e}")
    print(f"== {folder.name}: {ok} ok, {fail} failed ==")

def main(root):
    root = Path(root).expanduser().resolve()
    subs = [p for p in root.iterdir() if p.is_dir()]
    if not subs:
        print(f"No subfolders in {root}")
        return
    for sub in sorted(subs):
        refresh_folder(sub)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python3 refresh_fastas.py <root_folder_with_subdirs>")
    main(sys.argv[1])
