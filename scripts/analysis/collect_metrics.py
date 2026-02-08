"""
PyMOL plugin: Collect RMSD and Rosetta energy metrics for loaded structures.

Calculates RMSD of each loaded object to a reference structure and extracts
Rosetta total_score from PDB REMARK lines or an external scorefile.

Usage (inside PyMOL):
    run collect_metrics.py
    collect_metrics ref=crystal, sel=polymer and name CA, do_fit=1, scorefile=relax.fasc, out=metrics.tsv
"""
from pymol import cmd
import os, re

def _parse_rosetta_total_from_pdbstr(pdbstr: str):
    if not pdbstr:
        return None
    m = re.search(r'(?:^|\n)REMARK.*?(?:SCORE|score).*?total_score\s+([+-]?\d+(?:\.\d+)?)', pdbstr)
    if m:
        try: return float(m.group(1))
        except: pass
    lines = [l[6:].strip() for l in pdbstr.splitlines() if l.startswith('REMARK')]
    hdr = None
    for i, txt in enumerate(lines):
        toks = txt.split()
        if 'total_score' in toks and ('score' in toks or 'fa_atr' in toks):
            hdr = toks
            for j in range(i+1, len(lines)):
                dtoks = lines[j].split()
                if not dtoks:
                    continue
                try:
                    idx = hdr.index('total_score')
                    return float(dtoks[idx])
                except Exception:
                    continue
            break
    return None

def _parse_score_sc(path):
    if not path or not os.path.isfile(path):
        return {}
    mapping = {}
    header = None
    with open(path) as fh:
        for ln in fh:
            if not ln.startswith('SCORE:'):
                continue
            toks = ln.split()
            if header is None or ('description' in toks and not toks[1].replace('.','',1).lstrip('+-').isdigit()):
                header = toks[1:]
                continue
            if header:
                data = toks[1:]
                if len(data) >= len(header):
                    cols = dict(zip(header, data))
                    desc = cols.get('description') or cols.get('decoy') or cols.get('tag') or cols.get('desc')
                    score = cols.get('total_score') or cols.get('score')
                    if desc and score:
                        try: mapping[desc] = float(score)
                        except: pass
    return mapping

def collect_metrics(ref='crystal', sel='polymer and name CA', do_fit=1, scorefile=None, out='metrics_to_crystal.tsv'):
    do_fit = bool(int(do_fit)) if isinstance(do_fit, str) else bool(do_fit)
    scoremap = _parse_score_sc(scorefile) if scorefile else {}

    objs = [o for o in cmd.get_names('objects') if cmd.count_atoms(o) > 0]
    if ref not in objs:
        print(f'ERROR: reference object "{ref}" not found among objects: {objs}')
        return

    rows = []
    for obj in objs:
        if obj == ref:
            rmsd, pairs = float('nan'), 0
        else:
            if do_fit:
                rmsd, pairs, *_ = cmd.align(f'{obj} and {sel}', f'{ref} and {sel}',
                                            cycles=0, transform=1, object=None, quiet=1)
            else:
                rmsd = cmd.rms_cur(f'{obj} and {sel}', f'{ref} and {sel}', quiet=1)
                pairs = 0

        energy = scoremap.get(obj, None)
        if energy is None:
            pdbstr = cmd.get_pdbstr(obj)
            energy = _parse_rosetta_total_from_pdbstr(pdbstr)

        rows.append(dict(object=obj, rmsd=rmsd, pairs=pairs, energy=energy))

    print(f'Object\tRMSD_to_{ref}\tPairs\tEnergy(total_score)')
    def _rkey(r):
        return (1, 0) if r["object"]==ref else (0, r["rmsd"] if r["rmsd"]==r["rmsd"] else 9e9)
    for r in sorted(rows, key=_rkey):
        r_rmsd = 'NA' if r['rmsd']!=r['rmsd'] else f'{r["rmsd"]:.3f}'
        r_energy = 'NA' if r['energy'] is None else f'{r["energy"]:.3f}'
        print(f'{r["object"]}\t{r_rmsd}\t{r["pairs"]}\t{r_energy}')

    with open(out, 'w') as fh:
        fh.write(f'object\trmsd_to_{ref}\tpairs\tenergy_total_score\n')
        for r in rows:
            r_rmsd = '' if r['rmsd']!=r['rmsd'] else f'{r["rmsd"]:.6f}'
            r_energy = '' if r['energy'] is None else f'{r["energy"]:.6f}'
            fh.write(f'{r["object"]}\t{r_rmsd}\t{r["pairs"]}\t{r_energy}\n')
    print('Saved:', os.path.abspath(out))

cmd.extend('collect_metrics', collect_metrics)
