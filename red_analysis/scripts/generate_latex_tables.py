#!/usr/bin/env python3
"""
generate_latex_tables.py — Red Analysis Pipeline

Generates publication-ready LaTeX tables for the paper.

Tables:
  table1: Pre-Rosetta structural metrics by source (TM, RMSD, GDT)
  table2: MolProbity validation by source
  table3: AMBER effect on TM-score + MolProbity (the dual-effect table)
  table4: Rosetta protocol comparison
  table5: Blue-Green reproducibility summary
  table6: Outlier targets
"""

import os
import numpy as np
import pandas as pd
from scipy import stats

METRICS_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"
OUTDIR = "/data/p_csb_meiler/agarwm5/red_analysis/tables"

SOURCE_LABELS = {
    'af_relaxed': 'AF2 (relaxed)',
    'af_unrelaxed': 'AF2 (unrelaxed)',
    'boltz': 'Boltz-1',
    'amber_af': 'AMBER(AF)',
    'amber_boltz': 'AMBER(Boltz)',
    'crystal': 'Crystal',
}

SOURCE_ORDER = ['af_relaxed', 'af_unrelaxed', 'boltz', 'amber_af', 'amber_boltz']

PROTOCOL_LABELS = {
    'cartesian_beta': 'Cartesian $\\beta$',
    'cartesian_ref15': 'Cartesian REF15',
    'dualspace_beta': 'Dualspace $\\beta$',
    'dualspace_ref15': 'Dualspace REF15',
    'normal_beta': 'Normal $\\beta$',
    'normal_ref15': 'Normal REF15',
}

PROTOCOL_ORDER = ['cartesian_beta', 'normal_beta', 'normal_ref15',
                  'cartesian_ref15', 'dualspace_beta', 'dualspace_ref15']


def load_all():
    """Load all metric datasets."""
    data = {}

    # Pre-Rosetta TM-score
    path = os.path.join(METRICS_DIR, "combined_tmscore.tsv")
    if os.path.exists(path):
        tm = pd.read_csv(path, sep='\t')
        for col in ['tmscore', 'rmsd', 'gdtts', 'gdtha']:
            tm[col] = pd.to_numeric(tm[col], errors='coerce')
        data['tm'] = tm

    # MolProbity
    path = os.path.join(METRICS_DIR, "combined_molprobity.tsv")
    if os.path.exists(path):
        mp = pd.read_csv(path, sep='\t')
        for col in ['clashscore', 'rama_outliers', 'rama_favored', 'rota_outliers',
                     'molprobity_score', 'cbeta_outliers', 'rms_bonds', 'rms_angles']:
            mp[col] = pd.to_numeric(mp[col], errors='coerce')
        data['mp'] = mp

    # Rosetta TM-score
    path = os.path.join(METRICS_DIR, "combined_rosetta_tmscore.tsv")
    if os.path.exists(path):
        rt = pd.read_csv(path, sep='\t')
        for col in ['tmscore', 'rmsd', 'gdtts', 'gdtha']:
            rt[col] = pd.to_numeric(rt[col], errors='coerce')
        data['rosetta'] = rt

    return data


def table1_pre_rosetta_metrics(data):
    """Table 1: Pre-Rosetta structural metrics by source."""
    if 'tm' not in data:
        return

    tm = data['tm'][data['tm']['pipeline'] == 'blue']
    per_target = tm.groupby(['source', 'target'])[['tmscore', 'rmsd', 'gdtts']].mean().reset_index()

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Structural similarity to crystal reference (Blue pipeline, %d targets). Values are mean $\pm$ std across per-target averages.}" % per_target['target'].nunique(),
        r"\label{tab:pre-rosetta}",
        r"\begin{tabular}{lccc}",
        r"\toprule",
        r"Source & TM-score & RMSD (\AA) & GDT-TS \\",
        r"\midrule",
    ]

    for src in SOURCE_ORDER:
        s = per_target[per_target['source'] == src]
        if len(s) == 0:
            continue
        label = SOURCE_LABELS[src]
        tm_mean, tm_std = s['tmscore'].mean(), s['tmscore'].std()
        rmsd_mean, rmsd_std = s['rmsd'].mean(), s['rmsd'].std()
        gdt_mean, gdt_std = s['gdtts'].mean(), s['gdtts'].std()
        lines.append(f"  {label} & {tm_mean:.3f} $\\pm$ {tm_std:.3f} & "
                     f"{rmsd_mean:.2f} $\\pm$ {rmsd_std:.2f} & "
                     f"{gdt_mean:.3f} $\\pm$ {gdt_std:.3f} \\\\")

    lines.extend([
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ])

    tex = "\n".join(lines)
    path = os.path.join(OUTDIR, "table1_pre_rosetta.tex")
    with open(path, 'w') as f:
        f.write(tex)
    print(f"  Saved table1_pre_rosetta.tex")
    return tex


def table2_molprobity(data):
    """Table 2: MolProbity validation by source."""
    if 'mp' not in data:
        return

    mp = data['mp'][data['mp']['pipeline'] == 'blue']
    per_target = mp.groupby(['source', 'target'])[
        ['clashscore', 'rama_favored', 'rota_outliers', 'molprobity_score']
    ].mean().reset_index()

    n_targets = per_target['target'].nunique()

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{MolProbity validation metrics by structure source (Blue pipeline, %d targets). Values are mean $\pm$ std across per-target averages.}" % n_targets,
        r"\label{tab:molprobity}",
        r"\begin{tabular}{lcccc}",
        r"\toprule",
        r"Source & Clashscore$\downarrow$ & Rama Fav\%$\uparrow$ & Rota Out\%$\downarrow$ & MP Score$\downarrow$ \\",
        r"\midrule",
    ]

    # Crystal first (as reference)
    all_sources = ['crystal'] + SOURCE_ORDER
    for src in all_sources:
        s = per_target[per_target['source'] == src]
        if len(s) == 0:
            continue
        label = SOURCE_LABELS.get(src, src.capitalize())
        cs = f"{s['clashscore'].mean():.1f} $\\pm$ {s['clashscore'].std():.1f}"
        rf = f"{s['rama_favored'].mean():.1f} $\\pm$ {s['rama_favored'].std():.1f}"
        ro = f"{s['rota_outliers'].mean():.2f} $\\pm$ {s['rota_outliers'].std():.2f}"
        mp_s = f"{s['molprobity_score'].mean():.2f} $\\pm$ {s['molprobity_score'].std():.2f}"
        lines.append(f"  {label} & {cs} & {rf} & {ro} & {mp_s} \\\\")
        if src == 'crystal':
            lines.append(r"  \midrule")

    lines.extend([
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ])

    tex = "\n".join(lines)
    path = os.path.join(OUTDIR, "table2_molprobity.tex")
    with open(path, 'w') as f:
        f.write(tex)
    print(f"  Saved table2_molprobity.tex")
    return tex


def cliffs_delta(x, y):
    """Cliff's delta effect size."""
    n_x, n_y = len(x), len(y)
    if n_x == 0 or n_y == 0:
        return float('nan')
    more = sum(1 for xi in x for yi in y if xi > yi)
    less = sum(1 for xi in x for yi in y if xi < yi)
    return (more - less) / (n_x * n_y)


def interpret_d(d):
    """Interpret Cliff's delta magnitude."""
    d_abs = abs(d)
    if d_abs < 0.147:
        return "negligible"
    elif d_abs < 0.33:
        return "small"
    elif d_abs < 0.474:
        return "medium"
    else:
        return "large"


def table3_amber_dual_effect(data):
    """Table 3: AMBER's dual effect — TM-score (none) + MolProbity (large)."""
    if 'tm' not in data or 'mp' not in data:
        return

    tm = data['tm'][data['tm']['pipeline'] == 'blue']
    mp = data['mp'][data['mp']['pipeline'] == 'blue']

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{AMBER relaxation effect: structural accuracy (TM-score) vs.\ local geometry (MolProbity). $\Delta$ = after $-$ before; $d$ = Cliff's delta effect size.}",
        r"\label{tab:amber-dual}",
        r"\begin{tabular}{llrrrr}",
        r"\toprule",
        r"Comparison & Metric & $\Delta$ & $p$-value & $d$ & Effect \\",
        r"\midrule",
    ]

    pairs = [
        ('af_unrelaxed', 'amber_af', 'AF $\\to$ AMBER(AF)'),
        ('boltz', 'amber_boltz', 'Boltz $\\to$ AMBER(Boltz)'),
    ]

    for src_before, src_after, label in pairs:
        # TM-score
        tm_b = tm[tm['source'] == src_before].groupby('target')['tmscore'].mean()
        tm_a = tm[tm['source'] == src_after].groupby('target')['tmscore'].mean()
        common_tm = tm_b.index.intersection(tm_a.index)
        if len(common_tm) > 5:
            delta_tm = (tm_a.loc[common_tm] - tm_b.loc[common_tm]).mean()
            _, p_tm = stats.wilcoxon(tm_a.loc[common_tm], tm_b.loc[common_tm])
            d_tm = cliffs_delta(tm_a.loc[common_tm].values, tm_b.loc[common_tm].values)
            p_str = f"{p_tm:.1e}" if p_tm < 0.001 else f"{p_tm:.3f}"
            lines.append(f"  {label} & TM-score & {delta_tm:+.4f} & {p_str} & {d_tm:+.3f} & {interpret_d(d_tm)} \\\\")

        # MolProbity metrics
        for metric, metric_label in [('clashscore', 'Clashscore'), ('molprobity_score', 'MP Score')]:
            mp_b = mp[mp['source'] == src_before].groupby('target')[metric].mean()
            mp_a = mp[mp['source'] == src_after].groupby('target')[metric].mean()
            common_mp = mp_b.index.intersection(mp_a.index)
            if len(common_mp) > 5:
                delta = (mp_a.loc[common_mp] - mp_b.loc[common_mp]).mean()
                _, p = stats.wilcoxon(mp_a.loc[common_mp], mp_b.loc[common_mp])
                d = cliffs_delta(mp_a.loc[common_mp].values, mp_b.loc[common_mp].values)
                p_str = f"{p:.1e}" if p < 0.001 else f"{p:.3f}"
                lines.append(f"  & {metric_label} & {delta:+.2f} & {p_str} & {d:+.3f} & {interpret_d(d)} \\\\")

        lines.append(r"  \midrule")

    # Remove last \midrule
    lines[-1] = r"\bottomrule"

    lines.extend([
        r"\end{tabular}",
        r"\end{table}",
    ])

    tex = "\n".join(lines)
    path = os.path.join(OUTDIR, "table3_amber_dual_effect.tex")
    with open(path, 'w') as f:
        f.write(tex)
    print(f"  Saved table3_amber_dual_effect.tex")
    return tex


def table4_rosetta_protocols(data):
    """Table 4: Rosetta protocol comparison."""
    if 'rosetta' not in data:
        return

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Rosetta protocol comparison (mean TM-score $\pm$ std, per-target averages across all input types).}",
        r"\label{tab:rosetta-protocols}",
        r"\begin{tabular}{lcc}",
        r"\toprule",
        r"Protocol & Blue & Green \\",
        r"\midrule",
    ]

    rt = data['rosetta']
    for proto in PROTOCOL_ORDER:
        vals = []
        for pipe in ['blue', 'green']:
            pipe_data = rt[(rt['pipeline'] == pipe) & (rt['protocol'] == proto)]
            per_target = pipe_data.groupby('target')['tmscore'].mean()
            if len(per_target) > 0:
                vals.append(f"{per_target.mean():.3f} $\\pm$ {per_target.std():.3f}")
            else:
                vals.append("---")
        label = PROTOCOL_LABELS[proto]
        lines.append(f"  {label} & {vals[0]} & {vals[1]} \\\\")

    lines.extend([
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ])

    tex = "\n".join(lines)
    path = os.path.join(OUTDIR, "table4_rosetta_protocols.tex")
    with open(path, 'w') as f:
        f.write(tex)
    print(f"  Saved table4_rosetta_protocols.tex")
    return tex


def table5_reproducibility(data):
    """Table 5: Blue-Green reproducibility summary."""
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Blue-Green pipeline reproducibility. Pearson $r$ and mean absolute error (MAE) computed on per-target averages for matching source/target pairs.}",
        r"\label{tab:reproducibility}",
        r"\begin{tabular}{lccc}",
        r"\toprule",
        r"Metric & Pearson $r$ & MAE & $n$ \\",
        r"\midrule",
    ]

    # Pre-Rosetta TM-score
    if 'tm' in data:
        tm = data['tm']
        for metric, label in [('tmscore', 'Pre-Rosetta TM'), ('rmsd', 'Pre-Rosetta RMSD')]:
            blue = tm[tm['pipeline'] == 'blue'].groupby(['source', 'target'])[metric].mean()
            green = tm[tm['pipeline'] == 'green'].groupby(['source', 'target'])[metric].mean()
            common = blue.index.intersection(green.index)
            if len(common) > 10:
                bv = blue.loc[common].values
                gv = green.loc[common].values
                mask = np.isfinite(bv) & np.isfinite(gv)
                r, _ = stats.pearsonr(bv[mask], gv[mask])
                mae = np.mean(np.abs(bv[mask] - gv[mask]))
                fmt = '.4f' if 'TM' in label else '.3f'
                lines.append(f"  {label} & {r:.4f} & {mae:{fmt}} & {mask.sum()} \\\\")

    # Rosetta TM-score
    if 'rosetta' in data:
        rt = data['rosetta']
        blue = rt[rt['pipeline'] == 'blue'].groupby(['source', 'target'])['tmscore'].mean()
        green = rt[rt['pipeline'] == 'green'].groupby(['source', 'target'])['tmscore'].mean()
        common = blue.index.intersection(green.index)
        if len(common) > 10:
            bv = blue.loc[common].values
            gv = green.loc[common].values
            mask = np.isfinite(bv) & np.isfinite(gv)
            r, _ = stats.pearsonr(bv[mask], gv[mask])
            mae = np.mean(np.abs(bv[mask] - gv[mask]))
            lines.append(f"  Rosetta TM & {r:.4f} & {mae:.4f} & {mask.sum()} \\\\")

    # MolProbity
    if 'mp' in data:
        mp = data['mp']
        for metric, label in [('clashscore', 'Clashscore'), ('molprobity_score', 'MP Score')]:
            blue = mp[mp['pipeline'] == 'blue'].groupby(['source', 'target'])[metric].mean()
            green = mp[mp['pipeline'] == 'green'].groupby(['source', 'target'])[metric].mean()
            common = blue.index.intersection(green.index)
            if len(common) > 10:
                bv = blue.loc[common].values
                gv = green.loc[common].values
                mask = np.isfinite(bv) & np.isfinite(gv)
                r, _ = stats.pearsonr(bv[mask], gv[mask])
                mae = np.mean(np.abs(bv[mask] - gv[mask]))
                lines.append(f"  {label} & {r:.4f} & {mae:.2f} & {mask.sum()} \\\\")

    lines.extend([
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ])

    tex = "\n".join(lines)
    path = os.path.join(OUTDIR, "table5_reproducibility.tex")
    with open(path, 'w') as f:
        f.write(tex)
    print(f"  Saved table5_reproducibility.tex")
    return tex


def table6_outliers(data):
    """Table 6: Notable outlier targets."""
    if 'tm' not in data:
        return

    tm = data['tm'][data['tm']['pipeline'] == 'blue']
    per_target = tm.groupby(['target', 'source'])['tmscore'].mean().reset_index()
    best_tm = per_target.groupby('target')['tmscore'].max().reset_index()
    best_tm.columns = ['target', 'best_tmscore']
    outliers = best_tm[best_tm['best_tmscore'] < 0.85].sort_values('best_tmscore')

    if len(outliers) == 0:
        print("  No outliers found")
        return

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Outlier targets with best TM-score $<$ 0.85 across all prediction methods. These targets may have large conformational changes between bound and unbound states.}",
        r"\label{tab:outliers}",
        r"\begin{tabular}{lcccc}",
        r"\toprule",
        r"Target & Best TM & AF2 TM & Boltz TM & \# Chains \\",
        r"\midrule",
    ]

    for _, row in outliers.iterrows():
        target = row['target']
        best = row['best_tmscore']
        af_row = per_target[(per_target['target'] == target) & (per_target['source'] == 'af_relaxed')]
        boltz_row = per_target[(per_target['target'] == target) & (per_target['source'] == 'boltz')]
        af_tm = f"{af_row['tmscore'].values[0]:.3f}" if len(af_row) > 0 else "---"
        boltz_tm = f"{boltz_row['tmscore'].values[0]:.3f}" if len(boltz_row) > 0 else "---"
        lines.append(f"  {target} & {best:.3f} & {af_tm} & {boltz_tm} & --- \\\\")

    lines.extend([
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ])

    tex = "\n".join(lines)
    path = os.path.join(OUTDIR, "table6_outliers.tex")
    with open(path, 'w') as f:
        f.write(tex)
    print(f"  Saved table6_outliers.tex ({len(outliers)} outliers)")
    return tex


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("Loading all metric data...")
    data = load_all()
    for key, df in data.items():
        print(f"  {key}: {len(df)} rows, {df['target'].nunique()} targets")

    print("\nGenerating LaTeX tables...")
    table1_pre_rosetta_metrics(data)
    table2_molprobity(data)
    table3_amber_dual_effect(data)
    table4_rosetta_protocols(data)
    table5_reproducibility(data)
    table6_outliers(data)

    print("\nAll LaTeX tables generated!")


if __name__ == '__main__':
    main()
