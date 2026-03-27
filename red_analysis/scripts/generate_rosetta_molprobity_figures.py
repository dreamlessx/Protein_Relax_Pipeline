#!/usr/bin/env python3
"""
generate_rosetta_molprobity_figures.py — Red Analysis Pipeline

Phase 4 figures: does Rosetta improve local geometry?

Figures:
  fig16: Pre vs Post Rosetta clashscore (paired, per source)
  fig17: Rosetta protocol MolProbity comparison (box plots)
  fig18: AMBER vs Rosetta MolProbity comparison (the key figure)
  fig19: Rosetta MolProbity by protocol × source (heatmap)
  fig20: The FULL tradeoff: TM-score cost vs MolProbity gain

COMMENT: fig20 is potentially the single most important figure in the
paper. It shows the complete picture:
  - X axis: ΔTM-score (Rosetta - pre-Rosetta) → how much accuracy is lost
  - Y axis: ΔMolProbity (Rosetta - pre-Rosetta) → how much geometry improves
  - Each point is one target
  - Color by protocol
  - AMBER shown as reference band (ΔTM≈0, ΔMP negative)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

METRICS_DIR = "/data/p_csb_meiler/agarwm5/red_analysis/metrics"
FIGDIR = "/data/p_csb_meiler/agarwm5/red_analysis/figures"

SRC_COLORS = {
    'crystal': '#999999',
    'af_relaxed': '#0072B2',
    'af_unrelaxed': '#56B4E9',
    'boltz': '#009E73',
    'amber_af': '#E69F00',
    'amber_boltz': '#D55E00',
}
PROTO_COLORS = {
    'cartesian_beta': '#0072B2',
    'cartesian_ref15': '#56B4E9',
    'dualspace_beta': '#009E73',
    'dualspace_ref15': '#CC79A7',
    'normal_beta': '#E69F00',
    'normal_ref15': '#D55E00',
}
SOURCE_LABELS = {
    'crystal': 'Crystal',
    'af_relaxed': 'AF2 (relaxed)',
    'af_unrelaxed': 'AF2 (unrelaxed)',
    'boltz': 'Boltz-1',
    'amber_af': 'AMBER(AF)',
    'amber_boltz': 'AMBER(Boltz)',
}
PROTO_LABELS = {
    'cartesian_beta': 'Cart β',
    'cartesian_ref15': 'Cart REF15',
    'dualspace_beta': 'Dual β',
    'dualspace_ref15': 'Dual REF15',
    'normal_beta': 'Norm β',
    'normal_ref15': 'Norm REF15',
}
PROTOCOL_ORDER = ['cartesian_beta', 'cartesian_ref15', 'dualspace_beta',
                  'dualspace_ref15', 'normal_beta', 'normal_ref15']


def classify_src(s):
    s = s.lower()
    if s.startswith('crystal'): return 'crystal'
    elif s.startswith('amber_af'): return 'amber_af'
    elif s.startswith('amber_boltz'): return 'amber_boltz'
    elif s.startswith('af_relaxed'): return 'af_relaxed'
    elif s.startswith('af_unrelaxed'): return 'af_unrelaxed'
    elif s.startswith('boltz'): return 'boltz'
    return 'unknown'


def save_fig(fig, name):
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(FIGDIR, f"{name}.{fmt}"),
                   dpi=300 if fmt == 'png' else None, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {name}")


def fig16_pre_post_rosetta_clash(ros_data, pre_data):
    """Paired scatter: pre vs post Rosetta clashscore by source. Both pipelines."""
    sources = ['af_relaxed', 'af_unrelaxed', 'boltz']

    for pipe in ['blue', 'green']:
        pipe_label = pipe.capitalize()
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        for ax, src in zip(axes, sources):
            pre_src = pre_data[(pre_data['pipeline'] == pipe) & (pre_data['source'] == src)]
            pre_by_target = pre_src.groupby('target')['clashscore'].mean()

            post_src = ros_data[(ros_data['pipeline'] == pipe) & (ros_data['source'] == src)]
            if len(post_src) == 0:
                ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center')
                ax.set_title(SOURCE_LABELS.get(src, src))
                continue
            post_by_target = post_src.groupby('target')['clashscore'].mean()

            common = pre_by_target.index.intersection(post_by_target.index)
            if len(common) < 3:
                ax.text(0.5, 0.5, f'n={len(common)}', transform=ax.transAxes, ha='center')
                ax.set_title(SOURCE_LABELS.get(src, src))
                continue

            x = pre_by_target.loc[common].values
            y = post_by_target.loc[common].values
            mask = np.isfinite(x) & np.isfinite(y)

            ax.scatter(x[mask], y[mask], alpha=0.6, s=30, color=SRC_COLORS.get(src, '#888'),
                      edgecolors='black', linewidths=0.3)

            max_val = max(x[mask].max(), y[mask].max()) * 1.1
            ax.plot([0, max_val], [0, max_val], 'k--', linewidth=0.5, alpha=0.3, label='No change')

            improved = (y[mask] < x[mask]).sum()
            n_total = mask.sum()
            ax.text(0.95, 0.95, f'{improved}/{n_total} improved\n'
                    f'Pre: {np.mean(x[mask]):.1f}\nPost: {np.mean(y[mask]):.1f}',
                   transform=ax.transAxes, fontsize=8, ha='right', va='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

            ax.set_xlabel('Pre-Rosetta Clashscore')
            ax.set_ylabel('Post-Rosetta Clashscore')
            ax.set_title(SOURCE_LABELS.get(src, src))
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        fig.suptitle(f'Rosetta Effect on Clashscore ({pipe_label})\n(points below diagonal = improvement)',
                    fontsize=13, fontweight='bold')
        plt.tight_layout()
        save_fig(fig, f'fig16_pre_post_rosetta_clash_{pipe}')


def fig17_rosetta_protocol_mp(ros_data):
    """Box plots of MolProbity metrics by Rosetta protocol. Both pipelines."""
    metrics = [('clashscore', 'Clashscore'), ('molprobity_score', 'MP Score'),
               ('rota_outliers', 'Rota Outliers %'), ('rama_favored', 'Rama Favored %')]

    for pipe in ['blue', 'green']:
        pipe_label = pipe.capitalize()
        pipe_data = ros_data[ros_data['pipeline'] == pipe]
        if len(pipe_data) == 0:
            print(f"  SKIP fig17: no {pipe} data")
            continue

        per_target = pipe_data.groupby(['target', 'protocol'])[
            [m[0] for m in metrics]].mean().reset_index()

        fig, axes = plt.subplots(1, 4, figsize=(20, 5))

        for ax, (metric, label) in zip(axes, metrics):
            plot_data = []
            labels = []
            colors = []
            for proto in PROTOCOL_ORDER:
                vals = per_target[per_target['protocol'] == proto][metric].dropna().values
                if len(vals) > 0:
                    plot_data.append(vals)
                    labels.append(PROTO_LABELS.get(proto, proto))
                    colors.append(PROTO_COLORS.get(proto, '#888'))

            if plot_data:
                bp = ax.boxplot(plot_data, tick_labels=labels, patch_artist=True, widths=0.6,
                               medianprops={'color': 'black', 'linewidth': 1.5})
                for patch, color in zip(bp['boxes'], colors):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)

            ax.set_ylabel(label)
            ax.tick_params(axis='x', labelsize=7, rotation=30)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        fig.suptitle(f'Rosetta Protocol MolProbity Comparison\n({pipe_label} pipeline, per-target averages)',
                    fontsize=13, fontweight='bold')
        plt.tight_layout()
        save_fig(fig, f'fig17_rosetta_protocol_mp_{pipe}')


def fig18_amber_vs_rosetta(ros_data, pre_data):
    """AMBER vs Rosetta MolProbity comparison — the key figure. Both pipelines."""
    metrics = [('clashscore', 'Clashscore'), ('molprobity_score', 'MP Score')]

    for pipe in ['blue', 'green']:
        pipe_label = pipe.capitalize()
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        for ax, (metric, label) in zip(axes, metrics):
            methods = {}

            for src in ['af_relaxed', 'af_unrelaxed', 'boltz', 'amber_af', 'amber_boltz', 'crystal']:
                pre_src = pre_data[(pre_data['pipeline'] == pipe) & (pre_data['source'] == src)]
                by_t = pre_src.groupby('target')[metric].mean()
                methods[f'pre_{src}'] = by_t

            for src in ['af_relaxed', 'af_unrelaxed', 'boltz']:
                post_src = ros_data[(ros_data['pipeline'] == pipe) & (ros_data['source'] == src)]
                if len(post_src) > 0:
                    by_t = post_src.groupby('target')[metric].mean()
                    methods[f'ros_{src}'] = by_t

            bar_data = []
            bar_labels = []
            bar_colors = []

            display_order = [
                ('pre_crystal', 'Crystal', '#999999'),
                ('pre_af_unrelaxed', 'AF2 unrlx', '#56B4E9'),
                ('pre_boltz', 'Boltz', '#009E73'),
                ('pre_af_relaxed', 'AF2 rlx', '#0072B2'),
                ('pre_amber_af', 'AMBER(AF)', '#E69F00'),
                ('pre_amber_boltz', 'AMBER(Boltz)', '#D55E00'),
                ('ros_af_relaxed', 'Ros(AF)', '#0072B2'),
                ('ros_af_unrelaxed', 'Ros(AF unrlx)', '#56B4E9'),
                ('ros_boltz', 'Ros(Boltz)', '#009E73'),
            ]

            for key, lbl, color in display_order:
                if key in methods and len(methods[key]) > 0:
                    bar_data.append(methods[key].mean())
                    bar_labels.append(lbl)
                    bar_colors.append(color)

            if not bar_data:
                continue

            x = np.arange(len(bar_data))
            bars = ax.bar(x, bar_data, color=bar_colors, edgecolor='black', linewidth=0.5,
                         alpha=0.8)

            for i, (key, _, _) in enumerate(display_order):
                if i < len(bars) and key.startswith('ros_'):
                    bars[i].set_hatch('//')

            ax.set_xticks(x)
            ax.set_xticklabels(bar_labels, fontsize=7, rotation=45, ha='right')
            ax.set_ylabel(label)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            for i, v in enumerate(bar_data):
                ax.text(i, v + 0.02*max(bar_data), f'{v:.2f}', ha='center', fontsize=7)

        fig.suptitle(f'MolProbity: Pre-Rosetta vs AMBER vs Rosetta ({pipe_label})\n(hatched = Rosetta output)',
                    fontsize=13, fontweight='bold')
        plt.tight_layout()
        save_fig(fig, f'fig18_amber_vs_rosetta_mp_{pipe}')


def fig19_protocol_source_mp_heatmap(ros_data):
    """Protocol × Source heatmap for MolProbity clashscore. Both pipelines."""
    sources = ['af_relaxed', 'af_unrelaxed', 'boltz', 'amber_af', 'amber_boltz', 'crystal']
    protocols = PROTOCOL_ORDER

    for pipe in ['blue', 'green']:
        pipe_label = pipe.capitalize()
        pipe_data = ros_data[ros_data['pipeline'] == pipe]
        if len(pipe_data) == 0:
            print(f"  SKIP fig19: no {pipe} data")
            continue

        per_target = pipe_data.groupby(['source', 'target', 'protocol'])['clashscore'].mean().reset_index()

        matrix = pd.DataFrame(index=sources, columns=protocols, dtype=float)
        for src in sources:
            for proto in protocols:
                subset = per_target[(per_target['source'] == src) & (per_target['protocol'] == proto)]
                if len(subset) > 0:
                    matrix.loc[src, proto] = subset['clashscore'].mean()

        matrix = matrix.dropna(how='all')
        if matrix.empty:
            print(f"  SKIP fig19: no data for heatmap ({pipe})")
            continue

        fig, ax = plt.subplots(figsize=(10, 6))
        matrix_float = matrix.astype(float)
        sns.heatmap(matrix_float, annot=True, fmt='.2f', cmap='RdYlGn_r',
                   ax=ax, linewidths=0.5, cbar_kws={'label': 'Clashscore'},
                   xticklabels=[PROTO_LABELS.get(p, p) for p in matrix.columns],
                   yticklabels=[SOURCE_LABELS.get(s, s) for s in matrix.index])
        ax.set_title(f'Rosetta Clashscore by Protocol × Source\n({pipe_label} pipeline)',
                    fontsize=13, fontweight='bold')
        plt.tight_layout()
        save_fig(fig, f'fig19_rosetta_mp_heatmap_{pipe}')


def fig20_tradeoff(ros_mp, ros_tm, pre_mp, pre_tm):
    """THE FIGURE: TM-score cost vs MolProbity gain per target. Both pipelines."""
    for pipe in ['blue', 'green']:
        pipe_label = pipe.capitalize()
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        for ax, src in zip(axes, ['af_relaxed', 'af_unrelaxed']):
            # Pre-Rosetta TM and MP per target
            pre_tm_src = pre_tm[(pre_tm['pipeline'] == pipe) & (pre_tm['source'] == src)]
            pre_tm_by_t = pre_tm_src.groupby('target')['tmscore'].mean()

            pre_mp_src = pre_mp[(pre_mp['pipeline'] == pipe) & (pre_mp['source'] == src)]
            pre_mp_by_t = pre_mp_src.groupby('target')['clashscore'].mean()

            # Post-Rosetta (mean across protocols)
            post_tm_src = ros_tm[(ros_tm['pipeline'] == pipe) & (ros_tm['source'] == src)]
            if len(post_tm_src) == 0:
                ax.text(0.5, 0.5, 'No Rosetta TM data', transform=ax.transAxes, ha='center')
                ax.set_title(SOURCE_LABELS.get(src, src))
                continue
            post_tm_by_t = post_tm_src.groupby('target')['tmscore'].mean()

            post_mp_src = ros_mp[(ros_mp['pipeline'] == pipe) & (ros_mp['source'] == src)]
            if len(post_mp_src) == 0:
                ax.text(0.5, 0.5, 'No Rosetta MP data', transform=ax.transAxes, ha='center')
                ax.set_title(SOURCE_LABELS.get(src, src))
                continue
            post_mp_by_t = post_mp_src.groupby('target')['clashscore'].mean()

            # Find targets with all 4 values
            common = (pre_tm_by_t.index
                      .intersection(pre_mp_by_t.index)
                      .intersection(post_tm_by_t.index)
                      .intersection(post_mp_by_t.index))

            if len(common) < 3:
                ax.text(0.5, 0.5, f'n={len(common)} targets', transform=ax.transAxes, ha='center')
                ax.set_title(SOURCE_LABELS.get(src, src))
                continue

            delta_tm = post_tm_by_t.loc[common] - pre_tm_by_t.loc[common]
            delta_mp = post_mp_by_t.loc[common] - pre_mp_by_t.loc[common]

            mask = np.isfinite(delta_tm.values) & np.isfinite(delta_mp.values)

            ax.scatter(delta_tm.values[mask], delta_mp.values[mask],
                      alpha=0.6, s=40, color=SRC_COLORS.get(src, '#888'),
                      edgecolors='black', linewidths=0.3, label='Rosetta')

            # Reference lines
            ax.axhline(y=0, color='gray', linewidth=0.5, linestyle='-', alpha=0.5)
            ax.axvline(x=0, color='gray', linewidth=0.5, linestyle='-', alpha=0.5)

            # AMBER reference: add a star showing AMBER's effect (ΔTM≈0, ΔMP<0)
            amber_src = 'amber_af' if src == 'af_unrelaxed' else 'amber_af'
            pre_src_for_amber = 'af_unrelaxed' if src == 'af_unrelaxed' else 'af_relaxed'
            amber_tm = pre_tm[(pre_tm['pipeline'] == pipe) & (pre_tm['source'] == amber_src)]
            amber_tm_by_t = amber_tm.groupby('target')['tmscore'].mean()
            amber_mp = pre_mp[(pre_mp['pipeline'] == pipe) & (pre_mp['source'] == amber_src)]
            amber_mp_by_t = amber_mp.groupby('target')['clashscore'].mean()

            pre_for_amber_tm = pre_tm[(pre_tm['pipeline'] == pipe) & (pre_tm['source'] == pre_src_for_amber)]
            pre_for_amber_tm_by_t = pre_for_amber_tm.groupby('target')['tmscore'].mean()
            pre_for_amber_mp = pre_mp[(pre_mp['pipeline'] == pipe) & (pre_mp['source'] == pre_src_for_amber)]
            pre_for_amber_mp_by_t = pre_for_amber_mp.groupby('target')['clashscore'].mean()

            amber_common = (amber_tm_by_t.index
                           .intersection(amber_mp_by_t.index)
                           .intersection(pre_for_amber_tm_by_t.index)
                           .intersection(pre_for_amber_mp_by_t.index))

            if len(amber_common) >= 3:
                amber_dtm = amber_tm_by_t.loc[amber_common] - pre_for_amber_tm_by_t.loc[amber_common]
                amber_dmp = amber_mp_by_t.loc[amber_common] - pre_for_amber_mp_by_t.loc[amber_common]
                ax.scatter(amber_dtm.mean(), amber_dmp.mean(), marker='*', s=200,
                          color='#E69F00', edgecolors='black', linewidths=1, zorder=5,
                          label=f'AMBER mean')

            # Quadrant labels
            ax.text(0.02, 0.02, 'Better both\n(ideal)', transform=ax.transAxes,
                   fontsize=7, color='green', alpha=0.5, va='bottom')
            ax.text(0.98, 0.02, 'Better MP,\nworse TM', transform=ax.transAxes,
                   fontsize=7, color='blue', alpha=0.5, ha='right', va='bottom')

            ax.set_xlabel('ΔTM-score (post - pre)')
            ax.set_ylabel('ΔClashscore (post - pre)')
            ax.set_title(f'{SOURCE_LABELS.get(src, src)} → Rosetta\n(n={mask.sum()} targets)')
            ax.legend(fontsize=8, loc='upper right')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        fig.suptitle(f'Rosetta Tradeoff: TM-score Cost vs Clashscore Improvement ({pipe_label})\n'
                    '(bottom-left = ideal: better TM AND better clashscore)',
                    fontsize=13, fontweight='bold')
        plt.tight_layout()
        save_fig(fig, f'fig20_rosetta_tradeoff_{pipe}')


def main():
    os.makedirs(FIGDIR, exist_ok=True)

    print("Loading data...")

    # Rosetta MolProbity (Phase 4)
    ros_mp_path = os.path.join(METRICS_DIR, "combined_rosetta_molprobity.tsv")
    if not os.path.exists(ros_mp_path):
        # Try aggregating on the fly
        import glob as gl
        files = sorted(gl.glob(os.path.join(METRICS_DIR, "rosetta_molprobity_*.tsv")))
        if not files:
            print("ERROR: No Rosetta MolProbity data found")
            return
        dfs = []
        for f in files:
            try:
                df = pd.read_csv(f, sep='\t')
                if len(df) > 0:
                    dfs.append(df)
            except:
                pass
        ros_mp = pd.concat(dfs, ignore_index=True)
        ros_mp['source'] = ros_mp['src_type'].apply(classify_src)
    else:
        ros_mp = pd.read_csv(ros_mp_path, sep='\t')
        if 'source' not in ros_mp.columns:
            ros_mp['source'] = ros_mp['src_type'].apply(classify_src)

    for col in ['clashscore', 'molprobity_score', 'rota_outliers', 'rama_favored',
                'rms_bonds', 'rms_angles']:
        if col in ros_mp.columns:
            ros_mp[col] = pd.to_numeric(ros_mp[col], errors='coerce')
    print(f"  Rosetta MP: {len(ros_mp)} rows, {ros_mp['target'].nunique()} targets")

    # Pre-Rosetta MolProbity (Phase 3)
    pre_mp = pd.read_csv(os.path.join(METRICS_DIR, "combined_molprobity.tsv"), sep='\t')
    for col in ['clashscore', 'molprobity_score', 'rota_outliers', 'rama_favored']:
        pre_mp[col] = pd.to_numeric(pre_mp[col], errors='coerce')
    print(f"  Pre-Rosetta MP: {len(pre_mp)} rows")

    # Rosetta TM-score (Phase 2)
    ros_tm_path = os.path.join(METRICS_DIR, "combined_rosetta_tmscore.tsv")
    ros_tm = None
    if os.path.exists(ros_tm_path):
        ros_tm = pd.read_csv(ros_tm_path, sep='\t')
        ros_tm['tmscore'] = pd.to_numeric(ros_tm['tmscore'], errors='coerce')
        if 'source' not in ros_tm.columns:
            ros_tm['source'] = ros_tm['src_type'].apply(classify_src)
        print(f"  Rosetta TM: {len(ros_tm)} rows, {ros_tm['target'].nunique()} targets")

    # Pre-Rosetta TM-score (Phase 1)
    pre_tm = pd.read_csv(os.path.join(METRICS_DIR, "combined_tmscore.tsv"), sep='\t')
    pre_tm['tmscore'] = pd.to_numeric(pre_tm['tmscore'], errors='coerce')
    print(f"  Pre-Rosetta TM: {len(pre_tm)} rows")

    print("\nGenerating Phase 4 figures...")
    fig16_pre_post_rosetta_clash(ros_mp, pre_mp)
    fig17_rosetta_protocol_mp(ros_mp)
    fig18_amber_vs_rosetta(ros_mp, pre_mp)
    fig19_protocol_source_mp_heatmap(ros_mp)
    if ros_tm is not None:
        fig20_tradeoff(ros_mp, ros_tm, pre_mp, pre_tm)
    else:
        print("  SKIP fig20: no Rosetta TM data")

    print("\nAll Phase 4 figures generated!")


if __name__ == '__main__':
    main()
