"""
Per-replicate composition plots faceted by (treatment/strain-or-community) x
(timepoint), matching the layout of the existing reference figure
(plot_replicates_by_treatment.pdf): each cell is a small multiples panel with
one bar per replicate plus a bold "MEAN" bar on the right.

Two variants, controlled by --mode:
  relative : each bar normalized to 100% (this is what the reference PDF shows)
  depth    : bars keep their raw read count height, so low-depth wells read
             visibly shorter -- MEAN bar is the elementwise mean raw count per
             taxon across replicates (so it still sums to something meaningful:
             the average total depth for that group, broken down by taxon).

Produces one plot per (system, mode) combination:
  in vitro            -> facet cols = timepoint (6H/12H/18H/24H)
                          facet rows = STRAIN (S1, S1+C2) x TREATMENT (Buffer, delTFA, Tailocin)
  in planta trial 1/2 -> facet cols = dpi
                          facet rows = community (Buffer, SynCom1, SynCom1C) x treatment (Buf, TFA, Tail)

Only sample_type == 'experimental' rows are included, and (for in planta)
only truly_inoculated == True, since dropped/not-actually-inoculated wells
would misrepresent the treatment groups they're nominally labeled with.

Usage:
    python 04b_composition_by_treatment_timepoint.py \
        --community community_strain_abundance.csv \
        --in-vitro-metadata in_vitro_master_metadata.csv \
        --in-planta-metadata in_planta_master_metadata.csv \
        --output-dir composition_plots
"""

import argparse
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

TAXON_ORDER = ['Leaf13', 'Leaf132', 'Leaf179', 'Leaf201', 'Leaf21', 'Leaf233',
               'Leaf234', 'Leaf391', 'Leaf59', 'Leaf84', 'p25C2',
               'Leaf59+p25C2 (ambiguous)', 'Unassigned']

# Colors come from taxon_colors.py, shared with 08 so a strain is the same
# color in every figure: NatParksPalettes "DeathValley" (n = 11, continuous)
# for the 11 members, neutral grays for Unassigned / ambiguous reads.
from taxon_colors import taxon_color
COLORS = {t: taxon_color(t) for t in TAXON_ORDER}

IN_VITRO_TREATMENT_ORDER = ['Buffer', 'delTFA', 'Tailocin']
IN_PLANTA_TREATMENT_ORDER = ['Buf', 'TFA', 'Tail']


def relabel(s):
    s = str(s).strip()
    if s.lower() == 'unassigned':
        return 'Unassigned'
    if s == 'p25C2':
        return 'p25C2'
    if s == '59+p25C2':
        return 'Leaf59+p25C2 (ambiguous)'
    if re.match(r'^\d+\.0$', s):
        s = s[:-2]
    return f'Leaf{s}'


def build_pivots(comm, sample_uids):
    """Returns (raw_pivot, rel_pivot): sample_uid x taxon, in TAXON_ORDER,
    restricted to sample_uids. raw = counts, rel = proportions summing to 1."""
    sub = comm[comm['sample_uid'].isin(sample_uids)]
    raw = sub.pivot_table(index='sample_uid', columns='taxon', values='R2_count',
                           aggfunc='sum', fill_value=0)
    for t in TAXON_ORDER:
        if t not in raw.columns:
            raw[t] = 0
    raw = raw[TAXON_ORDER]
    raw = raw.reindex(sample_uids, fill_value=0)
    tot = raw.sum(axis=1).replace(0, np.nan)
    rel = raw.div(tot, axis=0).fillna(0)
    return raw, rel


def plot_grid(comm, cell_samples, row_labels, col_labels, mode, title, out_path,
              panel_w=1.6, panel_h=1.5, label_col_w=2.2):
    """
    cell_samples: dict (row_idx, col_idx) -> list of sample_uid (replicate order)
    row_labels: list of str, one per row (may be multi-line, joined with '\\n')
    col_labels: list of str, one per column
    mode: 'relative' or 'depth'
    """
    n_rows, n_cols = len(row_labels), len(col_labels)
    fig_w = label_col_w + panel_w * n_cols
    fig_h = 0.6 + panel_h * n_rows
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h), squeeze=False)

    # shared y max for depth mode so panels are visually comparable;
    # relative mode is always exactly 0-100 by construction.
    shared_max = 0.0
    if mode == 'depth':
        for uids in cell_samples.values():
            if not uids:
                continue
            raw, _ = build_pivots(comm, uids)
            mean_total = raw.mean(axis=0).sum()
            cell_max = max(raw.sum(axis=1).max(), mean_total)
            shared_max = max(shared_max, cell_max)

    for r in range(n_rows):
        for c in range(n_cols):
            ax = axes[r][c]
            uids = cell_samples.get((r, c), [])

            if not uids:
                ax.axis('off')
                continue

            raw, rel = build_pivots(comm, uids)
            n_reps = len(uids)
            x = np.arange(n_reps + 1)  # +1 for MEAN bar

            if mode == 'relative':
                bar_vals = rel.values * 100         # n_reps x n_taxa
                mean_vals = (rel.mean(axis=0) * 100).values
                ax.set_ylim(0, 100)
            else:
                bar_vals = raw.values
                mean_vals = raw.mean(axis=0).values
                ax.set_ylim(0, shared_max * 1.05 if shared_max > 0 else 1)

            bottom = np.zeros(n_reps)
            for ti, t in enumerate(TAXON_ORDER):
                vals = bar_vals[:, ti]
                ax.bar(x[:n_reps], vals, bottom=bottom, color=COLORS[t],
                       width=0.85, edgecolor='white', linewidth=0.2)
                bottom += vals

            bottom_mean = 0.0
            for ti, t in enumerate(TAXON_ORDER):
                v = mean_vals[ti]
                ax.bar(x[n_reps], v, bottom=bottom_mean, color=COLORS[t],
                       width=0.85, edgecolor='black', linewidth=0.8)
                bottom_mean += v

            xticks = [f'Rep {i+1}' for i in range(n_reps)] + ['MEAN']
            ax.set_xticks(x)
            ax.set_xticklabels(xticks, rotation=90, fontsize=5)
            for lbl, tick in zip(ax.get_xticklabels(), ax.get_xticklabels()):
                pass
            ax.get_xticklabels()[-1].set_fontweight('bold')
            ax.tick_params(axis='y', labelsize=6)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            if r == 0:
                ax.set_title(col_labels[c], fontsize=8)
            if c == 0:
                ax.set_ylabel(row_labels[r], fontsize=6.5, rotation=0,
                               ha='right', va='center', labelpad=8)

    ylab = 'Relative abundance (%)' if mode == 'relative' else 'Read count'
    fig.text(0.005, 0.5, ylab, va='center', rotation='vertical', fontsize=10)

    handles = [mpatches.Patch(color=COLORS[t], label=t) for t in TAXON_ORDER]
    fig.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5),
               fontsize=8, title='Taxon', frameon=False)
    fig.suptitle(title, fontsize=13, y=1.0)
    fig.tight_layout(rect=[0.03, 0, 1, 0.98])
    fig.savefig(out_path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f'  wrote {out_path}')


def check_coverage(comm, meta_subset, system_label):
    missing = meta_subset[~meta_subset['sample_uid'].isin(comm['sample_uid'])]
    if len(missing) > 0:
        print(f'  WARNING [{system_label}]: {len(missing)} experimental wells have no rows '
              f'in community_strain_abundance.csv (treated as dropped/missing replicates, '
              f'not an error): {sorted(missing["sample_uid"].tolist())}')


def plot_in_vitro(comm, iv, mode, out_dir):
    exp = iv[iv['sample_type'] == 'experimental'].copy()
    check_coverage(comm, exp, 'in vitro')
    exp = exp[exp['sample_uid'].isin(comm['sample_uid'])]

    strains = ['S1', 'S1+C2']
    treatments = IN_VITRO_TREATMENT_ORDER
    timepoints = ['6H', '12H', '18H', '24H']

    row_defs = [(s, t) for s in strains for t in treatments]
    row_labels = [f'treatment: {t}\nstrain: {s}' for s, t in row_defs]

    cell_samples = {}
    for r, (s, t) in enumerate(row_defs):
        for c, tp in enumerate(timepoints):
            sub = exp[(exp['STRAIN'] == s) & (exp['TREATMENT'] == t) & (exp['timepoint'] == tp)]
            cell_samples[(r, c)] = sorted(sub['sample_uid'].tolist())

    col_labels = [f'time_point: {tp.lower()}' for tp in timepoints]
    mode_label = 'relative abundance' if mode == 'relative' else 'read depth'
    title = f'In vitro: {mode_label} per replicate by treatment and time point'
    out_path = os.path.join(out_dir, f'invitro_{mode}.pdf')
    plot_grid(comm, cell_samples, row_labels, col_labels, mode, title, out_path)


def plot_in_planta(comm, ip, trial, mode, out_dir):
    exp = ip[(ip['sample_type'] == 'experimental') & (ip['trial'] == trial)].copy()
    exp = exp[exp['truly_inoculated'] == True]  # noqa: E712
    check_coverage(comm, exp, f'in planta trial {int(trial)}')
    exp = exp[exp['sample_uid'].isin(comm['sample_uid'])]

    communities = ['Buffer', 'SynCom1', 'SynCom1C']
    treatments = IN_PLANTA_TREATMENT_ORDER
    dpis = sorted(exp['dpi'].dropna().unique())

    row_defs = [(comm_, t) for comm_ in communities for t in treatments]
    row_labels = [f'treatment: {t}\ncommunity: {comm_}' for comm_, t in row_defs]

    cell_samples = {}
    for r, (comm_, t) in enumerate(row_defs):
        for c, dpi in enumerate(dpis):
            sub = exp[(exp['community'] == comm_) & (exp['treatment'] == t) & (exp['dpi'] == dpi)]
            cell_samples[(r, c)] = sorted(sub['sample_uid'].tolist())

    col_labels = [f'dpi: {dpi:g}' for dpi in dpis]
    mode_label = 'relative abundance' if mode == 'relative' else 'read depth'
    title = f'In planta trial {int(trial)}: {mode_label} per replicate by treatment and dpi'
    out_path = os.path.join(out_dir, f'inplanta_trial{int(trial)}_{mode}.pdf')
    plot_grid(comm, cell_samples, row_labels, col_labels, mode, title, out_path)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--community', required=True)
    p.add_argument('--in-vitro-metadata', required=True)
    p.add_argument('--in-planta-metadata', required=True)
    p.add_argument('--output-dir', required=True)
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    comm = pd.read_csv(args.community)
    comm['taxon'] = comm['strain_or_group'].astype(str).map(relabel)

    iv = pd.read_csv(args.in_vitro_metadata, dtype={'STRAIN': str})
    ip = pd.read_csv(args.in_planta_metadata)
    ip['trial'] = pd.to_numeric(ip['trial'], errors='coerce')

    for mode in ['relative', 'depth']:
        print(f'--- mode: {mode} ---')
        print('in vitro')
        plot_in_vitro(comm, iv, mode, args.output_dir)
        for trial in [1.0, 2.0]:
            print(f'in planta trial {int(trial)}')
            plot_in_planta(comm, ip, trial, mode, args.output_dir)

    print(f'\nAll plots written to {args.output_dir}/')


if __name__ == '__main__':
    main()
