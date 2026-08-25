"""
06_control_composition_plots.py

Generate the four control-well composition plots (relative abundance,
stacked bar) from community_strain_abundance.csv:
  1. in vitro t0 wells (S1, S1+C2)
  2. in vitro monoculture wells
  3. in planta t0 wells (t0 1, t0 1 L, t0 1C, t0 1C L)
  4. in planta monoculture wells

Design notes:
- Taxon coloring/order is fixed (TAXON_ORDER below) so it's identical across
  all four plots and across reruns.
- The ambiguous "59+p25C2" call_label from 03a is kept as its own category
  rather than folded into either Leaf59 or p25C2 -- collapsing it would hide
  a real signal (see the Leaf59/p25C2 R2-disambiguation note in the project
  log).
- in_planta_master_metadata.csv may or may not have a STRAIN column
  depending on when it was built relative to the 00a fix that adds it for
  monoculture rows. This script checks for STRAIN and falls back to
  normalizing growth_plate_raw itself if it's missing, so it doesn't hard-fail
  on an older metadata file -- but it prints a warning either way, since an
  old file elsewhere in the metadata could mean other things are stale too.

Usage:
    python 06_control_composition_plots.py \
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

_cmap = plt.get_cmap('tab20')
COLORS = {t: _cmap(i / len(TAXON_ORDER)) for i, t in enumerate(TAXON_ORDER)}
COLORS['Unassigned'] = '#B0B0B0'
COLORS['Leaf59+p25C2 (ambiguous)'] = '#444444'


def relabel(s):
    """Map a raw strain_or_group / STRAIN / growth_plate_raw value to a
    display taxon name matching TAXON_ORDER."""
    s = str(s).strip()
    if s.lower() == 'unassigned':
        return 'Unassigned'
    if s == 'p25C2':
        return 'p25C2'
    if s == '59+p25C2':
        return 'Leaf59+p25C2 (ambiguous)'
    if re.match(r'^\d+\.0$', s):  # stray float-formatted strain id, e.g. '21.0'
        s = s[:-2]
    return f'Leaf{s}'


def well_sort_key(w):
    """Sort by column then row (e.g. A09, B09, ... H09, A10, B10, ...),
    matching the well ordering used in the original reference plots."""
    row = w[0]
    col = int(w[1:])
    return (col, row)


def make_stacked_plot(comm, uid_to_well, group_order, group_of_uid, title,
                       out_path, well_label_map=None, figsize=(11, 5)):
    """comm: long-format df with sample_uid, taxon, R2_count, restricted to
    the sample_uids of interest.
    group_order: facet order (e.g. ['S1', 'S1+C2']).
    group_of_uid: dict sample_uid -> facet label.
    uid_to_well: dict sample_uid -> well label for x-axis.
    well_label_map: optional fn(well) -> display label (e.g. adds strain name).
    """
    missing_uids = set(group_of_uid) - set(comm['sample_uid'])
    if missing_uids:
        raise RuntimeError(
            f"{len(missing_uids)} sample_uids expected for '{title}' have no rows "
            f"in community_strain_abundance.csv: {sorted(missing_uids)}"
        )

    tot = comm.groupby('sample_uid')['R2_count'].transform('sum')
    comm = comm.assign(rel=comm['R2_count'] / tot)

    pivot = comm.pivot_table(index='sample_uid', columns='taxon', values='rel', fill_value=0)
    for t in TAXON_ORDER:
        if t not in pivot.columns:
            pivot[t] = 0.0
    pivot = pivot[TAXON_ORDER]

    fig, axes = plt.subplots(1, len(group_order), figsize=figsize, sharey=True)
    if len(group_order) == 1:
        axes = [axes]

    for ax, grp in zip(axes, group_order):
        uids_in_grp = [u for u in pivot.index if group_of_uid.get(u) == grp]
        uids_in_grp = sorted(uids_in_grp, key=lambda u: well_sort_key(uid_to_well[u]))
        sub = pivot.loc[uids_in_grp] * 100
        wells = [uid_to_well[u] for u in uids_in_grp]

        bottom = np.zeros(len(sub))
        x = np.arange(len(sub))
        for t in TAXON_ORDER:
            vals = sub[t].values
            ax.bar(x, vals, bottom=bottom, color=COLORS[t], width=0.85,
                   edgecolor='white', linewidth=0.3)
            bottom += vals

        ax.set_xticks(x)
        labels = [well_label_map(w) if well_label_map else w for w in wells]
        ax.set_xticklabels(labels, rotation=90, fontsize=7)
        ax.set_title(grp, fontsize=10)
        ax.set_ylim(0, 100)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    axes[0].set_ylabel('Relative abundance (%)')
    handles = [mpatches.Patch(color=COLORS[t], label=t) for t in TAXON_ORDER]
    fig.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5),
               fontsize=8, title='Taxon', frameon=False)
    fig.suptitle(title, fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f'  wrote {out_path}')


def plot_in_vitro_t0(comm, iv, out_dir):
    sub_meta = iv[iv['sample_type'].isin(['control_S1_t0', 'control_S1C_t0'])]
    uid_to_well = dict(zip(sub_meta['sample_uid'], sub_meta['well']))
    group_of_uid = dict(zip(
        sub_meta['sample_uid'],
        sub_meta['sample_type'].map({'control_S1_t0': 'S1', 'control_S1C_t0': 'S1+C2'})
    ))
    df = comm[comm['sample_uid'].isin(sub_meta['sample_uid'])]
    make_stacked_plot(
        df, uid_to_well, ['S1', 'S1+C2'], group_of_uid,
        'Relative abundance: in vitro t0 wells (leftover starting inoculum)',
        os.path.join(out_dir, '1_invitro_t0.pdf'), figsize=(6, 5),
    )


def plot_in_vitro_monoculture(comm, iv, out_dir):
    sub_meta = iv[iv['sample_type'] == 'monoculture']
    uid_to_well = dict(zip(sub_meta['sample_uid'], sub_meta['well']))
    group_of_uid = {u: 'Monocultures' for u in sub_meta['sample_uid']}
    df = comm[comm['sample_uid'].isin(sub_meta['sample_uid'])]
    strain_label = dict(zip(sub_meta['well'], sub_meta['STRAIN'].map(relabel)))
    make_stacked_plot(
        df, uid_to_well, ['Monocultures'], group_of_uid,
        'Relative abundance: in vitro individual-member (monoculture) control wells\n'
        '(expect ~100% one taxon per well)',
        os.path.join(out_dir, '2_invitro_monoculture.pdf'),
        well_label_map=lambda w: f'{w}\n({strain_label.get(w, "?")})',
        figsize=(9, 5),
    )


def plot_in_planta_t0(comm, ip, out_dir):
    sub_meta = ip[ip['sample_type'].isin(
        ['control_t01', 'control_t01C', 'control_t01CL', 'control_t01L']
    )]
    uid_to_well = dict(zip(sub_meta['sample_uid'], sub_meta['collection_well']))
    group_label_map = {'control_t01': 't0 1', 'control_t01L': 't0 1 L',
                        'control_t01C': 't0 1C', 'control_t01CL': 't0 1C L'}
    group_of_uid = dict(zip(sub_meta['sample_uid'], sub_meta['sample_type'].map(group_label_map)))
    df = comm[comm['sample_uid'].isin(sub_meta['sample_uid'])]
    make_stacked_plot(
        df, uid_to_well, ['t0 1', 't0 1 L', 't0 1C', 't0 1C L'], group_of_uid,
        'Relative abundance: in planta t0 wells (leftover starting inoculum)',
        os.path.join(out_dir, '3_inplanta_t0.pdf'), figsize=(11, 5),
    )


def plot_in_planta_monoculture(comm, ip, out_dir):
    sub_meta = ip[ip['sample_type'] == 'monoculture'].copy()
    uid_to_well = dict(zip(sub_meta['sample_uid'], sub_meta['collection_well']))
    group_of_uid = {u: 'Monocultures' for u in sub_meta['sample_uid']}
    df = comm[comm['sample_uid'].isin(sub_meta['sample_uid'])]

    if 'STRAIN' in sub_meta.columns and sub_meta['STRAIN'].notna().any():
        strain_label = dict(zip(sub_meta['collection_well'], sub_meta['STRAIN'].map(relabel)))
    else:
        print("  WARNING: in_planta_master_metadata.csv has no usable STRAIN column "
              "for monoculture rows -- falling back to growth_plate_raw. This metadata "
              "file predates the 00a STRAIN fix; consider rebuilding it.")
        strain_label = dict(zip(sub_meta['collection_well'], sub_meta['growth_plate_raw'].map(relabel)))

    make_stacked_plot(
        df, uid_to_well, ['Monocultures'], group_of_uid,
        'Relative abundance: in planta individual-member (monoculture) control wells\n'
        '(expect ~100% one taxon per well)',
        os.path.join(out_dir, '4_inplanta_monoculture.pdf'),
        well_label_map=lambda w: f'{w}\n({strain_label.get(w, "?")})',
        figsize=(9, 5),
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--community', required=True, help='community_strain_abundance.csv')
    p.add_argument('--in-vitro-metadata', required=True, help='in_vitro_master_metadata.csv')
    p.add_argument('--in-planta-metadata', required=True, help='in_planta_master_metadata.csv')
    p.add_argument('--output-dir', required=True)
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    comm = pd.read_csv(args.community)
    comm['taxon'] = comm['strain_or_group'].astype(str).map(relabel)

    iv = pd.read_csv(args.in_vitro_metadata, dtype={'STRAIN': str})
    ip = pd.read_csv(args.in_planta_metadata,
                      dtype={'STRAIN': str} if 'STRAIN' in pd.read_csv(args.in_planta_metadata, nrows=0).columns else None)

    print('Plot 1/4: in vitro t0 wells')
    plot_in_vitro_t0(comm, iv, args.output_dir)

    print('Plot 2/4: in vitro monoculture wells')
    plot_in_vitro_monoculture(comm, iv, args.output_dir)

    print('Plot 3/4: in planta t0 wells')
    plot_in_planta_t0(comm, ip, args.output_dir)

    print('Plot 4/4: in planta monoculture wells')
    plot_in_planta_monoculture(comm, ip, args.output_dir)

    print(f'\nAll plots written to {args.output_dir}/')


if __name__ == '__main__':
    main()
