"""
10_taxon_boxplot.py

Relative-abundance box plot for a single "headline" taxon (default: p25C2)
by treatment, faceted by community background x timepoint/dpi -- the
reader-friendly companion to the CLR effect-size forest plot from
08_tailocin_effect_plot.py.

DESIGN PRINCIPLE -- plot on relative abundance, test on CLR:
Relative abundance is what a reader intuitively understands ("what fraction
of the community is p25C2"), but it's compositional data (bounded, skewed,
not independent across taxa) and is NOT the right scale to run a fresh
significance test on. This script does NOT compute a new test from the
plotted numbers -- it reads the q-values straight out of the effects CSV
already produced by 08_tailocin_effect_plot.py (Welch's t-test on CLR values,
BH-FDR corrected across that script's full test family). That keeps the
box plot's significance annotation making the exact same statistical claim
as the forest plot, just visualized differently, rather than risking the two
disagreeing because they're secretly two different tests.

Shows all three treatments (Buffer/delTFA/Tailocin or Buf/TFA/Tail) per
panel, with two nested significance brackets: Tailocin vs Buffer (outer) and
Tailocin vs delTFA/TFA (inner) -- matching the two contrasts 08 already
computes. Faceted by timepoint (no separate "zoomed" panel -- full
time-course only, per request).

BACKGROUND SCOPE: defaults to S1+C2 (in vitro) and SynCom1C (in planta) only
-- the communities where the taxon is actually a deliberately-added member.
S1 and SynCom1 are the same kind of background in both systems (no real
p25C2 present, just contamination-level reads), so showing them alongside
their C2/C counterparts is a redundant near-null comparison for THIS taxon --
dropped by default, restorable via --invitro-backgrounds / --inplanta-backgrounds
for a taxon that's a real member of those communities instead.

Y-axis is log-scale by default (relative abundance commonly spans multiple
orders of magnitude across timepoints/treatments -- a linear axis gets
crushed by that range, and it pairs naturally with CLR being log-ratio in
the first place). Wells where the taxon has ZERO reads can't be placed on a
log axis -- they're plotted as hollow "not detected" markers at a fixed
visual floor just below the lowest real value, not silently dropped, so a
few big lift in y-axis, no zero point disappears from the row.

Usage:
    python 10_taxon_boxplot.py \
        --community community_strain_abundance.csv \
        --in-vitro-metadata in_vitro_master_metadata.csv \
        --in-planta-metadata in_planta_master_metadata.csv \
        --invitro-effects-csv effect_plots/invitro_tailocin_effects.csv \
        --inplanta-effects-csv effect_plots/inplanta_trial2_tailocin_effects.csv \
        --output-dir boxplots \
        [--taxon p25C2]
"""

import argparse
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

TAXON_COLOR = '#7a3ec2'


def numeric_sort_key(val):
    s = str(val)
    m = re.match(r'^(\d+(?:\.\d+)?)', s)
    return float(m.group(1)) if m else float('inf')


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


def significance_stars(q):
    if pd.isna(q):
        return 'n/a'
    if q < 0.001:
        return '***'
    if q < 0.01:
        return '**'
    if q < 0.05:
        return '*'
    return 'ns'


def compute_relative_abundance(comm, sample_uids, taxon):
    """Returns a df: sample_uid, rel_abundance_pct, raw_count, detected (bool).
    Denominator is each sample's TOTAL reads (all taxa incl. unassigned) --
    matches the definition used throughout this pipeline's earlier
    diagnostics, not just the CLR-parts total."""
    sub = comm[comm['sample_uid'].isin(sample_uids)].copy()
    sub['taxon'] = sub['strain_or_group'].astype(str).map(relabel)
    total = sub.groupby('sample_uid')['R2_count'].sum().rename('total_reads')
    taxon_count = sub[sub['taxon'] == taxon].groupby('sample_uid')['R2_count'].sum().rename('raw_count')
    out = pd.concat([total, taxon_count], axis=1)
    out['raw_count'] = out['raw_count'].fillna(0)
    out = out.reindex(sample_uids, fill_value=0)
    # samples with genuinely zero total reads shouldn't happen post-QC, but
    # guard against div-by-zero producing inf rather than a loud failure
    out['total_reads'] = out['total_reads'].replace(0, np.nan)
    out['rel_abundance_pct'] = (out['raw_count'] / out['total_reads'] * 100)
    out['detected'] = out['raw_count'] > 0
    return out.reset_index().rename(columns={'index': 'sample_uid'})


def draw_significance_bracket(ax, x1, x2, y, height, text):
    """Draws a bracket from x1 to x2 at height y (in data coords), with a
    centered label. height controls the tick length down at each end."""
    ax.plot([x1, x1, x2, x2], [y - height, y, y, y - height], color='black', linewidth=0.9)
    ax.text((x1 + x2) / 2, y, text, ha='center', va='bottom', fontsize=8)


def plot_taxon_boxplot(rel_df, effects_df, taxon, background_col, time_col, treatment_col,
                        backgrounds, treatment_order, contrasts, title, out_path,
                        log_scale=True, figsize_per_panel=(2.6, 3.4)):
    """
    rel_df: sample_uid, rel_abundance_pct, raw_count, detected, plus
        background_col/time_col/treatment_col merged in.
    effects_df: the 08-script output (background, timepoint, taxon, contrast,
        q_value, significant, ...) -- used ONLY for the bracket annotations,
        never for the plotted values.
    contrasts: list of (treatment_a, treatment_b) tuples, drawn as nested
        brackets in the order given (first = innermost/lowest).
    """
    timepoints = sorted(rel_df[time_col].dropna().unique(), key=numeric_sort_key)
    n_rows, n_cols = len(backgrounds), len(timepoints)
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(figsize_per_panel[0] * n_cols, figsize_per_panel[1] * n_rows),
                              squeeze=False)

    # shared y-range across every timepoint panel: compute once from the
    # WHOLE dataset being plotted (not per-panel), so panels are directly
    # comparable at a glance and brackets sit at a consistent height instead
    # of each panel silently rescaling to its own local max.
    positive_vals = rel_df.loc[rel_df['detected'], 'rel_abundance_pct']
    floor = positive_vals.min() / 2 if len(positive_vals) else 0.01
    global_max = positive_vals.max() if len(positive_vals) else floor * 10
    n_contrasts = len(contrasts)
    bracket_headroom = global_max * (1.15 + 0.35 * n_contrasts)  # top bracket + label clearance
    if log_scale:
        ylim_shared = (floor * 0.6, bracket_headroom * 1.15)
    else:
        pad = global_max * 0.1
        ylim_shared = (0, bracket_headroom + pad)

    rng = np.random.default_rng(0)  # fixed seed -- jitter is reproducible across reruns

    for r, bg in enumerate(backgrounds):
        for c, tp in enumerate(timepoints):
            ax = axes[r][c]
            cell = rel_df[(rel_df[background_col] == bg) & (rel_df[time_col] == tp)]

            box_data, positions, present_treatments = [], [], []
            for i, treat in enumerate(treatment_order):
                vals = cell.loc[cell[treatment_col] == treat, 'rel_abundance_pct']
                vals_for_box = vals.where(vals > 0, floor)  # zeros -> floor, for box stats only
                if len(vals_for_box) == 0:
                    continue
                box_data.append(vals_for_box.values)
                positions.append(i + 1)
                present_treatments.append(treat)

            if not box_data:
                ax.axis('off')
                continue

            bp = ax.boxplot(box_data, positions=positions, widths=0.55, showfliers=False,
                             patch_artist=True)
            for patch in bp['boxes']:
                patch.set_facecolor(TAXON_COLOR)
                patch.set_alpha(0.25)
                patch.set_edgecolor(TAXON_COLOR)
            for element in ('whiskers', 'caps', 'medians'):
                for line in bp[element]:
                    line.set_color(TAXON_COLOR)

            # jittered individual points, hollow marker for non-detected (floored) ones
            for i, treat in zip(positions, present_treatments):
                sub = cell[cell[treatment_col] == treat]
                x_jit = i + rng.uniform(-0.12, 0.12, size=len(sub))
                detected = sub['detected'].values
                y = np.where(detected, sub['rel_abundance_pct'].values, floor)
                ax.scatter(x_jit[detected], y[detected], color=TAXON_COLOR, s=14, alpha=0.75,
                           edgecolor='none', zorder=3)
                ax.scatter(x_jit[~detected], y[~detected], facecolor='white', edgecolor=TAXON_COLOR,
                           s=14, alpha=0.9, linewidth=0.8, zorder=3)

            if log_scale:
                ax.set_yscale('log')
            ax.set_ylim(ylim_shared)

            # significance brackets: nested, contrasts[0] innermost (lowest).
            # Positioned off the SHARED global_max (not this panel's own local
            # max) so bracket height is consistent across every panel --
            # comparable at a glance, and never silently drifts because one
            # panel happens to have a lower ceiling than its neighbors.
            step = (global_max * 0.35) if log_scale else (global_max * 0.12 + 1e-9)
            for k, (treat_a, treat_b) in enumerate(contrasts):
                if treat_a not in treatment_order or treat_b not in treatment_order:
                    continue
                erow = effects_df[
                    (effects_df['background'] == bg) & (effects_df['timepoint'] == tp) &
                    (effects_df['taxon'] == taxon) & (effects_df['contrast'] == f'{treat_a} - {treat_b}')
                ]
                q = erow['q_value'].iloc[0] if len(erow) else np.nan
                label = significance_stars(q)
                xa, xb = treatment_order.index(treat_a) + 1, treatment_order.index(treat_b) + 1
                y_bracket = global_max * (1.15 + 0.35 * k) if log_scale else global_max + step * (k + 1)
                draw_significance_bracket(ax, min(xa, xb), max(xa, xb), y_bracket,
                                           height=y_bracket * 0.06 if log_scale else step * 0.3,
                                           text=label)

            ax.set_xticks(range(1, len(treatment_order) + 1))
            ax.set_xticklabels(treatment_order, rotation=30, ha='right', fontsize=7)
            ax.tick_params(axis='y', labelsize=7, labelleft=(c == 0))
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            if r == 0:
                tp_label = f'{tp:g}' if isinstance(tp, (int, float, np.floating)) else str(tp)
                ax.set_title(tp_label, fontsize=9)
            if c == n_cols - 1:
                ax.text(1.05, 0.5, bg, transform=ax.transAxes, fontsize=9,
                        rotation=-90, va='center', ha='left')
            if c == 0:
                ax.set_ylabel(f'{taxon} relative abundance (%)', fontsize=8)

    fig.suptitle(title, fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f'  wrote {out_path}')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--community', required=True)
    p.add_argument('--in-vitro-metadata', required=True)
    p.add_argument('--in-planta-metadata', required=True)
    p.add_argument('--invitro-effects-csv', required=True,
                    help="invitro_tailocin_effects.csv from 08_tailocin_effect_plot.py")
    p.add_argument('--inplanta-effects-csv', required=True,
                    help="inplanta_trial2_tailocin_effects.csv from 08_tailocin_effect_plot.py")
    p.add_argument('--output-dir', required=True)
    p.add_argument('--taxon', default='p25C2')
    p.add_argument('--invitro-backgrounds', default='S1+C2',
                    help="Comma-separated community_background values to plot for in vitro. "
                         "Default: 'S1+C2' only -- S1 has no deliberately-added p25C2 (same "
                         "role as SynCom1 in planta: a background contamination-level read, "
                         "not a real inoculated member), so it's dropped by default rather than "
                         "shown as a redundant near-null comparison. Pass 'S1,S1+C2' to restore it.")
    p.add_argument('--inplanta-backgrounds', default='SynCom1C',
                    help="Comma-separated community_background values to plot for in planta. "
                         "Default: 'SynCom1C' only, same reasoning as --invitro-backgrounds -- "
                         "SynCom1 has no deliberately-added p25C2. Pass 'SynCom1,SynCom1C' to restore it.")
    p.add_argument('--linear-scale', action='store_true',
                    help="Use a linear y-axis instead of the default log scale.")
    args = p.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    log_scale = not args.linear_scale

    comm = pd.read_csv(args.community)

    # ---------------- in vitro ----------------
    print('in vitro')
    iv = pd.read_csv(args.in_vitro_metadata, dtype={'STRAIN': str})
    exp_iv = iv[iv['sample_type'] == 'experimental']
    exp_iv = exp_iv[exp_iv['sample_uid'].isin(comm['sample_uid'])]

    rel_iv = compute_relative_abundance(comm, exp_iv['sample_uid'].tolist(), args.taxon)
    rel_iv = rel_iv.merge(
        exp_iv[['sample_uid', 'STRAIN', 'TREATMENT', 'timepoint']].rename(
            columns={'STRAIN': 'community_background', 'TREATMENT': 'treatment'}),
        on='sample_uid', how='left')
    rel_iv.to_csv(os.path.join(args.output_dir, f'invitro_{args.taxon}_relative_abundance.csv'), index=False)

    effects_iv = pd.read_csv(args.invitro_effects_csv)
    invitro_backgrounds = [b.strip() for b in args.invitro_backgrounds.split(',')]
    plot_taxon_boxplot(
        rel_iv, effects_iv, args.taxon,
        background_col='community_background', time_col='timepoint', treatment_col='treatment',
        backgrounds=invitro_backgrounds, treatment_order=['Buffer', 'delTFA', 'Tailocin'],
        contrasts=[('Tailocin', 'delTFA'), ('Tailocin', 'Buffer')],
        title=f'In vitro: {args.taxon} relative abundance by treatment',
        out_path=os.path.join(args.output_dir, f'invitro_{args.taxon}_boxplot.pdf'),
        log_scale=log_scale,
    )

    # ---------------- in planta trial 2 ----------------
    print('in planta trial 2')
    ip = pd.read_csv(args.in_planta_metadata)
    ip['trial'] = pd.to_numeric(ip['trial'], errors='coerce')
    exp_ip = ip[(ip['sample_type'] == 'experimental') & (ip['trial'] == 2.0)]
    exp_ip = exp_ip[exp_ip['truly_inoculated'] == True]  # noqa: E712
    exp_ip = exp_ip[exp_ip['sample_uid'].isin(comm['sample_uid'])]

    rel_ip = compute_relative_abundance(comm, exp_ip['sample_uid'].tolist(), args.taxon)
    rel_ip = rel_ip.merge(
        exp_ip[['sample_uid', 'community', 'treatment', 'dpi']].rename(
            columns={'community': 'community_background'}),
        on='sample_uid', how='left')
    rel_ip.to_csv(os.path.join(args.output_dir, f'inplanta_trial2_{args.taxon}_relative_abundance.csv'), index=False)

    effects_ip = pd.read_csv(args.inplanta_effects_csv)
    inplanta_backgrounds = [b.strip() for b in args.inplanta_backgrounds.split(',')]
    plot_taxon_boxplot(
        rel_ip, effects_ip, args.taxon,
        background_col='community_background', time_col='dpi', treatment_col='treatment',
        backgrounds=inplanta_backgrounds, treatment_order=['Buf', 'TFA', 'Tail'],
        contrasts=[('Tail', 'TFA'), ('Tail', 'Buf')],
        title=f'In planta trial 2: {args.taxon} relative abundance by treatment',
        out_path=os.path.join(args.output_dir, f'inplanta_trial2_{args.taxon}_boxplot.pdf'),
        log_scale=log_scale,
    )

    print(f'\nDone. Outputs in {args.output_dir}/')


if __name__ == '__main__':
    main()