"""
06_tailocin_effect_plot.py

Per-member tailocin effect plot from CLR-transformed data (output of
05a_clr_transform.py): for each community background x timepoint/dpi group,
estimate the CLR difference between Tailocin and each comparison treatment
(Buffer, and the delTFA/TFA structural-mutant control) for every retained
taxon, with a 95% CI and BH-FDR-corrected significance.

This does NOT use PERMANOVA (that's the community-level test in
07a/07b) -- this is a per-member univariate comparison, one Welch's t-test
per (background, timepoint, taxon, contrast), matching the "which specific
members respond to tailocin" question the reference figure is answering.

FDR correction (Benjamini-Hochberg) is applied ONCE across every test shown
in a given output plot (all backgrounds x timepoints x taxa x contrasts
together) -- i.e. the family of comparisons is "everything in this figure",
matching the single FDR<0.05 threshold used in the legend. If you later
pare this down to one timepoint, rerun rather than reusing an old q-value --
the correction is a function of how many tests are in the family.

Currently plots in vitro and in planta trial 2 only (per request) -- trial 1
CLR output exists (from 05a_clr_transform.py) but is not plotted here since
it wasn't asked for; add a call to plot_effects() with trial 1's CLR file to
include it later.

NOTE ON CLR BASIS: this script does not itself decide how CLR was computed --
it just reads whatever's in the clr_value column of the input file. For the
S1+C2 / SynCom1C backgrounds, 05a_clr_transform.py's dual-CLR handling means
p25C2's clr_value comes from a CLR pass that includes p25C2 as a part, while
every other taxon's clr_value comes from a separate pass with p25C2 excluded
from parts. Both live in the same clr_value column, distinguished only by
the clr_basis column carried through into the results CSV -- so a p25C2 row
and a Leaf21 row in the same background/timepoint are NOT on the same CLR
scale as each other.

Produces, per system, a main plot (Tailocin vs Buffer only -- the comparison
that actually matters for the figure), a supplemental plot (adds the
delTFA/TFA structural-mutant contrast alongside it), and a collapsed plot
(Tailocin vs Buffer only, all timepoints/dpi pooled into one overall estimate
per member per background -- for a quick "does tailocin affect this member
at all, averaged across the whole time course" view, separate from the
timepoint-resolved main plot). All three are built from the CLR data
directly; the collapsed view runs its own t-tests and its own FDR
correction (a different, smaller test family than the timepoint-resolved
plots) rather than summarizing their per-timepoint results. Use
--skip-supplemental-plot to skip the delTFA/TFA supplemental plot only --
the collapsed plot is always produced.

Usage:
    python 06_tailocin_effect_plot.py \
        --clr-invitro clr_out/clr_invitro_long.csv \
        --clr-inplanta-trial2 clr_out/clr_inplanta_trial2_long.csv \
        --output-dir effect_plots \
        [--skip-supplemental-plot]
"""

import argparse
import os
import re

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt


def numeric_sort_key(val):
    """Sort timepoints/dpi correctly regardless of type -- '6H' < '12H' < '18H'
    < '24H' (not alphabetical), and plain numeric dpi values sort as numbers."""
    s = str(val)
    m = re.match(r'^(\d+(?:\.\d+)?)', s)
    return float(m.group(1)) if m else float('inf')

# Fixed taxon display order (top to bottom in each panel), matching the
# existing reference figure where possible; any taxon present in the data
# but not listed here is appended at the end (so a newly-retained taxon
# never silently disappears from the plot).
TAXON_ORDER_PREFERRED = ['Leaf21', 'Leaf234', 'Leaf59', 'Leaf233', 'Leaf179',
                          'Leaf13', 'Leaf132', 'Leaf201', 'Leaf84', 'p25C2']

# --- Styling --------------------------------------------------------------
# Color encodes the CONTRAST. Each contrast gets a paired dark/light set from
# NatParksPalettes "DeathValley": the dark color for confident taxa, the
# light color (at reduced opacity) for low-confidence taxa (low_confidence_flag
# from 05). Fill (solid vs hollow) still encodes FDR significance; marker
# shape is kept as a second cue for contrast so the plots stay readable in
# grayscale / for color-vision deficiency.
CONTRAST_STYLES = {
    0: {'dark': '#60485B', 'light': '#B47E83', 'marker': 'o'},  # Tailocin - Buffer
    1: {'dark': '#132F5B', 'light': '#435F90', 'marker': '^'},  # Tailocin - delTFA/TFA
}
LOW_CONFIDENCE_ALPHA = 0.55


def point_style(contrast_idx, flagged=False):
    """(color, alpha, marker) for one point: dark/opaque for confident taxa,
    light/faded for low-confidence taxa, within the contrast's color pair."""
    st = CONTRAST_STYLES[contrast_idx % len(CONTRAST_STYLES)]
    if flagged:
        return st['light'], LOW_CONFIDENCE_ALPHA, st['marker']
    return st['dark'], 1.0, st['marker']


def legend_handles(contrast_labels):
    """Legend: per contrast, a confident (dark) and low-confidence (light)
    entry; then solid/hollow for FDR, drawn in neutral gray."""
    handles = []
    for i, lbl in enumerate(contrast_labels):
        for flagged, suffix in [(False, ''), (True, ' (low-confidence taxon)')]:
            color, alpha, marker = point_style(i, flagged)
            handles.append(plt.Line2D([0], [0], color=color, marker=marker, markerfacecolor=color,
                                      markeredgecolor=color, markersize=8, linestyle='-',
                                      alpha=alpha, label=lbl + suffix))
    handles.append(plt.Line2D([0], [0], color='dimgray', marker='o', markerfacecolor='dimgray',
                              markersize=8, linestyle='none', label='FDR < 0.05 (solid)'))
    handles.append(plt.Line2D([0], [0], color='dimgray', marker='o', markerfacecolor='white',
                              markersize=8, linestyle='none', label='n.s. (hollow)'))
    return handles

# Font sizes -- edit here to rescale every label in every figure at once.
FONT = {
    'tick': 11,         # x/y tick labels (taxon names, CLR values)
    'col_title': 13,    # timepoint/dpi column headers, background titles
    'row_label': 12,    # background labels on the right edge
    'axis_label': 13,   # shared x-axis label
    'legend': 11,
    'suptitle': 15,
}


# p25C2 exclusion rule: p25C2 was only ever actually inoculated into the
# "+C2"/"C" variant of each system's community. It must be excluded from
# every OTHER background's results even when amplicon reads assigned to
# p25C2 show up there anyway (e.g. index-hopping, or edit-distance overlap
# with Leaf59 -- see project notes on the two being nearly identical in the
# amplicon window). This is a hard exclusion, not the flag-based fading used
# for contamination-floor taxa elsewhere in this script: presence of reads
# doesn't mean presence in the design, and no fade/significance call should
# ever be computed for it in these backgrounds in the first place.
P25C2_INOCULATED_BACKGROUNDS = {'S1+C2', 'SynCom1C'}


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR correction. Returns array of q-values, same
    order as input. Implemented directly (no statsmodels dependency)."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]  # enforce monotonicity
    q_final = np.empty(n)
    q_final[order] = np.clip(q, 0, 1)
    return q_final


def welch_ttest_ci(a, b, alpha=0.05):
    """Welch's t-test, mean(a) - mean(b), with (1-alpha) CI. Returns
    (estimate, ci_low, ci_high, p_value) or all-NaN if either group has <2
    observations (can't estimate variance)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan, np.nan, np.nan

    mean_diff = a.mean() - b.mean()
    var_a, var_b = a.var(ddof=1), b.var(ddof=1)
    na, nb = len(a), len(b)
    se = np.sqrt(var_a / na + var_b / nb)
    if se == 0:
        return mean_diff, mean_diff, mean_diff, 0.0 if mean_diff != 0 else 1.0

    # Welch-Satterthwaite degrees of freedom
    df = (var_a / na + var_b / nb) ** 2 / (
        (var_a / na) ** 2 / (na - 1) + (var_b / nb) ** 2 / (nb - 1)
    )
    t_crit = stats.t.ppf(1 - alpha / 2, df)
    ci_low = mean_diff - t_crit * se
    ci_high = mean_diff + t_crit * se
    t_stat = mean_diff / se
    p_value = 2 * stats.t.sf(abs(t_stat), df)
    return mean_diff, ci_low, ci_high, p_value


def compute_effects(clr_long, background_col, time_col, treatment_col,
                     backgrounds, contrasts):
    """
    clr_long: sample_uid, taxon, clr_value, plus background_col/time_col/treatment_col
    contrasts: list of (treatment_a, treatment_b) tuples; estimate = mean(a) - mean(b)
    Returns a results df: background, timepoint, taxon, contrast_label, estimate,
    ci_low, ci_high, p_value, q_value, significant
    """
    rows = []
    taxa = sorted(clr_long['taxon'].unique())
    timepoints = sorted(clr_long[time_col].dropna().unique(), key=numeric_sort_key)

    for bg in backgrounds:
        for tp in timepoints:
            grp = clr_long[(clr_long[background_col] == bg) & (clr_long[time_col] == tp)]
            taxa_this_bg = [t for t in taxa if not (t == 'p25C2' and bg not in P25C2_INOCULATED_BACKGROUNDS)]
            for taxon in taxa_this_bg:
                tsub = grp[grp['taxon'] == taxon]
                clr_basis = tsub['clr_basis'].iloc[0] if 'clr_basis' in tsub.columns and len(tsub) > 0 else 'standard'
                flagged = bool(tsub['low_confidence_flag'].iloc[0]) if 'low_confidence_flag' in tsub.columns and len(tsub) > 0 else False
                for treat_a, treat_b in contrasts:
                    a = tsub[tsub[treatment_col] == treat_a]['clr_value']
                    b = tsub[tsub[treatment_col] == treat_b]['clr_value']
                    est, lo, hi, p = welch_ttest_ci(a, b)
                    rows.append({
                        'background': bg, 'timepoint': tp, 'taxon': taxon,
                        'contrast': f'{treat_a} - {treat_b}',
                        'clr_basis': clr_basis,
                        'flagged': flagged,
                        'n_a': len(a), 'n_b': len(b),
                        'estimate': est, 'ci_low': lo, 'ci_high': hi, 'p_value': p,
                    })

    results = pd.DataFrame(rows)
    valid = results['p_value'].notna()
    results['q_value'] = np.nan
    results.loc[valid, 'q_value'] = bh_fdr(results.loc[valid, 'p_value'].values)
    results['significant'] = results['q_value'] < 0.05
    return results



def plot_effects(results, contrasts_to_plot, title, out_path, taxon_order=None):
    """contrasts_to_plot: list of (treatment_a, treatment_b) tuples to actually
    draw -- can be a subset of what's in `results` (e.g. just Tailocin-Buffer
    for a main figure, even though `results` also has Tailocin-delTFA computed
    for supplemental use). Dodge is only applied when >1 contrast is drawn."""
    contrast_labels = [f'{a} - {b}' for a, b in contrasts_to_plot]
    results = results[results['contrast'].isin(contrast_labels)]

    backgrounds = list(dict.fromkeys(results['background']))  # preserve order
    timepoints = sorted(results['timepoint'].unique(), key=numeric_sort_key)

    all_taxa = sorted(results.loc[results['estimate'].notna(), 'taxon'].unique())
    if taxon_order is None:
        ordered = [t for t in TAXON_ORDER_PREFERRED if t in all_taxa]
        ordered += [t for t in all_taxa if t not in ordered]
        taxon_order = ordered
    y_pos = {t: i for i, t in enumerate(taxon_order)}
    n_taxa = len(taxon_order)

    n_rows, n_cols = len(backgrounds), len(timepoints)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.2 * n_cols + 0.8, 0.55 * n_taxa * n_rows + 1.6),
                              squeeze=False)

    finite_vals = pd.concat([results['ci_low'], results['ci_high'], results['estimate']]).dropna()
    if len(finite_vals) > 0:
        xpad = 0.08 * (finite_vals.max() - finite_vals.min() + 1e-9)
        xlim = (finite_vals.min() - xpad, finite_vals.max() + xpad)
    else:
        xlim = (-1, 1)

    dodge = 0.15 if len(contrast_labels) > 1 else 0.0
    for r, bg in enumerate(backgrounds):
        for c, tp in enumerate(timepoints):
            ax = axes[r][c]
            cell = results[(results['background'] == bg) & (results['timepoint'] == tp)]

            for ci, clabel in enumerate(contrast_labels):
                sub = cell[cell['contrast'] == clabel]
                offset = (ci - (len(contrast_labels) - 1) / 2) * dodge
                for _, row in sub.iterrows():
                    if row['taxon'] not in y_pos or pd.isna(row['estimate']):
                        continue
                    y = y_pos[row['taxon']] + offset
                    color, point_alpha, marker = point_style(ci, row.get('flagged', False))
                    ax.plot([row['ci_low'], row['ci_high']], [y, y],
                            color=color, linewidth=1.8, zorder=2, alpha=point_alpha)
                    facecolor = color if row['significant'] else 'white'
                    ax.scatter([row['estimate']], [y], marker=marker, facecolor=facecolor,
                               edgecolor=color, s=55, zorder=3, linewidth=1.4, alpha=point_alpha)

            ax.axvline(0, color='gray', linewidth=0.9, linestyle='--', zorder=1)
            ax.set_yticks(range(n_taxa))
            ax.set_yticklabels(taxon_order, fontsize=FONT['tick'])
            ax.tick_params(axis='y', labelleft=(c == 0))
            ax.tick_params(axis='x', labelsize=FONT['tick'])
            ax.set_xlim(xlim)
            ax.set_ylim(-0.5, n_taxa - 0.5)
            ax.invert_yaxis()
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            if r == 0:
                tp_label = f'{tp:g}' if isinstance(tp, (int, float, np.floating)) else str(tp)
                ax.set_title(tp_label, fontsize=FONT['col_title'])
            if c == n_cols - 1:
                ax.text(1.03, 0.5, bg, transform=ax.transAxes, fontsize=FONT['row_label'],
                        rotation=-90, va='center', ha='left')

    # Label/legend/title offsets are set in INCHES and converted to figure
    # fractions, so the spacing stays the same however tall the figure is
    fig_h = fig.get_size_inches()[1]
    fig.tight_layout(rect=[0, 0, 1, 1])
    fig.text(0.5, -0.35 / fig_h, 'Estimated CLR effect (contrast ± 95% CI)', ha='center',
             va='top', fontsize=FONT['axis_label'])

    handles = legend_handles(contrast_labels)
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.75 / fig_h),
               ncol=len(handles) // 2, fontsize=FONT['legend'], frameon=False)
    fig.suptitle(title, fontsize=FONT['suptitle'], y=1 + 0.25 / fig_h, va='bottom')
    fig.savefig(out_path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f'  wrote {out_path}')


def compute_effects_collapsed(clr_long, background_col, treatment_col, backgrounds, contrast,
                               replicate_id_col=None):
    """
    Same idea as compute_effects(), but pools every timepoint/dpi together
    within each (background, taxon) group before running the t-test -- gives
    a single "what's the overall effect of tailocin on this member, averaged
    across the whole time course" estimate per background, instead of one
    estimate per timepoint.

    replicate_id_col: REQUIRED for any system with repeated measures on the
    same physical unit across timepoints (in vitro: pass 'BIOREP_ID'). Since
    each biorep_id splits into one subculture per treatment, grouping by
    BIOREP_ID within a single treatment arm is exactly grouping by
    subculture -- so this averages each subculture's repeated timepoint
    measurements into ONE value before the t-test, giving n = number of
    independent subcultures, not n = number of subculture-timepoint rows.
    Skipping this (leaving it None) for a repeated-measures system silently
    pseudoreplicates: it treats correlated repeated measurements from the
    same physical unit as if they were independent replicates, which
    understates the true variance and makes the p-value anti-conservative
    (falsely confident), not just underpowered. This is the same category of
    problem 09a's split-plot restricted permutation exists to prevent at the
    community level -- this function needs the analogous care at the
    per-taxon level.
    Leave as None for destructive-sampling systems (in planta) where each row
    already IS an independent physical unit and no aggregation is needed.

    Only ONE contrast (e.g. Tailocin vs Buffer) since that's the whole point
    of the collapsed view here -- pass a single (treatment_a, treatment_b)
    tuple, not a list.

    Returns a results df: background, taxon, contrast, estimate, ci_low,
    ci_high, p_value, q_value, significant, n_a, n_b (n_a/n_b are counts of
    independent replicates actually tested, i.e. subcultures when
    replicate_id_col is given, not raw rows) -- FDR corrected across this
    collapsed family only (backgrounds x taxa), separate from any
    timepoint-resolved family computed elsewhere.
    """
    treat_a, treat_b = contrast
    rows = []
    taxa = sorted(clr_long['taxon'].unique())

    for bg in backgrounds:
        grp = clr_long[clr_long[background_col] == bg]
        taxa_this_bg = [t for t in taxa if not (t == 'p25C2' and bg not in P25C2_INOCULATED_BACKGROUNDS)]
        for taxon in taxa_this_bg:
            tsub = grp[grp['taxon'] == taxon]
            clr_basis = tsub['clr_basis'].iloc[0] if 'clr_basis' in tsub.columns and len(tsub) > 0 else 'standard'
            flagged = bool(tsub['low_confidence_flag'].iloc[0]) if 'low_confidence_flag' in tsub.columns and len(tsub) > 0 else False

            def get_group_values(treat):
                sub = tsub[tsub[treatment_col] == treat]
                if replicate_id_col and replicate_id_col in sub.columns:
                    # collapse each independent unit's repeated timepoint
                    # rows to one mean CLR value BEFORE testing
                    return sub.groupby(replicate_id_col)['clr_value'].mean()
                return sub['clr_value']

            a = get_group_values(treat_a)
            b = get_group_values(treat_b)
            est, lo, hi, p = welch_ttest_ci(a, b)
            rows.append({
                'background': bg, 'taxon': taxon,
                'contrast': f'{treat_a} - {treat_b}',
                'clr_basis': clr_basis,
                'flagged': flagged,
                'n_a': len(a), 'n_b': len(b),
                'estimate': est, 'ci_low': lo, 'ci_high': hi, 'p_value': p,
            })

    results = pd.DataFrame(rows)
    valid = results['p_value'].notna()
    results['q_value'] = np.nan
    results.loc[valid, 'q_value'] = bh_fdr(results.loc[valid, 'p_value'].values)
    results['significant'] = results['q_value'] < 0.05
    return results



def plot_effects_collapsed(results, title, out_path, taxon_order=None):
    """One row of panels (no timepoint facet), one column per background,
    single contrast, no dodge needed since there's only one series per row."""
    backgrounds = list(dict.fromkeys(results['background']))  # preserve order

    all_taxa = sorted(results.loc[results['estimate'].notna(), 'taxon'].unique())
    if taxon_order is None:
        ordered = [t for t in TAXON_ORDER_PREFERRED if t in all_taxa]
        ordered += [t for t in all_taxa if t not in ordered]
        taxon_order = ordered
    y_pos = {t: i for i, t in enumerate(taxon_order)}
    n_taxa = len(taxon_order)

    n_cols = len(backgrounds)
    fig, axes = plt.subplots(1, n_cols, figsize=(3.4 * n_cols + 0.8, 0.55 * n_taxa + 1.8), squeeze=False)
    axes = axes[0]

    finite_vals = pd.concat([results['ci_low'], results['ci_high'], results['estimate']]).dropna()
    if len(finite_vals) > 0:
        xpad = 0.08 * (finite_vals.max() - finite_vals.min() + 1e-9)
        xlim = (finite_vals.min() - xpad, finite_vals.max() + xpad)
    else:
        xlim = (-1, 1)

    for c, bg in enumerate(backgrounds):
        ax = axes[c]
        cell = results[results['background'] == bg]
        for _, row in cell.iterrows():
            if row['taxon'] not in y_pos or pd.isna(row['estimate']):
                continue
            y = y_pos[row['taxon']]
            color, point_alpha, marker = point_style(0, row.get('flagged', False))
            ax.plot([row['ci_low'], row['ci_high']], [y, y], color=color, linewidth=2.0,
                    zorder=2, alpha=point_alpha)
            facecolor = color if row['significant'] else 'white'
            ax.scatter([row['estimate']], [y], marker=marker, facecolor=facecolor,
                       edgecolor=color, s=65, zorder=3, linewidth=1.5, alpha=point_alpha)

        ax.axvline(0, color='gray', linewidth=0.9, linestyle='--', zorder=1)
        ax.set_yticks(range(n_taxa))
        ax.set_yticklabels(taxon_order, fontsize=FONT['tick'])
        ax.tick_params(axis='y', labelleft=(c == 0))
        ax.tick_params(axis='x', labelsize=FONT['tick'])
        ax.set_xlim(xlim)
        ax.set_ylim(-0.5, n_taxa - 0.5)
        ax.invert_yaxis()
        ax.set_title(bg, fontsize=FONT['col_title'])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    contrast_label = results['contrast'].iloc[0] if len(results) else ''
    fig_h = fig.get_size_inches()[1]
    fig.tight_layout(rect=[0, 0, 1, 1])
    fig.text(0.5, -0.35 / fig_h, f'Estimated CLR effect, {contrast_label}, all timepoints pooled (± 95% CI)',
              ha='center', va='top', fontsize=FONT['axis_label'])
    handles = legend_handles([contrast_label])
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.75 / fig_h),
               ncol=len(handles) // 2, fontsize=FONT['legend'], frameon=False)
    fig.suptitle(title, fontsize=FONT['suptitle'], y=1 + 0.25 / fig_h, va='bottom')
    fig.savefig(out_path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f'  wrote {out_path}')


def find_raw_timepoint(clr_long, time_col, target_numeric):
    """Find the actual raw value in clr_long[time_col] whose leading numeric
    portion equals target_numeric (e.g. target_numeric=24 matches a raw value
    of '24H', 24, '24.0', etc, regardless of exactly how it's stored) -- so
    callers don't need to know the exact string/dtype used upstream, just the
    number. Raises if nothing matches, rather than silently filtering to an
    empty frame."""
    candidates = [v for v in clr_long[time_col].dropna().unique()
                  if numeric_sort_key(v) == target_numeric]
    if not candidates:
        available = sorted(clr_long[time_col].dropna().unique(), key=numeric_sort_key)
        raise ValueError(f"No {time_col} value matching {target_numeric} found. "
                          f"Available values: {available}")
    return candidates[0]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--clr-invitro', required=True)
    p.add_argument('--clr-inplanta-trial2', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--skip-supplemental-plot', action='store_true',
                    help="Only produce the main (Tailocin-Buffer only) plot; skip the "
                         "supplemental plot with delTFA/TFA also drawn. The results CSV "
                         "always includes both contrasts either way.")
    args = p.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # --- in vitro ---
    print('in vitro')
    clr_iv = pd.read_csv(args.clr_invitro)
    contrasts_iv = [('Tailocin', 'Buffer'), ('Tailocin', 'delTFA')]
    results_iv = compute_effects(
        clr_iv, background_col='community_background', time_col='timepoint',
        treatment_col='treatment', backgrounds=['S1', 'S1+C2'], contrasts=contrasts_iv,
    )
    results_iv.to_csv(os.path.join(args.output_dir, 'invitro_tailocin_effects.csv'), index=False)
    plot_effects(results_iv, [contrasts_iv[0]],
                 'In vitro: estimated tailocin effect per member (Tailocin vs Buffer)',
                 os.path.join(args.output_dir, 'invitro_tailocin_effect_plot_main.pdf'))
    if not args.skip_supplemental_plot:
        plot_effects(results_iv, contrasts_iv,
                     'In vitro: estimated tailocin effect per member (supplemental: + delTFA)',
                     os.path.join(args.output_dir, 'invitro_tailocin_effect_plot_supplemental.pdf'))

    results_iv_collapsed = compute_effects_collapsed(
        clr_iv, background_col='community_background', treatment_col='treatment',
        backgrounds=['S1', 'S1+C2'], contrast=contrasts_iv[0], replicate_id_col='BIOREP_ID',
    )
    results_iv_collapsed.to_csv(os.path.join(args.output_dir, 'invitro_tailocin_effects_collapsed.csv'), index=False)
    plot_effects_collapsed(results_iv_collapsed,
                            'In vitro: overall tailocin effect per member (all timepoints pooled)',
                            os.path.join(args.output_dir, 'invitro_tailocin_effect_plot_collapsed.pdf'))

    # --- in planta trial 2 ---
    print('in planta trial 2')
    clr_ip = pd.read_csv(args.clr_inplanta_trial2)
    contrasts_ip = [('Tail', 'Buf'), ('Tail', 'TFA')]
    # Buffer community is excluded from the background rows here: there's no
    # inoculated member community to estimate a per-member tailocin effect
    # for in wells that were never inoculated with the SynCom.
    results_ip = compute_effects(
        clr_ip, background_col='community_background', time_col='dpi',
        treatment_col='treatment', backgrounds=['SynCom1', 'SynCom1C'], contrasts=contrasts_ip,
    )
    results_ip.to_csv(os.path.join(args.output_dir, 'inplanta_trial2_tailocin_effects.csv'), index=False)
    plot_effects(results_ip, [contrasts_ip[0]],
                 'In planta trial 2: estimated tailocin effect per member (Tailocin vs Buffer)',
                 os.path.join(args.output_dir, 'inplanta_trial2_tailocin_effect_plot_main.pdf'))
    if not args.skip_supplemental_plot:
        plot_effects(results_ip, contrasts_ip,
                     'In planta trial 2: estimated tailocin effect per member (supplemental: + TFA)',
                     os.path.join(args.output_dir, 'inplanta_trial2_tailocin_effect_plot_supplemental.pdf'))

    # NOTE: no replicate_id_col here (unlike in vitro) -- in planta is
    # destructive sampling, each row is already an independent seedling, no
    # repeated-measures aggregation needed before pooling across dpi.
    results_ip_collapsed = compute_effects_collapsed(
        clr_ip, background_col='community_background', treatment_col='treatment',
        backgrounds=['SynCom1', 'SynCom1C'], contrast=contrasts_ip[0],
    )
    results_ip_collapsed.to_csv(os.path.join(args.output_dir, 'inplanta_trial2_tailocin_effects_collapsed.csv'), index=False)
    plot_effects_collapsed(results_ip_collapsed,
                            'In planta trial 2: overall tailocin effect per member (all dpi pooled)',
                            os.path.join(args.output_dir, 'inplanta_trial2_tailocin_effect_plot_collapsed.pdf'))

    print('standalone single-panel plots (main figures)')

    tp_24h = find_raw_timepoint(results_iv, 'timepoint', 24)
    for bg, fname_bg in [('S1', 'S1'), ('S1+C2', 'S1C2')]:
        sub = results_iv[(results_iv['background'] == bg) & (results_iv['timepoint'] == tp_24h)]
        sub.to_csv(os.path.join(args.output_dir, f'invitro_{fname_bg}_24H_tailocin_effects.csv'), index=False)
        plot_effects(sub, [contrasts_iv[0]],
                     f'In vitro, 24H: estimated tailocin effect per member (Tailocin vs Buffer) -- {bg}',
                     os.path.join(args.output_dir, f'invitro_{fname_bg}_24H_tailocin_effect_plot.pdf'))

    dpi3 = find_raw_timepoint(results_ip, 'timepoint', 3)
    for bg in ['SynCom1', 'SynCom1C']:
        sub = results_ip[(results_ip['background'] == bg) & (results_ip['timepoint'] == dpi3)]
        sub.to_csv(os.path.join(args.output_dir, f'inplanta_trial2_{bg}_3dpi_tailocin_effects.csv'), index=False)
        plot_effects(sub, [contrasts_ip[0]],
                     f'In planta trial 2, 3 dpi: estimated tailocin effect per member (Tailocin vs Buffer) -- {bg}',
                     os.path.join(args.output_dir, f'inplanta_trial2_{bg}_3dpi_tailocin_effect_plot.pdf'))

    print(f'\nDone. Outputs in {args.output_dir}/')


if __name__ == '__main__':
    main()
