"""
08_tailocin_effect_plot.py

Per-member tailocin effect plot from CLR-transformed data (output of
07_clr_transform.py): for each community background x timepoint/dpi group,
estimate the CLR difference between Tailocin and each comparison treatment
(Buffer, and the delTFA/TFA structural-mutant control) for every retained
taxon, with a 95% CI and BH-FDR-corrected significance.

This deliberately does NOT use PERMANOVA (that's the community-level test in
05a/05b) -- this is a per-member univariate comparison, one Welch's t-test
per (background, timepoint, taxon, contrast), matching the "which specific
members respond to tailocin" question the reference figure is answering.

FDR correction (Benjamini-Hochberg) is applied ONCE across every test shown
in a given output plot (all backgrounds x timepoints x taxa x contrasts
together) -- i.e. the family of comparisons is "everything in this figure",
matching the single FDR<0.05 threshold used in the legend. If you later
pare this down to one timepoint, rerun rather than reusing an old q-value --
the correction is a function of how many tests are in the family.

Currently plots in vitro and in planta trial 2 only (per request) -- trial 1
CLR output exists (from 07_clr_transform.py) but is not plotted here since
it wasn't asked for; add a call to plot_effects() with trial 1's CLR file to
include it later.

NOTE ON CLR BASIS: this script does not itself decide how CLR was computed --
it just reads whatever's in the clr_value column of the input file. For the
S1+C2 / SynCom1C backgrounds, 07_clr_transform.py's dual-CLR handling means
p25C2's clr_value comes from a CLR pass that includes p25C2 as a part, while
every other taxon's clr_value comes from a separate pass with p25C2 excluded
from parts. Both live in the same clr_value column, distinguished only by
the clr_basis column carried through into the results CSV -- so a p25C2 row
and a Leaf21 row in the same background/timepoint are NOT on the same CLR
scale as each other. That's intentional (see 07's docstring), but worth
remembering if you ever want to compare magnitudes across taxa within a
panel rather than just within a taxon's own row across timepoints/contrasts.

Usage:
    python 08_tailocin_effect_plot.py \
        --clr-invitro clr_out/clr_invitro_long.csv \
        --clr-inplanta-trial2 clr_out/clr_inplanta_trial2_long.csv \
        --output-dir effect_plots
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

CONTRAST_COLORS = {0: '#1f77b4', 1: '#d62728'}  # contrast index -> color


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
            for taxon in taxa:
                tsub = grp[grp['taxon'] == taxon]
                clr_basis = tsub['clr_basis'].iloc[0] if 'clr_basis' in tsub.columns and len(tsub) > 0 else 'standard'
                for treat_a, treat_b in contrasts:
                    a = tsub[tsub[treatment_col] == treat_a]['clr_value']
                    b = tsub[tsub[treatment_col] == treat_b]['clr_value']
                    est, lo, hi, p = welch_ttest_ci(a, b)
                    rows.append({
                        'background': bg, 'timepoint': tp, 'taxon': taxon,
                        'contrast': f'{treat_a} - {treat_b}',
                        'clr_basis': clr_basis,
                        'n_a': len(a), 'n_b': len(b),
                        'estimate': est, 'ci_low': lo, 'ci_high': hi, 'p_value': p,
                    })

    results = pd.DataFrame(rows)
    valid = results['p_value'].notna()
    results['q_value'] = np.nan
    results.loc[valid, 'q_value'] = bh_fdr(results.loc[valid, 'p_value'].values)
    results['significant'] = results['q_value'] < 0.05
    return results


def plot_effects(results, contrasts, title, out_path, taxon_order=None):
    contrast_labels = [f'{a} - {b}' for a, b in contrasts]
    backgrounds = list(dict.fromkeys(results['background']))  # preserve order
    timepoints = sorted(results['timepoint'].unique(), key=numeric_sort_key)

    all_taxa = sorted(results['taxon'].unique())
    if taxon_order is None:
        ordered = [t for t in TAXON_ORDER_PREFERRED if t in all_taxa]
        ordered += [t for t in all_taxa if t not in ordered]
        taxon_order = ordered
    y_pos = {t: i for i, t in enumerate(taxon_order)}
    n_taxa = len(taxon_order)

    # NOTE: sharey=True is deliberately NOT used here. With shared y-axes,
    # matplotlib shares the tick-label objects across all subplots in the
    # row/column group, so the LAST subplot's set_yticklabels() call silently
    # overwrites the labels on every other subplot sharing that axis --
    # this cost real debugging time the first pass (every panel came out
    # blank). Instead, ylim is matched manually below and only column-0
    # panels get labelleft=True.
    n_rows, n_cols = len(backgrounds), len(timepoints)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(2.6 * n_cols, 0.45 * n_taxa * n_rows + 1),
                              squeeze=False)

    # shared x-range across every panel (computed manually since sharex is
    # avoided for the same tick-label-sharing reason as sharey above)
    finite_vals = pd.concat([results['ci_low'], results['ci_high'], results['estimate']]).dropna()
    if len(finite_vals) > 0:
        xpad = 0.08 * (finite_vals.max() - finite_vals.min() + 1e-9)
        xlim = (finite_vals.min() - xpad, finite_vals.max() + xpad)
    else:
        xlim = (-1, 1)

    dodge = 0.15
    for r, bg in enumerate(backgrounds):
        for c, tp in enumerate(timepoints):
            ax = axes[r][c]
            cell = results[(results['background'] == bg) & (results['timepoint'] == tp)]

            for ci, clabel in enumerate(contrast_labels):
                sub = cell[cell['contrast'] == clabel]
                color = CONTRAST_COLORS[ci % len(CONTRAST_COLORS)]
                offset = (ci - (len(contrast_labels) - 1) / 2) * dodge
                for _, row in sub.iterrows():
                    if row['taxon'] not in y_pos or pd.isna(row['estimate']):
                        continue
                    y = y_pos[row['taxon']] + offset
                    ax.plot([row['ci_low'], row['ci_high']], [y, y],
                            color=color, linewidth=1.2, zorder=2)
                    facecolor = color if row['significant'] else 'white'
                    ax.scatter([row['estimate']], [y], color=color, facecolor=facecolor,
                               edgecolor=color, s=22, zorder=3, linewidth=1.0)

            ax.axvline(0, color='gray', linewidth=0.8, linestyle='--', zorder=1)
            ax.set_yticks(range(n_taxa))
            ax.set_yticklabels(taxon_order, fontsize=7)
            ax.tick_params(axis='y', labelleft=(c == 0))
            ax.tick_params(axis='x', labelsize=7)
            ax.set_xlim(xlim)
            ax.set_ylim(-0.5, n_taxa - 0.5)
            ax.invert_yaxis()
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            if r == 0:
                tp_label = f'{tp:g}' if isinstance(tp, (int, float, np.floating)) else str(tp)
                ax.set_title(tp_label, fontsize=9)
            if c == n_cols - 1:
                ax.text(1.02, 0.5, bg, transform=ax.transAxes, fontsize=8,
                        rotation=-90, va='center', ha='left')

    fig.text(0.5, -0.06, 'Estimated CLR effect (contrast ± 95% CI)', ha='center', fontsize=10)

    handles = [plt.Line2D([0], [0], color=CONTRAST_COLORS[i], marker='o',
                           markerfacecolor=CONTRAST_COLORS[i], linestyle='-', label=lbl)
               for i, lbl in enumerate(contrast_labels)]
    handles.append(plt.Line2D([0], [0], color='gray', marker='o', markerfacecolor='gray',
                               linestyle='none', label='FDR < 0.05 (solid)'))
    handles.append(plt.Line2D([0], [0], color='gray', marker='o', markerfacecolor='white',
                               linestyle='none', label='n.s. (hollow)'))
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.11),
               ncol=2, fontsize=8, frameon=False)
    fig.suptitle(title, fontsize=12, y=1.02)
    fig.tight_layout(rect=[0, 0.08, 1, 0.97])
    fig.savefig(out_path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f'  wrote {out_path}')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--clr-invitro', required=True)
    p.add_argument('--clr-inplanta-trial2', required=True)
    p.add_argument('--output-dir', required=True)
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
    plot_effects(results_iv, contrasts_iv,
                 'In vitro: estimated tailocin effect per member',
                 os.path.join(args.output_dir, 'invitro_tailocin_effect_plot.pdf'))

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
    plot_effects(results_ip, contrasts_ip,
                 'In planta trial 2: estimated tailocin effect per member',
                 os.path.join(args.output_dir, 'inplanta_trial2_tailocin_effect_plot.pdf'))

    print(f'\nDone. Outputs in {args.output_dir}/')


if __name__ == '__main__':
    main()
