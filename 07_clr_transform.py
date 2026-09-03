"""
07_clr_transform.py

CLR (centered log-ratio) transformation of strain-level composition data,
run separately for each trial (in vitro; in planta trial 1; in planta trial
2), since they have different experimental designs and different sample
sets. One invocation produces all three outputs.

Design:
- "Unassigned" reads and the ambiguous "59+p25C2" call are NOT included as
  CLR parts -- they don't represent a real, identifiable community member,
  so putting them in the composition would treat noise/ambiguity as if it
  were a taxon with its own abundance. They ARE reported in a coverage
  diagnostic (what fraction of each sample's reads are usable as real parts)
  so a sample that's mostly unassigned doesn't silently get treated as
  equivalent to one that's mostly resolved.
- FLAG, DON'T DROP (as of this version): earlier versions of this script
  excluded taxa from the CLR parts entirely if they failed a
  prevalence/magnitude check. That's now understood to be the wrong move for
  a KNOWN, deliberately-inoculated synthetic community: you don't need
  statistical evidence a taxon is "real" when you know for a fact it was
  added to the inoculum -- that's a design fact, not a hypothesis to test.
  Dropping a real member also corrupts the CLR geometry for every OTHER
  taxon too, since the geometric mean is computed over whatever parts remain.
  So every taxon that appears at all (prevalence > 0) is now KEPT as a CLR
  part, always -- what changed is that two boolean flag columns
  (low_prevalence_flag, below_noise_floor_flag, combined into
  low_confidence_flag) are carried through into the output CSV instead of
  silently removing the taxon. Downstream scripts and the manuscript can
  caveat or exclude on this basis themselves; the CLR computation itself
  always uses the full, real compositional structure.
  EXCEPTION: the dual-CLR taxon (p25C2) is NEVER flagged low_confidence
  (everything is already always kept now, so there's no separate "keep"
  exemption to make -- just the flag exemption). Its abundance is legitimately, expectedly driven
  down by the Tailocin treatment itself -- that's the entire point of the
  experiment, not a sign of an unreliable measurement. Comparing its
  Tailocin-depleted abundance against a contamination floor and flagging it
  "low confidence" would be actively misleading: low abundance there is the
  real biological effect, not noise. See run_dual_clr for exactly where this
  exemption is applied.
- Two criteria feed the flag (a taxon is flagged if EITHER fails):
    1. Prevalence: non-zero in >= --prevalence-threshold fraction of samples.
    2. Magnitude: MEDIAN relative abundance (raw count / that sample's TOTAL
       reads, including unassigned), among samples where present, must
       clear EITHER --blank-noise-floor-multiplier times that taxon's own
       median relative abundance in extraction_blank wells for this trial,
       OR exceed the blank level by at least --absolute-margin-pct
       percentage points. Both a ratio and an absolute-margin criterion are
       needed together -- confirmed empirically that a ratio-only version
       incorrectly penalizes abundant, obviously-real taxa (blank
       contamination scales with a taxon's true overall abundance, so a
       dominant taxon can never clear "Nx its own blank level" once it
       already makes up a large share of the community), while a flat
       magnitude-only version let a mean-skewed low-abundance taxon (Leaf13)
       through undetected. MEDIAN (not mean) is used on both sides of this
       comparison for the same reason it was chosen everywhere else in this
       pipeline: mean is too sensitive to a handful of high-outlier wells.
    See compute_blank_noise_floor() for why blank-derived and per-taxon
    (not a flat percentage borrowed from natural-community literature):
    extraction blanks in this dataset are NOT clean (median ~12,000-15,000
    reads each, composition that mirrors the real community roughly
    proportional to each taxon's overall abundance -- consistent with index
    hopping / cross-well contamination during library prep). A flat
    threshold can't account for contamination scaling with abundance; each
    taxon needs its own empirically-measured floor.
- Zero handling: zeros are replaced with a pseudocount equal to half the
  smallest non-zero raw count observed anywhere in that trial's part matrix
  (a standard, dependency-free multiplicative replacement -- avoids pulling
  in zCompositions/statsmodels just for this).
- --exclude-dropped removes truly_inoculated == False wells (in planta) or
  any wells otherwise flagged, before the prevalence filter and geometric
  mean are computed, so a known-bad well can't skew either.
- --exclude-from-parts TAXON can be used to drop a specific taxon (e.g.
  p25C2) from the CLR computation entirely, for analyses that want CLR
  computed on the remaining member basis (kept separate from the prevalence
  filter, since this is an intentional exclusion, not a data-driven one).
- DUAL-CLR BACKGROUNDS: for community backgrounds where p25C2 was
  deliberately added as an inoculated member (in vitro S1+C2, in planta
  SynCom1C), two CLR passes are run and merged into the final output:
    1. a "standard" pass with p25C2 included as a part -- only p25C2's own
       row is kept from this pass, since that's the version that correctly
       reflects p25C2's own change in the community it's actually part of.
    2. a "p25C2_excluded" pass with p25C2 removed from parts before the
       geometric mean is computed -- every OTHER taxon's row is kept from
       this pass, so their CLR values reflect the community with p25C2's own
       change factored out, rather than being pulled around by it.
  Every output row carries a `clr_basis` column recording which pass it came
  from ('standard', 'standard_p25C2_only', or 'p25C2_excluded') so this
  choice is traceable downstream rather than silently baked in. Backgrounds
  not listed as dual-CLR just get the single standard pass, unchanged.
  Configurable via --invitro-dual-clr-backgrounds / --inplanta-dual-clr-backgrounds
  (defaults: S1+C2 and SynCom1C) and --dual-clr-taxon (default p25C2).

Output: one long-format CSV per trial (sample_uid, taxon, clr_value, plus
relevant metadata columns for downstream modeling/plotting), plus a
part-coverage diagnostic CSV per trial.

Usage:
    python 07_clr_transform.py \
        --community community_strain_abundance.csv \
        --in-vitro-metadata in_vitro_master_metadata.csv \
        --in-planta-metadata in_planta_master_metadata.csv \
        --output-dir clr_output \
        [--exclude-from-parts p25C2] \
        [--prevalence-threshold 0.1]
"""

import argparse
import os
import re

import numpy as np
import pandas as pd

NON_MEMBER_LABELS = {'unassigned', '59+p25c2'}  # lowercased for matching


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


def compute_blank_noise_floor(comm, blank_sample_uids):
    """
    Per-taxon MEDIAN relative abundance (of that sample's total reads, incl.
    unassigned -- same denominator convention used everywhere else) across
    extraction_blank wells. This is the empirical contamination floor a real
    taxon's signal needs to clear -- see build_clr()'s docstring for why
    median (not mean) and why blank-derived (not a flat percentage).
    Returns a Series indexed by relabeled taxon name. Empty Series if no
    blank wells / no blank read data is available for this trial.
    """
    if not blank_sample_uids:
        return pd.Series(dtype=float)
    sub = comm[comm['sample_uid'].isin(blank_sample_uids)].copy()
    if len(sub) == 0:
        return pd.Series(dtype=float)
    sub['taxon'] = sub['strain_or_group'].astype(str).map(relabel)
    total = sub.groupby('sample_uid')['R2_count'].transform('sum')
    sub = sub.assign(rel_pct=(sub['R2_count'] / total.replace(0, np.nan)) * 100)
    return sub.groupby('taxon')['rel_pct'].median()


def build_clr(comm, sample_uids, exclude_from_parts=None, prevalence_threshold=0.1,
              blank_noise_floor=None, noise_floor_multiplier=3.0,
              fallback_min_median_relative_abundance=0.0, absolute_margin_pct=1.0,
              never_flag_taxa=None):
    """
    comm: long-format df with sample_uid, strain_or_group (raw label), R2_count
    sample_uids: samples to include (already filtered by caller: system,
        trial, sample_type, truly_inoculated, etc.)
    blank_noise_floor: Series (taxon -> median blank relative abundance %)
        from compute_blank_noise_floor(), or None/empty to use the flat
        fallback threshold instead.
    noise_floor_multiplier / absolute_margin_pct: see module docstring --
        a taxon is flagged low_confidence if it fails BOTH the ratio test
        (median >= multiplier x blank floor) AND the absolute-margin test
        (median - blank floor >= absolute_margin_pct percentage points).
    fallback_min_median_relative_abundance: flat threshold (fraction, e.g.
        0.0005 = 0.05%) used ONLY for taxa missing from blank_noise_floor.
    never_flag_taxa: taxa (e.g. the dual-CLR taxon, p25C2) that are NEVER
        flagged low_confidence regardless of what the prevalence/magnitude
        checks say -- see module docstring for why (their low abundance
        under treatment is the real effect being measured, not a
        data-quality concern).
    Returns:
        clr_long: sample_uid, taxon, clr_value, low_prevalence_flag,
                  below_noise_floor_flag, low_confidence_flag
        coverage: sample_uid, total_reads, part_reads, part_fraction
        kept_taxa: ALL taxa with any nonzero reads (nothing is dropped by
                   prevalence/magnitude anymore -- only truly never-observed
                   taxa are excluded, since there's no data to include)
        flag_detail: DataFrame, one row per taxon, with the underlying
                     numbers behind each flag.
    """
    exclude_from_parts = set(exclude_from_parts or [])
    never_flag_taxa = set(never_flag_taxa or [])
    blank_noise_floor = blank_noise_floor if blank_noise_floor is not None else pd.Series(dtype=float)
    sub = comm[comm['sample_uid'].isin(sample_uids)].copy()

    sub['is_member'] = ~sub['strain_or_group'].astype(str).str.lower().isin(NON_MEMBER_LABELS)
    sub['taxon'] = sub['strain_or_group'].astype(str).map(relabel)

    total_reads = sub.groupby('sample_uid')['R2_count'].sum().rename('total_reads')
    part_reads = sub[sub['is_member']].groupby('sample_uid')['R2_count'].sum().rename('part_reads')
    coverage = pd.concat([total_reads, part_reads], axis=1).fillna(0)
    coverage['part_fraction'] = (coverage['part_reads'] / coverage['total_reads']).fillna(0)
    coverage = coverage.reindex(sample_uids, fill_value=0).reset_index().rename(columns={'index': 'sample_uid'})

    parts = sub[sub['is_member']]
    wide = parts.pivot_table(index='sample_uid', columns='taxon', values='R2_count',
                              aggfunc='sum', fill_value=0)
    wide = wide.reindex(sample_uids, fill_value=0)

    if exclude_from_parts:
        drop_cols = [c for c in wide.columns if c in exclude_from_parts or
                     c.replace('Leaf', '') in exclude_from_parts]
        if drop_cols:
            wide = wide.drop(columns=drop_cols)
            print(f'    excluded from parts (--exclude-from-parts): {drop_cols}')

    # The only thing that still actually removes a column: a taxon with
    # literally zero reads anywhere in this sample set -- no compositional
    # information to include, a data-availability fact, not a QC judgment.
    never_observed = wide.columns[(wide == 0).all(axis=0)].tolist()
    if never_observed:
        print(f'    never observed in this sample set (excluded, not flagged -- no data to include): {never_observed}')
        wide = wide.drop(columns=never_observed)

    n_samples = len(wide)
    prevalence = (wide > 0).sum(axis=0) / max(n_samples, 1)
    low_prevalence_flag = prevalence < prevalence_threshold

    total_reads_aligned = coverage.set_index('sample_uid')['total_reads'].reindex(wide.index)
    rel_abund = wide.div(total_reads_aligned.replace(0, np.nan), axis=0) * 100
    present_mask = wide > 0
    median_rel_when_present = rel_abund.where(present_mask).median(axis=0).fillna(0)

    required_threshold = {}
    threshold_basis = {}
    for taxon in wide.columns:
        if taxon in blank_noise_floor.index and pd.notna(blank_noise_floor[taxon]):
            required_threshold[taxon] = blank_noise_floor[taxon] * noise_floor_multiplier
            threshold_basis[taxon] = f'{noise_floor_multiplier:g}x blank floor ({blank_noise_floor[taxon]:.4f}%)'
        else:
            required_threshold[taxon] = fallback_min_median_relative_abundance * 100
            threshold_basis[taxon] = f'flat fallback ({fallback_min_median_relative_abundance*100:.3f}%, no blank data for this taxon)'
    required_threshold = pd.Series(required_threshold)
    passes_ratio = median_rel_when_present >= required_threshold
    blank_floor_aligned = pd.Series({t: blank_noise_floor.get(t, 0.0) for t in wide.columns})
    passes_absolute_margin = (median_rel_when_present - blank_floor_aligned) >= absolute_margin_pct
    below_noise_floor_flag = ~(passes_ratio | passes_absolute_margin)

    low_confidence_flag = low_prevalence_flag | below_noise_floor_flag

    # never_flag_taxa exemption: force flags False regardless of the
    # computed checks -- for the dual-CLR taxon (p25C2), whose expected
    # treatment-driven depletion should never read as a data-quality flag.
    for t in never_flag_taxa:
        if t in low_confidence_flag.index:
            low_prevalence_flag[t] = False
            below_noise_floor_flag[t] = False
            low_confidence_flag[t] = False

    flagged = low_confidence_flag[low_confidence_flag].index.tolist()
    if flagged:
        for t in flagged:
            reasons = []
            if low_prevalence_flag.get(t, False):
                reasons.append(f'prevalence {prevalence[t]:.1%} < {prevalence_threshold:.0%}')
            if below_noise_floor_flag.get(t, False):
                reasons.append(f'median {median_rel_when_present[t]:.4f}% when present vs blank floor '
                                f'{blank_floor_aligned[t]:.4f}% [{threshold_basis.get(t, "n/a")}]')
            print(f'    FLAGGED low_confidence: {t} -- {"; ".join(reasons)}')
    if never_flag_taxa:
        exempted = sorted(t for t in never_flag_taxa if t in wide.columns)
        if exempted:
            print(f'    never-flagged (exempt from low_confidence regardless of checks): {exempted}')

    flag_detail = pd.DataFrame({
        'prevalence': prevalence, 'median_when_present_pct': median_rel_when_present,
        'blank_floor_pct': blank_floor_aligned, 'required_threshold_pct': required_threshold,
        'low_prevalence_flag': low_prevalence_flag, 'below_noise_floor_flag': below_noise_floor_flag,
        'low_confidence_flag': low_confidence_flag,
    })

    kept_taxa = wide.columns.tolist()  # nothing dropped by prevalence/magnitude anymore

    nonzero_vals = wide.values[wide.values > 0]
    if len(nonzero_vals) == 0:
        raise RuntimeError('No non-zero counts remain after filtering -- check inputs.')
    pseudocount = nonzero_vals.min() / 2.0
    replaced = wide.replace(0, pseudocount)

    log_vals = np.log(replaced.values)
    geo_mean_log = log_vals.mean(axis=1, keepdims=True)
    clr_vals = log_vals - geo_mean_log
    clr_wide = pd.DataFrame(clr_vals, index=wide.index, columns=wide.columns)

    clr_long = clr_wide.reset_index().melt(id_vars='sample_uid', var_name='taxon', value_name='clr_value')
    clr_long = clr_long.merge(
        flag_detail[['low_prevalence_flag', 'below_noise_floor_flag', 'low_confidence_flag']]
        .reset_index().rename(columns={'index': 'taxon'}),
        on='taxon', how='left')
    return clr_long, coverage, kept_taxa, flag_detail


def taxon_passes_noise_floor(comm, sample_uids, taxon, prevalence_threshold,
                              blank_noise_floor, noise_floor_multiplier,
                              fallback_min_median_relative_abundance, absolute_margin_pct):
    """
    Same prevalence + (ratio OR absolute-margin) logic as build_clr's magnitude
    filter, but for a single named taxon over a specific, caller-chosen
    sample_uids subset -- used to give the dual-CLR taxon a background-scoped
    check instead of the trial-wide one. See run_dual_clr for why: a taxon
    deliberately present in only ONE background (e.g. p25C2 in S1+C2, not S1)
    gets its real signal diluted below the noise floor if judged against the
    WHOLE trial's pooled samples, since its near-absent contamination-level
    presence in the other background(s) drags the pooled statistic down.
    Confirmed empirically: without this, p25C2 was dropped from the pipeline
    entirely.
    Returns (passes: bool, median_when_present: float, blank_floor: float).
    """
    sub = comm[comm['sample_uid'].isin(sample_uids)].copy()
    sub['taxon'] = sub['strain_or_group'].astype(str).map(relabel)
    total = sub.groupby('sample_uid')['R2_count'].transform('sum')
    sub = sub.assign(rel_pct=(sub['R2_count'] / total.replace(0, np.nan)) * 100)

    tsub = sub[sub['taxon'] == taxon]
    n_samples = len(sample_uids)
    prevalence = tsub['sample_uid'].nunique() / max(n_samples, 1)
    if prevalence < prevalence_threshold:
        return False, 0.0, blank_noise_floor.get(taxon, 0.0)

    present = tsub[tsub['R2_count'] > 0]
    median_when_present = present['rel_pct'].median() if len(present) else 0.0
    floor = blank_noise_floor.get(taxon, np.nan)
    if pd.notna(floor):
        required = floor * noise_floor_multiplier
        passes = (median_when_present >= required) or (median_when_present - floor >= absolute_margin_pct)
    else:
        passes = median_when_present >= fallback_min_median_relative_abundance * 100
    return passes, median_when_present, (floor if pd.notna(floor) else 0.0)


def run_dual_clr(comm, sample_uids, meta_cols, dual_backgrounds, dual_taxon,
                  exclude_from_parts, prevalence_threshold, blank_noise_floor,
                  noise_floor_multiplier, fallback_min_median_relative_abundance,
                  absolute_margin_pct, background_col, treatment_col='treatment',
                  baseline_treatment=None):
    """
    Runs two CLR passes over the same sample set:
      - standard: dual_taxon included as a part
      - excluded: dual_taxon removed from parts before geomean
    and merges them so that, for samples whose background is in
    dual_backgrounds: dual_taxon's row comes from the standard pass, every
    other taxon's row comes from the excluded pass. Samples with backgrounds
    NOT in dual_backgrounds get the standard pass in full, unchanged.

    dual_taxon is ALWAYS exempt from the low_confidence flag (passed via
    never_flag_taxa into build_clr) -- its abundance is legitimately,
    expectedly driven down by the very treatment under study, so a
    contamination-floor comparison would misleadingly flag the real effect
    as a data-quality problem. This is unconditional, not a "check first,
    then decide" gate -- there's nothing to decide, the exemption always
    applies. The diagnostic check below (scoped to dual_backgrounds and, if
    given, restricted to baseline_treatment) is printed purely for
    transparency, so the underlying numbers are visible in the log even
    though they no longer drive any decision. Both scoping restrictions
    were found necessary when this WAS still a drop/keep decision (an
    earlier version of this pipeline): pooling across backgrounds diluted
    p25C2's real S1+C2 signal below the noise floor entirely, and even
    pooling across treatments within S1+C2 dragged it down further via the
    Tailocin-driven depletion that's the whole point of the experiment.
    Kept here as a printed diagnostic since it's still informative context,
    just no longer a gate.

    Returns the merged long-format df (with clr_basis column), plus the two
    raw per-pass long-format dfs (for optional diagnostic output) and the
    coverage diagnostic from the standard pass.
    """
    never_flag = {dual_taxon} if dual_backgrounds else set()
    if dual_backgrounds:
        scope_mask = meta_cols[background_col].isin(dual_backgrounds)
        scope_label = str(dual_backgrounds)
        if baseline_treatment is not None and treatment_col in meta_cols.columns:
            scope_mask = scope_mask & (meta_cols[treatment_col] == baseline_treatment)
            scope_label = f'{dual_backgrounds}, {treatment_col}=={baseline_treatment!r} only'
        dual_bg_uids = meta_cols.loc[scope_mask, 'sample_uid'].tolist()
        passes, med, floor = taxon_passes_noise_floor(
            comm, dual_bg_uids, dual_taxon, prevalence_threshold, blank_noise_floor,
            noise_floor_multiplier, fallback_min_median_relative_abundance, absolute_margin_pct)
        verdict = 'would PASS anyway' if passes else 'would FAIL on its own -- exemption is doing real work here'
        print(f'    {dual_taxon} diagnostic (scoped to {scope_label}, n={len(dual_bg_uids)}, '
              f'informational only -- always exempt from low_confidence regardless): '
              f'median {med:.4f}% vs blank floor {floor:.4f}% -- {verdict}')

    standard_long, coverage, kept, flag_detail = build_clr(
        comm, sample_uids, exclude_from_parts, prevalence_threshold,
        blank_noise_floor, noise_floor_multiplier, fallback_min_median_relative_abundance,
        absolute_margin_pct, never_flag_taxa=never_flag)
    excluded_long, _, kept_excl, flag_detail_excl = build_clr(
        comm, sample_uids, list(exclude_from_parts) + [dual_taxon], prevalence_threshold,
        blank_noise_floor, noise_floor_multiplier, fallback_min_median_relative_abundance,
        absolute_margin_pct)

    standard_long = standard_long.merge(meta_cols, on='sample_uid', how='left')
    excluded_long = excluded_long.merge(meta_cols, on='sample_uid', how='left')

    is_dual = standard_long[background_col].isin(dual_backgrounds)
    keep_standard = standard_long[
        ((standard_long['taxon'] == dual_taxon) & is_dual) | (~is_dual)
    ].copy()
    keep_standard['clr_basis'] = np.where(
        keep_standard[background_col].isin(dual_backgrounds),
        f'standard_{dual_taxon}_only', 'standard')

    is_dual_excl = excluded_long[background_col].isin(dual_backgrounds)
    keep_excluded = excluded_long[
        (excluded_long['taxon'] != dual_taxon) & is_dual_excl
    ].copy()
    keep_excluded['clr_basis'] = f'{dual_taxon}_excluded'

    merged = pd.concat([keep_standard, keep_excluded], ignore_index=True)
    if dual_backgrounds:
        print(f'    dual-CLR applied for background(s) {dual_backgrounds}: '
              f'{dual_taxon} from standard pass, all other taxa from '
              f'{dual_taxon}-excluded pass')
    return merged, standard_long, excluded_long, coverage


def run_in_vitro(comm, iv, out_dir, exclude_from_parts, prevalence_threshold,
                  noise_floor_multiplier, fallback_min_median_relative_abundance,
                  absolute_margin_pct, dual_backgrounds, dual_taxon):
    print('  in vitro')
    exp = iv[iv['sample_type'] == 'experimental'].copy()
    exp = exp[exp['sample_uid'].isin(comm['sample_uid'])]

    blanks = iv[iv['sample_type'] == 'extraction_blank']
    blanks = blanks[blanks['sample_uid'].isin(comm['sample_uid'])]
    blank_noise_floor = compute_blank_noise_floor(comm, blanks['sample_uid'].tolist())
    if len(blank_noise_floor) == 0:
        print('    WARNING: no extraction_blank read data found for in vitro -- falling back to '
              f'flat {fallback_min_median_relative_abundance*100:.3f}% threshold for all taxa')
    else:
        print(f'    blank noise floor computed from {len(blanks)} extraction_blank wells '
              f'({noise_floor_multiplier:g}x median blank level required to pass)')

    meta_cols = exp[['sample_uid', 'STRAIN', 'TREATMENT', 'timepoint', 'BIOREP_ID', 'well',
                      'sequencing_plate']].rename(
        columns={'STRAIN': 'community_background', 'TREATMENT': 'treatment'})

    out, standard_full, excluded_full, coverage = run_dual_clr(
        comm, exp['sample_uid'].tolist(), meta_cols, dual_backgrounds, dual_taxon,
        exclude_from_parts, prevalence_threshold, blank_noise_floor, noise_floor_multiplier,
        fallback_min_median_relative_abundance, absolute_margin_pct, background_col='community_background',
        treatment_col='treatment', baseline_treatment='Buffer')

    out_path = os.path.join(out_dir, 'clr_invitro_long.csv')
    cov_path = os.path.join(out_dir, 'clr_invitro_part_coverage.csv')
    out.to_csv(out_path, index=False)
    coverage.to_csv(cov_path, index=False)
    print(f'    {out["taxon"].nunique()} taxa, {out["sample_uid"].nunique()} samples -> {out_path}')

    # Single-basis tables (NOT merged/dual) -- these are what PERMANOVA/Aitchison
    # distance work needs: every taxon in the matrix must be on the SAME CLR
    # basis, which the merged dual-basis table above deliberately is not (see
    # its docstring). standard = all taxa incl. dual_taxon as a part (Fig 3 /
    # direct-effect basis). {dual_taxon}_excluded = dual_taxon removed from
    # parts entirely for every sample, not just dual_backgrounds (Fig 4 /
    # indirect-effect basis, matches the old pipeline's --exclude_suffix run).
    standard_path = os.path.join(out_dir, 'clr_invitro_standard_long.csv')
    excluded_path = os.path.join(out_dir, f'clr_invitro_{dual_taxon}_excluded_long.csv')
    standard_full.to_csv(standard_path, index=False)
    excluded_full.to_csv(excluded_path, index=False)
    print(f'    single-basis tables -> {standard_path}, {excluded_path}')

    low_cov = coverage[coverage['part_fraction'] < 0.5]
    if len(low_cov) > 0:
        print(f'    NOTE: {len(low_cov)} samples are <50% real-member reads '
              f'(mostly unassigned/ambiguous) -- see {cov_path}')


def run_in_planta(comm, ip, trial, out_dir, exclude_from_parts, prevalence_threshold,
                   noise_floor_multiplier, fallback_min_median_relative_abundance,
                   absolute_margin_pct, dual_backgrounds, dual_taxon):
    print(f'  in planta trial {int(trial)}')
    exp = ip[(ip['sample_type'] == 'experimental') & (ip['trial'] == trial)].copy()
    exp = exp[exp['truly_inoculated'] == True]  # noqa: E712
    exp = exp[exp['sample_uid'].isin(comm['sample_uid'])]

    # NOTE: blank floor is computed from ALL of this trial's extraction
    # blanks, not filtered to `trial` on the metadata's own trial column --
    # blanks don't have a meaningful community/treatment/dpi of their own,
    # and pooling gives a more stable per-taxon median than splitting an
    # already-small blank set further. If blank contamination turns out to
    # differ meaningfully by trial, this would need revisiting.
    blanks = ip[ip['sample_type'] == 'extraction_blank']
    blanks = blanks[blanks['sample_uid'].isin(comm['sample_uid'])]
    blank_noise_floor = compute_blank_noise_floor(comm, blanks['sample_uid'].tolist())
    if len(blank_noise_floor) == 0:
        print('    WARNING: no extraction_blank read data found for in planta -- falling back to '
              f'flat {fallback_min_median_relative_abundance*100:.3f}% threshold for all taxa')
    else:
        print(f'    blank noise floor computed from {len(blanks)} extraction_blank wells '
              f'({noise_floor_multiplier:g}x median blank level required to pass)')

    meta_cols = exp[['sample_uid', 'community', 'treatment', 'dpi', 'collection_well',
                      'growth_plate_raw', 'collection_plate', 'trial']].rename(
        columns={'community': 'community_background', 'growth_plate_raw': 'growth_plate',
                 'collection_plate': 'seq_plate'})

    out, standard_full, excluded_full, coverage = run_dual_clr(
        comm, exp['sample_uid'].tolist(), meta_cols, dual_backgrounds, dual_taxon,
        exclude_from_parts, prevalence_threshold, blank_noise_floor, noise_floor_multiplier,
        fallback_min_median_relative_abundance, absolute_margin_pct, background_col='community_background',
        treatment_col='treatment', baseline_treatment='Buf')

    out_path = os.path.join(out_dir, f'clr_inplanta_trial{int(trial)}_long.csv')
    cov_path = os.path.join(out_dir, f'clr_inplanta_trial{int(trial)}_part_coverage.csv')
    out.to_csv(out_path, index=False)
    coverage.to_csv(cov_path, index=False)
    print(f'    {out["taxon"].nunique()} taxa, {out["sample_uid"].nunique()} samples -> {out_path}')

    standard_path = os.path.join(out_dir, f'clr_inplanta_trial{int(trial)}_standard_long.csv')
    excluded_path = os.path.join(out_dir, f'clr_inplanta_trial{int(trial)}_{dual_taxon}_excluded_long.csv')
    standard_full.to_csv(standard_path, index=False)
    excluded_full.to_csv(excluded_path, index=False)
    print(f'    single-basis tables -> {standard_path}, {excluded_path}')

    low_cov = coverage[coverage['part_fraction'] < 0.5]
    if len(low_cov) > 0:
        print(f'    NOTE: {len(low_cov)} samples are <50% real-member reads '
              f'(mostly unassigned/ambiguous) -- see {cov_path}')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--community', required=True)
    p.add_argument('--in-vitro-metadata', required=True)
    p.add_argument('--in-planta-metadata', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--exclude-from-parts', nargs='*', default=[],
                    help="Taxon labels to exclude from CLR parts entirely, in ADDITION to the "
                         "dual-CLR taxon handling below (rarely needed).")
    p.add_argument('--prevalence-threshold', type=float, default=0.1,
                    help="Minimum fraction of samples (within each trial) a taxon must be "
                         "non-zero in to be kept as a CLR part. Default 0.1.")
    p.add_argument('--blank-noise-floor-multiplier', type=float, default=3.0,
                    help="A taxon's median relative abundance (in real samples, among samples "
                         "where it's detected) must be >= this many times its own median relative "
                         "abundance in extraction_blank wells to be kept as a CLR part. Default 3.0.")
    p.add_argument('--min-median-relative-abundance', type=float, default=0.0005,
                    help="Flat fallback threshold (fraction, e.g. 0.0005 = 0.05%%), used ONLY for "
                         "taxa with no extraction_blank data available (or if no blanks exist for "
                         "a trial at all). Default 0.0005 (0.05%%).")
    p.add_argument('--absolute-margin-pct', type=float, default=1.0,
                    help="A taxon also passes the magnitude filter if its median real-sample "
                         "abundance exceeds its blank floor by at least this many PERCENTAGE "
                         "POINTS, regardless of the ratio test -- protects abundant, obviously-real "
                         "taxa from a ratio test that becomes impossible to clear once a taxon "
                         "already dominates the community (blank contamination scales with true "
                         "abundance too). Default 1.0 (one percentage point).")
    p.add_argument('--dual-clr-taxon', default='p25C2',
                    help="Taxon that gets the two-pass CLR treatment for the backgrounds "
                         "listed below. Default p25C2.")
    p.add_argument('--invitro-dual-clr-backgrounds', nargs='*', default=['S1+C2'],
                    help="in vitro STRAIN/community_background values that get the dual-CLR "
                         "treatment. Default: S1+C2. Pass nothing (empty list) to disable.")
    p.add_argument('--inplanta-dual-clr-backgrounds', nargs='*', default=['SynCom1C'],
                    help="in planta community values that get the dual-CLR treatment. "
                         "Default: SynCom1C. Pass nothing (empty list) to disable.")
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    comm = pd.read_csv(args.community)
    iv = pd.read_csv(args.in_vitro_metadata, dtype={'STRAIN': str})
    ip = pd.read_csv(args.in_planta_metadata)
    ip['trial'] = pd.to_numeric(ip['trial'], errors='coerce')

    print('Running CLR transform for all trials...')
    run_in_vitro(comm, iv, args.output_dir, args.exclude_from_parts, args.prevalence_threshold,
                 args.blank_noise_floor_multiplier, args.min_median_relative_abundance,
                 args.absolute_margin_pct, args.invitro_dual_clr_backgrounds, args.dual_clr_taxon)
    run_in_planta(comm, ip, 1.0, args.output_dir, args.exclude_from_parts, args.prevalence_threshold,
                  args.blank_noise_floor_multiplier, args.min_median_relative_abundance,
                  args.absolute_margin_pct, args.inplanta_dual_clr_backgrounds, args.dual_clr_taxon)
    run_in_planta(comm, ip, 2.0, args.output_dir, args.exclude_from_parts, args.prevalence_threshold,
                  args.blank_noise_floor_multiplier, args.min_median_relative_abundance,
                  args.absolute_margin_pct, args.inplanta_dual_clr_backgrounds, args.dual_clr_taxon)

    print(f'\nDone. Outputs in {args.output_dir}/')


if __name__ == '__main__':
    main()