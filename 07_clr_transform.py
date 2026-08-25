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
- Two filters gate whether a taxon is kept as a CLR part, both computed per
  trial (not globally across trials):
    1. Prevalence: non-zero in >= --prevalence-threshold fraction of samples
       (presence/absence only).
    2. Magnitude: mean relative abundance (raw count / that sample's TOTAL
       reads, including unassigned -- not just the parts total) >=
       --min-mean-relative-abundance, averaged only over the samples where
       the taxon is actually present (zeros excluded from that average, so a
       taxon isn't penalized for being absent in some samples, only for
       being negligible in the samples where it IS detected).
    Both must pass. The prevalence-only filter was letting through taxa that
    are technically "present" in most samples but only ever at a handful of
    stray reads (single-digit-to-low-double-digit counts against
    100,000+-read libraries) -- the magnitude filter is what actually catches
    that, prevalence alone does not. A taxon dropped by either filter is
    reported separately in the printed log so it's clear which criterion
    caught it.
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


def build_clr(comm, sample_uids, exclude_from_parts=None, prevalence_threshold=0.1,
              min_mean_relative_abundance=0.0):
    """
    comm: long-format df with sample_uid, strain_or_group (raw label), R2_count
    sample_uids: samples to include (already filtered by caller: system,
        trial, sample_type, truly_inoculated, etc.)
    min_mean_relative_abundance: a taxon must ALSO average at least this
        fraction of a sample's TOTAL reads (across the samples where it's
        actually detected, i.e. zeros excluded from the average) to be kept
        as a part -- catches taxa that are "present" in most samples but only
        ever at a handful of stray reads. 0.0 disables this filter (prevalence
        alone still applies).
    Returns:
        clr_long: sample_uid, taxon, clr_value
        coverage: sample_uid, total_reads, part_reads, part_fraction
                  (part_fraction low = sample is mostly unassigned/ambiguous/
                  excluded -- worth checking before trusting its CLR values)
        kept_taxa, dropped_taxa: for the printed report
    """
    exclude_from_parts = set(exclude_from_parts or [])
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

    n_samples = len(wide)
    prevalence = (wide > 0).sum(axis=0) / max(n_samples, 1)
    passes_prevalence = prevalence >= prevalence_threshold

    # relative abundance denominator = each sample's TOTAL reads (matches the
    # definition used when this was diagnosed by hand -- count / everything
    # in the well, including unassigned -- not count / parts-only total,
    # which would inflate the percentage for samples that are mostly
    # unassigned).
    total_reads_aligned = coverage.set_index('sample_uid')['total_reads'].reindex(wide.index)
    rel_abund = wide.div(total_reads_aligned.replace(0, np.nan), axis=0)
    # mean over samples where the taxon is actually present, zeros excluded
    # from the average (a taxon isn't penalized for legitimate absence
    # elsewhere, only for being negligible where it IS detected)
    present_mask = wide > 0
    mean_rel_when_present = rel_abund.where(present_mask).mean(axis=0).fillna(0)
    passes_magnitude = mean_rel_when_present >= min_mean_relative_abundance

    kept_taxa = prevalence[passes_prevalence & passes_magnitude].index.tolist()
    dropped_by_prevalence = prevalence[~passes_prevalence].index.tolist()
    dropped_by_magnitude = prevalence[passes_prevalence & ~passes_magnitude].index.tolist()
    if dropped_by_prevalence:
        print(f'    dropped by prevalence filter (<{prevalence_threshold:.0%} of samples): '
              f'{dropped_by_prevalence}')
    if dropped_by_magnitude:
        detail = ', '.join(f'{t} ({mean_rel_when_present[t]*100:.3f}% mean when present)'
                            for t in dropped_by_magnitude)
        print(f'    dropped by magnitude filter (<{min_mean_relative_abundance*100:.3f}% mean '
              f'relative abundance when present): {detail}')
    dropped_taxa = dropped_by_prevalence + dropped_by_magnitude
    wide = wide[kept_taxa]

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
    return clr_long, coverage, kept_taxa, dropped_taxa


def run_dual_clr(comm, sample_uids, meta_cols, dual_backgrounds, dual_taxon,
                  exclude_from_parts, prevalence_threshold, min_mean_relative_abundance,
                  background_col):
    """
    Runs two CLR passes over the same sample set:
      - standard: dual_taxon included as a part
      - excluded: dual_taxon removed from parts before geomean
    and merges them so that, for samples whose background is in
    dual_backgrounds: dual_taxon's row comes from the standard pass, every
    other taxon's row comes from the excluded pass. Samples with backgrounds
    NOT in dual_backgrounds get the standard pass in full, unchanged.

    Returns the merged long-format df (with clr_basis column), plus the two
    raw per-pass long-format dfs (for optional diagnostic output) and the
    coverage diagnostic from the standard pass.
    """
    standard_long, coverage, kept, dropped = build_clr(
        comm, sample_uids, exclude_from_parts, prevalence_threshold, min_mean_relative_abundance)
    excluded_long, _, kept_excl, dropped_excl = build_clr(
        comm, sample_uids, list(exclude_from_parts) + [dual_taxon], prevalence_threshold,
        min_mean_relative_abundance)

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
                  min_mean_relative_abundance, dual_backgrounds, dual_taxon):
    print('  in vitro')
    exp = iv[iv['sample_type'] == 'experimental'].copy()
    exp = exp[exp['sample_uid'].isin(comm['sample_uid'])]

    meta_cols = exp[['sample_uid', 'STRAIN', 'TREATMENT', 'timepoint', 'BIOREP_ID', 'well',
                      'sequencing_plate']].rename(
        columns={'STRAIN': 'community_background', 'TREATMENT': 'treatment'})

    out, standard_full, excluded_full, coverage = run_dual_clr(
        comm, exp['sample_uid'].tolist(), meta_cols, dual_backgrounds, dual_taxon,
        exclude_from_parts, prevalence_threshold, min_mean_relative_abundance,
        background_col='community_background')

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
                   min_mean_relative_abundance, dual_backgrounds, dual_taxon):
    print(f'  in planta trial {int(trial)}')
    exp = ip[(ip['sample_type'] == 'experimental') & (ip['trial'] == trial)].copy()
    exp = exp[exp['truly_inoculated'] == True]  # noqa: E712
    exp = exp[exp['sample_uid'].isin(comm['sample_uid'])]

    meta_cols = exp[['sample_uid', 'community', 'treatment', 'dpi', 'collection_well',
                      'growth_plate_raw', 'collection_plate', 'trial']].rename(
        columns={'community': 'community_background', 'growth_plate_raw': 'growth_plate',
                 'collection_plate': 'seq_plate'})

    out, standard_full, excluded_full, coverage = run_dual_clr(
        comm, exp['sample_uid'].tolist(), meta_cols, dual_backgrounds, dual_taxon,
        exclude_from_parts, prevalence_threshold, min_mean_relative_abundance,
        background_col='community_background')

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
    p.add_argument('--min-mean-relative-abundance', type=float, default=0.0005,
                    help="Minimum mean relative abundance (fraction of a sample's TOTAL reads, "
                         "e.g. 0.0005 = 0.05%%), averaged only over samples where the taxon is "
                         "actually detected, required to keep a taxon as a CLR part -- catches "
                         "taxa that are technically 'present' in most samples but only ever at a "
                         "handful of stray reads. Default 0.0005 (0.05%%). Set to 0 to disable "
                         "(prevalence-only filtering, the old behavior).")
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
                 args.min_mean_relative_abundance, args.invitro_dual_clr_backgrounds, args.dual_clr_taxon)
    run_in_planta(comm, ip, 1.0, args.output_dir, args.exclude_from_parts, args.prevalence_threshold,
                  args.min_mean_relative_abundance, args.inplanta_dual_clr_backgrounds, args.dual_clr_taxon)
    run_in_planta(comm, ip, 2.0, args.output_dir, args.exclude_from_parts, args.prevalence_threshold,
                  args.min_mean_relative_abundance, args.inplanta_dual_clr_backgrounds, args.dual_clr_taxon)

    print(f'\nDone. Outputs in {args.output_dir}/')


if __name__ == '__main__':
    main()
