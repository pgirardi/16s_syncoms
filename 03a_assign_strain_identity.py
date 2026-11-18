"""
Assign strain identity to community-sample ASVs using the monoculture
reference built by build_monoculture_reference_fasta.py.

Design:
- R1 and R2 were denoised INDEPENDENTLY -- an R1 ASV and R2 ASV from the same
  well are NOT linked read-by-read, only at the sample level.
- R2 is PRIMARY for strain identity: it's the direction with enough SNPs to
  resolve every strain, including Leaf59 vs p25C2 (identical/near-identical
  in the R1 window, 4 distinguishing SNPs in R2).
- R1 is used as an independent QC cross-check: for any strain group R1 can't 
  distinguish (an "ambiguous group"), we verify the R1 combined count for that 
  group roughly matches the sum of the R2-resolved individual strain counts in
  the same sample -- large mismatches get flagged
- Cross-primer-pair matching is valid (confirmed): 341F1-F6 / 799F1-F2 all
  target the same conserved amplicon region, so a reference built from one
  primer pair's reads is comparable to ASVs from any other primer pair after
  trimming.

Matching: exact sequence match first (reference already includes known
operon-copy variants). Falls back to closest match within
FUZZY_MATCH_MAX_DISTANCE for non-exact hits -- flagged as "close_match".
Anything beyond that distance to every reference strain is reported as 
unassigned rather than forced into a guess.

Usage:
    python assign_strain_identity.py \
        --metadata combined_master_metadata.csv \
        --seqtab-r1 group1_R1_long.csv group2_R1_long.csv ... \
        --seqtab-r2 group1_R2_long.csv group2_R2_long.csv ... \
        --ref-r1 monoculture_reference_R1.csv \
        --ref-r2 monoculture_reference_R2.csv \
        --output-dir strain_assignment
"""

import argparse
import os
from collections import defaultdict

import pandas as pd

FUZZY_MATCH_MAX_DISTANCE = 3  # same threshold used when building the reference
QC_MISMATCH_FLAG_RATIO = 0.20  # flag if |R1 combined - R2 summed| exceeds this
                                # fraction of the R1 combined count
SELF_MATCH_FLAG_THRESHOLD = 0.90  # monoculture wells below this self-match
                                   # fraction get flagged for review


def edit_distance(a: str, b: str) -> int:
    """Full (unbounded) Levenshtein distance -- kept for cases where the exact
    value matters. NOT used in the hot matching loop below; see
    edit_distance_within(), which is what match_asv() actually calls."""
    if a == b:
        return 0
    la, lb = len(a), len(b)
    prev = list(range(lb + 1))
    for i in range(1, la + 1):
        curr = [i] + [0] * lb
        ai = a[i - 1]
        for j in range(1, lb + 1):
            cost = 0 if ai == b[j - 1] else 1
            curr[j] = min(prev[j] + 1, curr[j - 1] + 1, prev[j - 1] + cost)
        prev = curr
    return prev[lb]


def edit_distance_within(a: str, b: str, max_dist: int) -> int:
    """
    Banded Levenshtein distance: only answers 'is the distance <= max_dist',
    returning max_dist + 1 (a sentinel meaning 'too far, exact value doesn't
    matter') if not. Restricting the DP to a diagonal band of width
    2*max_dist+1 -- instead of the full L x L grid -- turns each comparison
    from O(L^2) into roughly O(L * max_dist), with an immediate length-based
    bailout before any DP at all.
    """
    la, lb = len(a), len(b)
    if abs(la - lb) > max_dist:
        return max_dist + 1  # length alone rules it out -- skip the DP entirely

    INF = max_dist + 1
    prev = [INF] * (lb + 1)
    for j in range(0, min(lb, max_dist) + 1):
        prev[j] = j
    for i in range(1, la + 1):
        curr = [INF] * (lb + 1)
        lo = max(0, i - max_dist)
        hi = min(lb, i + max_dist)
        if lo == 0:
            curr[0] = i
        ai = a[i - 1]
        for j in range(max(1, lo), hi + 1):
            cost = 0 if ai == b[j - 1] else 1
            sub = prev[j] + 1 if prev[j] < INF else INF
            ins = curr[j - 1] + 1 if curr[j - 1] < INF else INF
            diag = prev[j - 1] + cost if prev[j - 1] < INF else INF
            curr[j] = min(sub, ins, diag)
        row_min = min((v for v in curr[lo:hi + 1]), default=INF)
        if row_min > max_dist:
            return INF  # whole band already exceeds threshold -- bail out early
        prev = curr
    return prev[lb]


def load_reference(ref_csv_path: str):
    """
    Returns:
      exact_lookup: {sequence: set of STRAINs} -- a set because reference
        variants could (rarely) coincide across strains; size >1 means this
        exact sequence is ambiguous between those strains.
      all_ref_seqs: list of (STRAIN, sequence) for fuzzy fallback search.
    """
    ref = pd.read_csv(ref_csv_path)
    exact_lookup = defaultdict(set)
    for _, row in ref.iterrows():
        exact_lookup[row["sequence"]].add(row["STRAIN"])
    all_ref_seqs = list(zip(ref["STRAIN"], ref["sequence"]))
    return dict(exact_lookup), all_ref_seqs


def match_asv(sequence: str, exact_lookup: dict, all_ref_seqs: list):
    """
    Returns (call, match_type, distance) where:
      call = frozenset of STRAIN(s) -- singleton if unambiguous, multiple if
             genuinely ambiguous at this sequence
      match_type = 'exact' | 'close_match' | 'unassigned'
      distance = 0 for exact, else the edit distance to the nearest ref seq(s)
    """
    if sequence in exact_lookup:
        return frozenset(exact_lookup[sequence]), "exact", 0

    # Fuzzy fallback: find the minimum distance to ANY reference sequence,
    # then collect all strains achieving that minimum (handles ties honestly
    # instead of picking one arbitrarily). Uses the banded distance function
    # (fast bail-out past FUZZY_MATCH_MAX_DISTANCE) since we only need to know
    # "is this within threshold", not the exact distance for far-away candidates
    best_dist = FUZZY_MATCH_MAX_DISTANCE + 1  # sentinel: "nothing within threshold yet"
    best_strains = set()
    for strain, ref_seq in all_ref_seqs:
        d = edit_distance_within(sequence, ref_seq, FUZZY_MATCH_MAX_DISTANCE)
        if d < best_dist:
            best_dist = d
            best_strains = {strain}
        elif d == best_dist and d <= FUZZY_MATCH_MAX_DISTANCE:
            best_strains.add(strain)

    if best_dist <= FUZZY_MATCH_MAX_DISTANCE:
        return frozenset(best_strains), "close_match", best_dist
    return frozenset(), "unassigned", None  # exact distance not computed past threshold -- not needed for the decision


def load_and_concat_seqtabs(paths: list) -> pd.DataFrame:
    dfs = []
    for p in paths:
        df = pd.read_csv(p)
        required = {"sample_uid", "sequence", "count"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"{p} missing columns: {missing}")
        dfs.append(df)
    combined = pd.concat(dfs, ignore_index=True)

    # A given sample_uid should only appear in ONE input file (one demux
    # group). If it shows up with different ASV sets across files that's
    # fine (expected -- same sample listed once, many ASV rows); but if the
    # exact same (sample_uid, sequence) pair appears more than once across
    # files, something's duplicated upstream -- flag
    dupes = combined[combined.duplicated(subset=["sample_uid", "sequence"], keep=False)]
    if len(dupes) > 0:
        raise RuntimeError(
            f"{len(dupes)} duplicate (sample_uid, sequence) rows found across "
            f"the provided seqtab files -- likely the same demux group's file "
            f"was passed in twice, or two files unexpectedly cover overlapping "
            f"samples:\n{dupes[['sample_uid','sequence','count']].head(20).to_string()}"
        )
    return combined


def assign_direction(seqtab: pd.DataFrame, exact_lookup: dict, all_ref_seqs: list, direction: str):
    """Returns (per_asv_calls_df, sample_strain_abundance_df)."""
    unique_seqs = seqtab["sequence"].unique()
    print(f"  Matching {len(unique_seqs)} unique {direction} ASVs against reference...")

    match_cache = {}
    for seq in unique_seqs:
        match_cache[seq] = match_asv(seq, exact_lookup, all_ref_seqs)

    seqtab = seqtab.copy()
    seqtab["call"] = seqtab["sequence"].map(lambda s: match_cache[s][0])
    seqtab["match_type"] = seqtab["sequence"].map(lambda s: match_cache[s][1])
    seqtab["match_distance"] = seqtab["sequence"].map(lambda s: match_cache[s][2])
    seqtab["n_strains_in_call"] = seqtab["call"].map(len)
    seqtab["call_label"] = seqtab["call"].map(
        lambda s: "+".join(sorted(s)) if s else "unassigned"
    )

    n_exact = (seqtab["match_type"] == "exact").sum()
    n_close = (seqtab["match_type"] == "close_match").sum()
    n_unassigned = (seqtab["match_type"] == "unassigned").sum()
    print(f"    {n_exact} exact, {n_close} close_match, {n_unassigned} unassigned "
          f"(row counts, not unique ASVs)")

    # Sample x call_label abundance (call_label may be a single strain, an
    # ambiguous multi-strain group like "59+p25C2", or "unassigned").
    abundance = (
        seqtab.groupby(["sample_uid", "call_label"], as_index=False)["count"].sum()
    )
    return seqtab, abundance


def resolve_ambiguous_groups(r1_abundance: pd.DataFrame, r2_abundance: pd.DataFrame):
    """
    Final strain-abundance table, per sample:
    - Any call_label from R2 that's a single strain (no '+') is used directly.
    - Any R2 call_label that's still an ambiguous multi-strain group (R2 failed
      to resolve it further) is left as that group, flagged, since it means R2
      itself couldn't separate those strains for this particular ASV either.
    - R1 combined counts for each ambiguous group are cross-checked against the
      SUM of R2-resolved individual-strain counts for that same group/sample --
      large mismatches are flagged, not silently ignored.

    Returns (final_abundance_df, qc_mismatch_df)
    """
    final = r2_abundance.rename(columns={"call_label": "strain_or_group", "count": "R2_count"})

    # Build R1 ambiguous-group totals per sample: any R1 call_label containing
    # '+' represents a group R1 alone can't split.
    r1_groups = r1_abundance[r1_abundance["call_label"].str.contains(r"\+", na=False)].copy()
    r1_groups = r1_groups.rename(columns={"call_label": "group_label", "count": "R1_combined_count"})

    qc_rows = []
    for _, row in r1_groups.iterrows():
        sample_uid, group_label, r1_count = row["sample_uid"], row["group_label"], row["R1_combined_count"]
        member_strains = set(group_label.split("+"))

        # Sum R2-resolved counts for the individual member strains in this sample
        r2_sub = final[
            (final["sample_uid"] == sample_uid) &
            (final["strain_or_group"].isin(member_strains))
        ]
        r2_sum = r2_sub["R2_count"].sum()

        denom = r1_count if r1_count > 0 else 1
        rel_diff = abs(r1_count - r2_sum) / denom
        flagged = rel_diff > QC_MISMATCH_FLAG_RATIO

        qc_rows.append({
            "sample_uid": sample_uid,
            "ambiguous_group": group_label,
            "R1_combined_count": r1_count,
            "R2_resolved_sum": r2_sum,
            "relative_difference": round(rel_diff, 4),
            "flagged_mismatch": flagged,
        })
        if flagged:
            print(
                f"  *** QC FLAG: sample {sample_uid}, group {group_label}: "
                f"R1 combined={r1_count} vs R2 resolved sum={r2_sum} "
                f"({rel_diff:.1%} relative difference) ***"
            )

    qc_df = pd.DataFrame(qc_rows)
    return final, qc_df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--metadata", required=True)
    p.add_argument("--seqtab-r1", nargs="+", required=True)
    p.add_argument("--seqtab-r2", nargs="+", required=True)
    p.add_argument("--ref-r1", required=True)
    p.add_argument("--ref-r2", required=True)
    p.add_argument("--output-dir", required=True)
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Loading reference sequences...")
    exact_r1, allref_r1 = load_reference(args.ref_r1)
    exact_r2, allref_r2 = load_reference(args.ref_r2)
    print(f"  R1 reference: {len(allref_r1)} sequences, R2 reference: {len(allref_r2)} sequences")

    print("\nLoading metadata...")
    # dtype=str on STRAIN is deliberate, some strains are numeric IDs
    meta = pd.read_csv(args.metadata, dtype={"STRAIN": str})
    if "sample_type" not in meta.columns:
        raise ValueError("Metadata missing 'sample_type' column")

    # ALL samples -- including monoculture -- go through the main assignment
    # pass now. Monoculture wells are the reference source, so a monoculture
    # well matching its own strain is circular and not informative on its
    # own; but off-target/contaminating/unassigned reads within a monoculture
    # well ARE informative and previously got dropped by excluding these
    # wells outright. They're kept flagged (via monoculture_uids below) so
    # they can still be pulled out separately for the self-check summary
    # without re-running the matching a second time.
    community_uids = set(meta["sample_uid"])
    monoculture_uids = set(meta[meta["sample_type"] == "monoculture"]["sample_uid"])
    print(f"  {len(community_uids)}/{len(meta)} metadata rows carried through the main assignment pass")
    print(f"  ({len(monoculture_uids)} of those are monoculture wells, flagged for the self-check summary below)")

    print("\nLoading and concatenating seqtabs...")
    seqtab_r1 = load_and_concat_seqtabs(args.seqtab_r1)
    seqtab_r2 = load_and_concat_seqtabs(args.seqtab_r2)

    n_r1_before, n_r2_before = len(seqtab_r1), len(seqtab_r2)
    seqtab_r1 = seqtab_r1[seqtab_r1["sample_uid"].isin(community_uids)]
    seqtab_r2 = seqtab_r2[seqtab_r2["sample_uid"].isin(community_uids)]
    print(f"  R1: {len(seqtab_r1)}/{n_r1_before} rows kept, {seqtab_r1['sample_uid'].nunique()} samples")
    print(f"  R2: {len(seqtab_r2)}/{n_r2_before} rows kept, {seqtab_r2['sample_uid'].nunique()} samples")

    print("\n=== R1 assignment (QC cross-check role) ===")
    r1_calls, r1_abundance = assign_direction(seqtab_r1, exact_r1, allref_r1, "R1")

    print("\n=== R2 assignment (primary strain identity) ===")
    r2_calls, r2_abundance = assign_direction(seqtab_r2, exact_r2, allref_r2, "R2")

    print("\n=== Resolving ambiguous groups + QC cross-check ===")
    final_abundance, qc_mismatch = resolve_ambiguous_groups(r1_abundance, r2_abundance)

    # --- Outputs ---
    final_path = os.path.join(args.output_dir, "community_strain_abundance.csv")
    final_abundance.to_csv(final_path, index=False)
    print(f"\nWrote final strain abundance table (now includes monoculture wells): {final_path}")

    r1_calls_path = os.path.join(args.output_dir, "r1_asv_calls_detail.csv")
    r2_calls_path = os.path.join(args.output_dir, "r2_asv_calls_detail.csv")
    r1_calls.to_csv(r1_calls_path, index=False)
    r2_calls.to_csv(r2_calls_path, index=False)
    print(f"Wrote per-ASV call detail: {r1_calls_path}, {r2_calls_path}")

    if len(qc_mismatch) > 0:
        qc_path = os.path.join(args.output_dir, "ambiguous_group_qc.csv")
        qc_mismatch.to_csv(qc_path, index=False)
        n_flagged = qc_mismatch["flagged_mismatch"].sum()
        print(f"Wrote ambiguous-group QC report: {qc_path} "
              f"({n_flagged}/{len(qc_mismatch)} flagged for mismatch)")

    unassigned_r2 = r2_calls[r2_calls["match_type"] == "unassigned"]
    if len(unassigned_r2) > 0:
        unassigned_summary = (
            unassigned_r2.groupby("sequence", as_index=False)["count"].sum()
            .sort_values("count", ascending=False)
        )
        unassigned_path = os.path.join(args.output_dir, "unassigned_r2_asvs.csv")
        unassigned_summary.to_csv(unassigned_path, index=False)
        print(f"\n*** {len(unassigned_summary)} distinct R2 ASVs unassigned to any "
              f"reference strain (total {unassigned_r2['count'].sum()} reads) -- "
              f"see {unassigned_path} for review. ***")

    # --- Monoculture self-check summary ---
    # Monoculture wells already went through the exact same R1/R2 pipeline as
    # every other sample above and are already in final_abundance. This 
    # pulls them back out and summarizes self-match fraction for a quick
    # pipeline sanity check -- it does NOT re-run assign_direction.
    if monoculture_uids:
        print("\n=== Monoculture self-check summary (subset of final_abundance above) ===")
        mono_final = final_abundance[final_abundance["sample_uid"].isin(monoculture_uids)].copy()

        mono_meta = meta[meta["sample_type"] == "monoculture"][["sample_uid", "STRAIN"]]
        mono_check = mono_final.merge(mono_meta, on="sample_uid", how="left")
        mono_check["self_match"] = mono_check["strain_or_group"] == mono_check["STRAIN"].astype(str)

        totals = mono_check.groupby("sample_uid")["R2_count"].transform("sum")
        mono_check["fraction_of_well"] = mono_check["R2_count"] / totals

        self_check_summary = (
            mono_check[mono_check["self_match"]]
            .groupby("sample_uid", as_index=False)
            .agg(expected_strain=("STRAIN", "first"),
                 self_match_fraction=("fraction_of_well", "sum"))
        )
        # Wells with zero self-matching rows (e.g. entirely unassigned/off-target)
        # won't appear in self_check_summary above -- add them back at 0% so
        # they don't just silently vanish from the report.
        missing_wells = monoculture_uids - set(self_check_summary["sample_uid"])
        if missing_wells:
            fill_rows = mono_meta[mono_meta["sample_uid"].isin(missing_wells)].rename(
                columns={"STRAIN": "expected_strain"}
            )
            fill_rows["self_match_fraction"] = 0.0
            self_check_summary = pd.concat([self_check_summary, fill_rows], ignore_index=True)

        low_self_match = self_check_summary[self_check_summary["self_match_fraction"] < SELF_MATCH_FLAG_THRESHOLD]
        if len(low_self_match) > 0:
            print(f"  *** {len(low_self_match)} monoculture wells self-assign to their "
                  f"own strain at <{SELF_MATCH_FLAG_THRESHOLD:.0%} -- pipeline/reference issue worth investigating: ***")
            print(low_self_match.to_string(index=False))
        else:
            print(f"  All {len(self_check_summary)} monoculture wells self-assigned "
                  f"cleanly (>={SELF_MATCH_FLAG_THRESHOLD:.0%} to their own strain).")

        mono_path = os.path.join(args.output_dir, "monoculture_self_check.csv")
        mono_check.to_csv(mono_path, index=False)
        print(f"  Wrote {mono_path}")


if __name__ == "__main__":
    main()
