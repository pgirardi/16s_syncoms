"""
Build reference FASTA sequences (both R1 and R2) for each strain from the
IN VITRO monoculture wells -- confirmed as the cleanest, most reliable source
for reference sequences (single strain per well, no community mixing).

Both directions are built because R1 alone can't reliably distinguish Leaf59
from p25C2 (only 2 SNPs across the R1 window); R2 carries 4 distinguishing
SNPs and is the intended basis for that disambiguation.

For each strain, in each direction independently: pools read counts across
all monoculture wells for that strain, then takes the single most abundant
ASV as the reference sequence. Flags (does NOT silently accept) any strain
whose top ASV doesn't clearly dominate -- that's a sign of cross-contamination
or an unexpectedly diverse/wrong well, and needs a human look before being
trusted as a reference.

Usage:
    python build_monoculture_reference_fasta.py \
        --metadata in_vitro_master_metadata.csv \
        --seqtab-r1 seqtab_R1_long.csv \
        --seqtab-r2 seqtab_R2_long.csv \
        --output-dir monoculture_reference

Expects seqtab CSVs in long format (sample_uid, sequence, count) -- i.e. the
output of export_seqtab_to_csv.R.
"""

import argparse
import os
import sys

import pandas as pd

DOMINANCE_FLAG_THRESHOLD = 0.90  # combined dominance (top + included operon-
                                  # copy variants) should carry at least this
                                  # fraction of a strain's total reads to be
                                  # trusted without manual review

NEAR_IDENTICAL_EDIT_DISTANCE = 3  # bp -- ASVs within this distance of the top
                                   # sequence are treated as real intragenomic
                                   # operon-copy variants of the SAME strain,
                                   # not noise/contamination, and are included
                                   # in the reference. Calibrated against real
                                   # data: observed R2 top-vs-second distances
                                   # were exactly 1bp for every strain, while
                                   # genuinely separate (noise/background)
                                   # ASVs sat 32-47bp away -- wide margin
                                   # either side of this threshold.


def edit_distance(a: str, b: str) -> int:
    """Simple Levenshtein distance, no external dependency. Fine for the
    small number of comparisons here (one per strain/direction)."""
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


def load_monoculture_strain_map(metadata_path: str) -> pd.DataFrame:
    meta = pd.read_csv(metadata_path)
    required = {"sample_uid", "sample_type", "STRAIN"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"Metadata missing required columns: {missing}")

    mono = meta[meta["sample_type"] == "monoculture"].copy()
    if mono.empty:
        raise ValueError(
            "No rows with sample_type == 'monoculture' found in metadata -- "
            "check you passed the in vitro metadata CSV, not in planta."
        )
    if mono["STRAIN"].isna().any():
        bad = mono[mono["STRAIN"].isna()]["sample_uid"].tolist()
        raise ValueError(
            f"{len(bad)} monoculture rows have no STRAIN value, can't build "
            f"a reference from them: {bad}"
        )
    return mono[["sample_uid", "STRAIN"]]


def build_reference_for_direction(
    strain_map: pd.DataFrame, seqtab_long: pd.DataFrame, direction: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns (reference_df, qc_df).
    reference_df: one row PER INCLUDED VARIANT (top ASV, plus any ASV within
        NEAR_IDENTICAL_EDIT_DISTANCE of it -- treated as the same strain's own
        operon-copy variants, not separate things).
    qc_df: one row per strain -- summary diagnostics including flagged rows.
    """
    required = {"sample_uid", "sequence", "count"}
    missing = required - set(seqtab_long.columns)
    if missing:
        raise ValueError(f"seqtab CSV ({direction}) missing columns: {missing}")

    merged = strain_map.merge(seqtab_long, on="sample_uid", how="left")

    missing_reads = merged[merged["sequence"].isna()]["sample_uid"].unique()
    if len(missing_reads) > 0:
        print(
            f"  WARNING ({direction}): {len(missing_reads)} monoculture wells have "
            f"NO reads in this seqtab at all (missing from DADA2 output entirely, "
            f"not just low-count): {list(missing_reads)}"
        )
    merged = merged.dropna(subset=["sequence"])

    # Pool counts across all wells for a given strain (usually 1 well/strain,
    # but this handles replicates correctly if present) and per unique ASV.
    pooled = (
        merged.groupby(["STRAIN", "sequence"], as_index=False)["count"].sum()
    )

    ref_rows = []
    qc_rows = []
    for strain, grp in pooled.groupby("STRAIN"):
        grp_sorted = grp.sort_values("count", ascending=False).reset_index(drop=True)
        total_reads = grp_sorted["count"].sum()
        top = grp_sorted.iloc[0]
        n_asvs = len(grp_sorted)

        # Walk down from the top ASV, including any ASV within
        # NEAR_IDENTICAL_EDIT_DISTANCE as a real operon-copy variant of this
        # strain. Stops at the first ASV that's too far away -- doesn't keep
        # searching past it, since count-ranking + a hard distance cutoff is
        # enough to separate "same strain, different rRNA copy" from
        # "unrelated background sequence" given the wide gap observed in
        # real data (1bp vs 32-47bp).
        included = [{"sequence": top["sequence"], "count": top["count"], "edit_distance_from_top": 0}]
        for i in range(1, n_asvs):
            cand = grp_sorted.iloc[i]
            dist = edit_distance(top["sequence"], cand["sequence"])
            if dist <= NEAR_IDENTICAL_EDIT_DISTANCE:
                included.append({
                    "sequence": cand["sequence"], "count": cand["count"],
                    "edit_distance_from_top": dist,
                })

        combined_count = sum(v["count"] for v in included)
        combined_dominance = combined_count / total_reads if total_reads > 0 else 0.0
        flagged = combined_dominance < DOMINANCE_FLAG_THRESHOLD

        for rank, variant in enumerate(included, start=1):
            ref_rows.append({
                "STRAIN": strain,
                "variant_rank": rank,
                "sequence": variant["sequence"],
                "seq_length": len(variant["sequence"]),
                "read_count": variant["count"],
                "edit_distance_from_top": variant["edit_distance_from_top"],
                "total_reads_all_asvs": total_reads,
                "combined_dominance_fraction": round(combined_dominance, 4),
                "flagged_low_dominance": flagged,
            })

        # Distance from top to the very next-ranked ASV regardless of whether
        # it got included -- kept as a diagnostic column so you can still see
        # what's just outside the cutoff, not just what made it in.
        next_dist = None
        next_count = 0
        if n_asvs > 1:
            nxt = grp_sorted.iloc[1]
            next_dist = edit_distance(top["sequence"], nxt["sequence"])
            next_count = nxt["count"]

        qc_rows.append({
            "STRAIN": strain,
            "direction": direction,
            "total_reads": total_reads,
            "n_distinct_asvs": n_asvs,
            "n_variants_included": len(included),
            "combined_dominance_fraction": round(combined_dominance, 4),
            "flagged_low_dominance": flagged,
            "top_vs_next_edit_distance": next_dist,
            "next_asv_count": next_count,
        })
        if flagged:
            print(
                f"  *** FLAG ({direction}) strain {strain!r}: top ASV(s) within "
                f"{NEAR_IDENTICAL_EDIT_DISTANCE}bp of each other only cover "
                f"{combined_dominance:.1%} of {total_reads} reads across "
                f"{n_asvs} distinct ASVs -- review before trusting as reference "
                f"(possible cross-contamination or mixed well). ***"
            )

    return pd.DataFrame(ref_rows), pd.DataFrame(qc_rows)


def write_fasta(reference_df: pd.DataFrame, out_path: str, direction: str):
    variants_per_strain = reference_df.groupby("STRAIN")["variant_rank"].transform("max")
    with open(out_path, "w") as f:
        for (_, row), n_variants in zip(reference_df.iterrows(), variants_per_strain):
            suffix = f"_v{row['variant_rank']}" if n_variants > 1 else ""
            header = f">{row['STRAIN']}_{direction}{suffix}"
            f.write(f"{header}\n{row['sequence']}\n")
    print(f"Wrote {len(reference_df)} sequences to {out_path} "
          f"({reference_df['STRAIN'].nunique()} strains)")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--metadata", required=True)
    p.add_argument("--seqtab-r1", required=True)
    p.add_argument("--seqtab-r2", required=True)
    p.add_argument("--output-dir", required=True)
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    strain_map = load_monoculture_strain_map(args.metadata)
    print(f"Found {len(strain_map)} monoculture wells covering "
          f"{strain_map['STRAIN'].nunique()} distinct strains.")

    seqtab_r1 = pd.read_csv(args.seqtab_r1)
    seqtab_r2 = pd.read_csv(args.seqtab_r2)

    all_qc = []
    for direction, seqtab in [("R1", seqtab_r1), ("R2", seqtab_r2)]:
        print(f"\n=== Building {direction} reference ===")
        ref_df, qc_df = build_reference_for_direction(strain_map, seqtab, direction)

        fasta_path = os.path.join(args.output_dir, f"monoculture_reference_{direction}.fasta")
        write_fasta(ref_df, fasta_path, direction)

        ref_csv_path = os.path.join(args.output_dir, f"monoculture_reference_{direction}.csv")
        ref_df.to_csv(ref_csv_path, index=False)

        all_qc.append(qc_df)

    qc_combined = pd.concat(all_qc, ignore_index=True)
    qc_path = os.path.join(args.output_dir, "monoculture_reference_qc.csv")
    qc_combined.to_csv(qc_path, index=False)
    print(f"\nWrote combined QC report: {qc_path}")

    n_flagged = qc_combined["flagged_low_dominance"].sum()
    n_multi_variant = (qc_combined["n_variants_included"] > 1).sum()
    print(f"\n{n_multi_variant}/{len(qc_combined)} strain/direction combinations "
          f"included more than one operon-copy variant in the reference.")
    if n_flagged > 0:
        print(
            f"\n*** {n_flagged} strain/direction combinations still flagged for low "
            f"combined dominance even after including near-identical variants -- "
            f"review monoculture_reference_qc.csv, these look more like genuine "
            f"contamination than operon copy variation. ***"
        )
    else:
        print("\nAll strains passed the combined-dominance check cleanly in both directions.")


if __name__ == "__main__":
    main()
