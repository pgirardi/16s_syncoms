"""
Build master sample metadata for the IN PLANTA trials (1 & 2) from
growth_coll_plate_mapping.xlsx.

Replaces 00_build_master_metadata.R for the in planta portion.
"""

import re
import openpyxl
import pandas as pd

from primer_lookup import lookup_in_planta_primer, fastq_paths

XLSX_PATH = "growth_coll_plate_mapping.xlsx"

FASTQ_BASE_DIR_TEMPLATE = "../../demulti_data_trial{seq_run}"


def load_sheet_as_df(path, sheet_name):
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb[sheet_name]
    rows = list(ws.iter_rows(values_only=True))
    header, data = rows[0], rows[1:]
    df = pd.DataFrame(data, columns=header)
    return df.dropna(how="all")


def classify_sample_type(growth_plate, growth_well):
    """Classify a collection_map row based on what's in growth_plate/growth_well."""
    gp = str(growth_plate).strip() if growth_plate is not None else "NA"
    gw = str(growth_well).strip() if growth_well is not None else "NA"

    if gp == "NA" and gw == "NA":
        return "empty"
    if gp == "Extraction Blank":
        return "extraction_blank"
    if re.match(r"^Plate\d+$", gp):
        return "experimental"
    if gp == "p25C2":
        return "monoculture"
    if re.match(r"^\d+(\.\d+)?$", gp):
        return "monoculture"
    if gp.startswith("t0"):
        # 't0 1' -> t0, 't0 1 L' -> t0_L, 't0 1C' -> t0C, 't0 1C L' -> t0C_L
        norm = gp.replace(" ", "")
        return f"control_{norm}"
    return "unclassified"  # deliberately NOT silently dropped -- surfaced in QC report


def normalize_strain_id(value):
    """Normalize a growth_plate value into a STRAIN id matching the format used
    everywhere else in the pipeline (monoculture reference, community_strain_
    abundance.csv, in vitro's STRAIN column)"""
    s = str(value).strip()
    if re.match(r"^\d+\.0$", s):
        s = s[:-2]
    return s


def build():
    growth_layout = load_sheet_as_df(XLSX_PATH, "growth_layout")
    collection_map = load_sheet_as_df(XLSX_PATH, "collection_map")

    growth_layout.columns = [str(c).strip() for c in growth_layout.columns]
    collection_map.columns = [str(c).strip() for c in collection_map.columns]

    records = []
    lookup_failures = []

    for _, row in collection_map.iterrows():
        coll_plate = str(row["collection_plate"]).strip()
        coll_well = str(row["collection_well"]).strip()
        growth_plate = row["growth_plate"]
        growth_well = row["growth_well"]
        trial = row["trial"]

        sample_type = classify_sample_type(growth_plate, growth_well)
        strain = normalize_strain_id(growth_plate) if sample_type == "monoculture" else pd.NA

        # Primer/seq_run lookup -- hard fail (not silent NA) on any miss.
        try:
            primer_info = lookup_in_planta_primer(coll_plate, coll_well)
        except (KeyError, ValueError) as e:
            lookup_failures.append({
                "collection_plate": coll_plate,
                "collection_well": coll_well,
                "error": str(e),
            })
            continue

        r1, r2 = fastq_paths(
            primer_info["primer_pair"], primer_info["seq_run"], coll_well,
            FASTQ_BASE_DIR_TEMPLATE,
        )

        rec = {
            "sample_uid": f"{coll_plate}_{coll_well}",
            "collection_plate": coll_plate,
            "collection_well": coll_well,
            "plate_label": primer_info["plate_label"],
            "primer_pair": primer_info["primer_pair"],
            "seq_run": primer_info["seq_run"],
            "sample_type": sample_type,
            "STRAIN": strain,
            "growth_plate_raw": growth_plate,
            "growth_well_raw": growth_well,
            "trial": trial,
            "fastq_R1": r1,
            "fastq_R2": r2,
        }
        records.append(rec)

    df = pd.DataFrame(records)

    gl = growth_layout.rename(columns={
        "growth_plate": "growth_plate_raw",
        "growth_well": "growth_well_raw",
    })

    dropped_col = [c for c in gl.columns if str(c).strip().lower() == "nan"]
    if dropped_col:
        gl = gl.rename(columns={dropped_col[0]: "dropped_flag_raw"})
    else:
        gl["dropped_flag_raw"] = None

    df["trial"] = pd.to_numeric(df["trial"], errors="coerce")
    gl["trial"] = pd.to_numeric(gl["trial"], errors="coerce")

    df = df.merge(
        gl[["growth_plate_raw", "growth_well_raw", "trial", "community", "treatment",
            "dpi", "dropped_flag_raw"]],
        on=["growth_plate_raw", "growth_well_raw", "trial"],
        how="left",
    )

    # --- HARD VALIDATION (the stopifnot() the old pipeline was missing) ---
    if lookup_failures:
        msg = "\n".join(
            f"  {f['collection_plate']} {f['collection_well']}: {f['error']}"
            for f in lookup_failures
        )
        raise RuntimeError(f"Primer lookup failed for {len(lookup_failures)} wells:\n{msg}")

    exp_mask = df["sample_type"] == "experimental"
    missing_community = df[exp_mask & df["community"].isna()]
    if len(missing_community) > 0:
        raise RuntimeError(
            f"{len(missing_community)} rows classified 'experimental' but failed to "
            f"join to growth_layout (community is NaN):\n"
            f"{missing_community[['sample_uid','growth_plate_raw','growth_well_raw']]}"
        )

    unclassified = df[df["sample_type"] == "unclassified"]
    if len(unclassified) > 0:
        print("WARNING -- unclassified rows (need manual review):")
        print(unclassified[["sample_uid", "growth_plate_raw", "growth_well_raw"]])

    monoculture_missing_strain = df[
        (df["sample_type"] == "monoculture") & (df["STRAIN"].isna())
    ]
    if len(monoculture_missing_strain) > 0:
        raise RuntimeError(
            f"{len(monoculture_missing_strain)} rows classified 'monoculture' but "
            f"have no STRAIN value -- growth_plate did not resolve to a strain ID:\n"
            f"{monoculture_missing_strain[['sample_uid','growth_plate_raw']].to_string()}"
        )

    # --- dropped_flag semantics  ---
    # dropped_flag == "dropped" means: NOT actually inoculated with culture (plant
    # looked unhealthy, so inoculation was skipped) -- but the plant WAS still
    # collected and sequenced. So these are real samples with real fastq data,
    # but the community/treatment values from growth_layout describe the INTENDED
    # design, not what actually happened. They must not be silently included in
    # PERMANOVA/CLR as if inoculation occurred.
    #
    # We keep them in the metadata (real data, potentially useful e.g. as
    # uninoculated-plant reference points) but flag them explicitly so downstream
    # scripts can filter with one clear boolean rather than re-deriving this logic.
    df["truly_inoculated"] = (~(
        (df["dropped_flag_raw"] == "dropped") & (df["sample_type"] == "experimental")
    )).astype("boolean")  # pandas nullable boolean dtype, so pd.NA assignment below works
    # Only meaningful for experimental rows; leave NaN for everything else so it's
    # not misread as "not inoculated" for e.g. monoculture/control/blank rows.
    df.loc[df["sample_type"] != "experimental", "truly_inoculated"] = pd.NA

    n_not_inoculated = (df["truly_inoculated"] == False).sum()  # noqa: E712
    if n_not_inoculated > 0:
        print(
            f"\nNOTE: {n_not_inoculated} experimental-labeled rows are flagged "
            f"truly_inoculated=False (collected/sequenced but NOT actually "
            f"inoculated -- exclude from PERMANOVA/CLR via this column, don't drop "
            f"the rows outright)."
        )

    return df


if __name__ == "__main__":
    df = build()
    print(f"\nTotal rows: {len(df)}")
    print("\nBy sample_type:")
    print(df["sample_type"].value_counts())
    print("\nMonoculture STRAIN values:")
    print(df[df["sample_type"] == "monoculture"][["sample_uid", "STRAIN"]].to_string())
    print("\nBy plate_label:")
    print(df["plate_label"].value_counts())
    print("\nBy trial (experimental only):")
    print(df[df["sample_type"] == "experimental"]["trial"].value_counts())

    out_path = "./in_planta_master_metadata.csv"
    df.to_csv(out_path, index=False)
    print(f"\nWrote {out_path}")
    print(df.head(10).to_string())
