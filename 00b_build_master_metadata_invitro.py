"""
Build master sample metadata for the IN VITRO trial from in_vitro_plate_layout.xlsx.

Design:
- One 56-well treatment/replicate layout (Sample Metadata sheet, WELL A1-H7) is
  sampled REPEATEDLY (non-destructive -- same physical cultures over time) at
  4 timepoints: 6H, 12H, 18H, 24H.
- Each timepoint's aliquots are plated into a different well range, spread
  across 3 sequencing plates (12H is split across the tail of plate 1 and the
  head of plate 2).
- Plate 3 additionally carries monoculture controls (cols 9-10) and t0
  references (S1 t0 col 11, S1C t0 col 12) -- NOT part of the repeated 56-well
  design, so handled as a separate category.
- Every plate also has dedicated BLANKS well-columns.

Well-range boundaries below were derived from the merged-cell structure of the
'Sequencing Plate Layout' sheet in the source xlsx
"""

import re
import openpyxl
import pandas as pd

from primer_lookup import IN_VITRO_PRIMER_TABLE, fastq_paths

XLSX_PATH = "in_vitro_plate_layout.xlsx"
FASTQ_BASE_DIR_TEMPLATE = "../../demulti_data_trial{seq_run}"

SHEET_LAYOUT = "Sequencing Plate Layout"
SHEET_DESIGN = "Sample Metadata"

# (plate_num, timepoint_label, well_col_start, well_col_end) -- order matters for
# timepoints split across plates (12H: plate1 tail then plate2 head).
TIMEPOINT_BLOCKS = [
    (1, "6H",  1, 7),
    (1, "12H", 9, 12),
    (2, "12H", 1, 3),
    (2, "18H", 5, 11),
    (3, "24H", 1, 7),
]

BLANK_COLUMNS = [(1, 8), (2, 4), (2, 12), (3, 8)]

# Plate 3 special region: (well_col_start, well_col_end) -> category
PLATE3_SPECIAL_COLS = {
    (9, 10): "monoculture",
    (11, 11): "control_S1_t0",
    (12, 12): "control_S1C_t0",
}

ROWS = list("ABCDEFGH")


def load_design_table():
    """56-well treatment/replicate design from Sample Metadata sheet."""
    wb = openpyxl.load_workbook(XLSX_PATH, data_only=True)
    ws = wb[SHEET_DESIGN]
    rows = list(ws.iter_rows(values_only=True))
    header, data = rows[0], rows[1:]
    df = pd.DataFrame(data, columns=header).dropna(how="all")
    df.columns = [str(c).strip() for c in df.columns]

    def parse_design_well(w):
        m = re.match(r"^([A-H])(\d)$", str(w).strip())
        if not m:
            raise ValueError(f"Unexpected design WELL format: {w!r}")
        return m.group(1), int(m.group(2))

    parsed = df["WELL"].apply(parse_design_well)
    df["design_row"] = parsed.apply(lambda x: x[0])
    df["design_col"] = parsed.apply(lambda x: x[1])

    n_expected = 8 * 7
    if len(df) != n_expected:
        raise RuntimeError(
            f"Expected {n_expected} rows (8 rows x 7 cols) in Sample Metadata, got {len(df)}"
        )
    return df


def build_timepoint_col_map():
    """
    Flattens TIMEPOINT_BLOCKS into: {timepoint: {design_col (1-7): (plate, well_col)}}
    Concatenates blocks IN ORDER per timepoint (e.g. 12H: plate1 cols9-12 then
    plate2 cols1-3), so design_col 1-4 -> plate1, design_col 5-7 -> plate2.
    """
    by_timepoint = {}
    for plate, tp, col_start, col_end in TIMEPOINT_BLOCKS:
        by_timepoint.setdefault(tp, []).extend(
            (plate, c) for c in range(col_start, col_end + 1)
        )
    col_maps = {}
    for tp, flat in by_timepoint.items():
        if len(flat) != 7:
            raise RuntimeError(
                f"Timepoint {tp} spans {len(flat)} well-columns, expected 7 "
                f"(design has 7 columns) -- check TIMEPOINT_BLOCKS."
            )
        col_maps[tp] = {i + 1: flat[i] for i in range(7)}
    return col_maps


def find_plate_anchor_row(ws, plate_num):
    """Locate the 'PLATE N' header row in column B of the left grid (robust to
    row-shift edits, rather than hardcoding row numbers)."""
    for row in ws.iter_rows(min_col=2, max_col=2):
        cell = row[0]
        if cell.value == f"PLATE {plate_num}":
            return cell.row
    raise RuntimeError(f"Could not locate 'PLATE {plate_num}' header in {SHEET_LAYOUT}")


def read_plate_grid(ws, plate_num):
    """Returns {(row_letter, well_col 1-12): cell_value} for one plate's left grid."""
    anchor_row = find_plate_anchor_row(ws, plate_num)
    data_start = anchor_row + 2  # anchor row = 'PLATE N', +1 = col-number row, +2 = row A
    grid = {}
    for i, row_letter in enumerate(ROWS):
        excel_row = data_start + i
        row_label = ws.cell(row=excel_row, column=1).value
        if row_label != row_letter:
            raise RuntimeError(
                f"Plate {plate_num}: expected row label {row_letter!r} at excel "
                f"row {excel_row}, found {row_label!r} -- sheet layout assumption broken."
            )
        for well_col in range(1, 13):
            val = ws.cell(row=excel_row, column=well_col + 1).value
            grid[(row_letter, well_col)] = val
    return grid


def build():
    design_df = load_design_table()
    col_maps = build_timepoint_col_map()

    wb = openpyxl.load_workbook(XLSX_PATH, data_only=True)
    ws = wb[SHEET_LAYOUT]

    plate_grids = {p: read_plate_grid(ws, p) for p in (1, 2, 3)}

    records = []

    # --- 1. Repeated-measures design wells x 4 timepoints ---
    for _, drow in design_df.iterrows():
        drow_letter, dcol = drow["design_row"], drow["design_col"]
        for tp, colmap in col_maps.items():
            plate, well_col = colmap[dcol]
            well = f"{drow_letter}{well_col:02d}"
            plate_label = f"plate {plate}"
            primer_info = IN_VITRO_PRIMER_TABLE.get(plate_label)
            if primer_info is None:
                raise KeyError(f"No in vitro primer entry for {plate_label!r}")
            r1, r2 = fastq_paths(
                primer_info["primer_pair"], primer_info["seq_run"], well,
                FASTQ_BASE_DIR_TEMPLATE,
            )
            records.append({
                "sample_uid": f"invitro_plate{plate}_{well}",
                "sequencing_plate": plate_label,
                "well": well,
                "primer_pair": primer_info["primer_pair"],
                "seq_run": primer_info["seq_run"],
                "sample_type": "experimental",
                "timepoint": tp,
                "design_well": f"{drow_letter}{dcol}",
                "SAMPLE_ID": drow["SAMPLE_ID"],
                "STRAIN": drow["STRAIN"],
                "TREATMENT": drow["TREATMENT"],
                "REPLICATE": drow["REPLICATE"],
                "SAMPLE_DESC": drow["SAMPLE_DESC"],
                "BIOREP_ID": drow["BIOREP_ID"],
                "fastq_R1": r1,
                "fastq_R2": r2,
            })

    # --- 2. Blank wells ---
    for plate, well_col in BLANK_COLUMNS:
        plate_label = f"plate {plate}"
        primer_info = IN_VITRO_PRIMER_TABLE[plate_label]
        for row_letter in ROWS:
            well = f"{row_letter}{well_col:02d}"
            cell_val = plate_grids[plate][(row_letter, well_col)]
            r1, r2 = fastq_paths(
                primer_info["primer_pair"], primer_info["seq_run"], well,
                FASTQ_BASE_DIR_TEMPLATE,
            )
            records.append({
                "sample_uid": f"invitro_plate{plate}_{well}",
                "sequencing_plate": plate_label,
                "well": well,
                "primer_pair": primer_info["primer_pair"],
                "seq_run": primer_info["seq_run"],
                "sample_type": "extraction_blank",
                "timepoint": pd.NA,
                "design_well": pd.NA,
                "SAMPLE_ID": pd.NA, "STRAIN": pd.NA, "TREATMENT": pd.NA,
                "REPLICATE": pd.NA, "SAMPLE_DESC": cell_val, "BIOREP_ID": pd.NA,
                "fastq_R1": r1, "fastq_R2": r2,
            })

    # --- 3. Plate 3 monoculture / t0 special region (content read from sheet) ---
    plate3_grid = plate_grids[3]
    for (col_start, col_end), category in PLATE3_SPECIAL_COLS.items():
        primer_info = IN_VITRO_PRIMER_TABLE["plate 3"]
        for well_col in range(col_start, col_end + 1):
            for row_letter in ROWS:
                well = f"{row_letter}{well_col:02d}"
                cell_val = plate3_grid[(row_letter, well_col)]
                is_blank_cell = (cell_val == "BLANKS") or (cell_val is None)
                r1, r2 = fastq_paths(
                    primer_info["primer_pair"], primer_info["seq_run"], well,
                    FASTQ_BASE_DIR_TEMPLATE,
                )
                records.append({
                    "sample_uid": f"invitro_plate3_{well}",
                    "sequencing_plate": "plate 3",
                    "well": well,
                    "primer_pair": primer_info["primer_pair"],
                    "seq_run": primer_info["seq_run"],
                    "sample_type": "extraction_blank" if is_blank_cell else category,
                    "timepoint": pd.NA,
                    "design_well": pd.NA,
                    "SAMPLE_ID": pd.NA,
                    "STRAIN": cell_val if category == "monoculture" else pd.NA,
                    "TREATMENT": pd.NA, "REPLICATE": pd.NA,
                    "SAMPLE_DESC": str(cell_val), "BIOREP_ID": pd.NA,
                    "fastq_R1": r1, "fastq_R2": r2,
                })

    df = pd.DataFrame(records)

    # --- Hard validation ---
    dupes = df[df.duplicated(subset=["sequencing_plate", "well"], keep=False)]
    if len(dupes) > 0:
        raise RuntimeError(
            f"{len(dupes)} rows have duplicate (sequencing_plate, well) -- "
            f"well-range boundaries overlap somewhere:\n"
            f"{dupes[['sample_uid','sequencing_plate','well','sample_type']].to_string()}"
        )

    n_biorep_missing = df[df["sample_type"] == "experimental"]["BIOREP_ID"].isna().sum()
    if n_biorep_missing > 0:
        raise RuntimeError(
            f"{n_biorep_missing} experimental rows missing BIOREP_ID -- repeated-measures "
            f"blocking for PERMANOVA will break."
        )

    return df


if __name__ == "__main__":
    df = build()
    print(f"Total rows: {len(df)}")
    print("\nBy sample_type:")
    print(df["sample_type"].value_counts(dropna=False))
    print("\nExperimental rows by timepoint:")
    print(df[df["sample_type"] == "experimental"]["timepoint"].value_counts())
    print("\nBIOREP_ID count check (should be 56 unique design wells x 4 timepoints = 224):")
    exp = df[df["sample_type"] == "experimental"]
    print(f"  experimental rows: {len(exp)}")
    print(f"  unique BIOREP_ID values: {exp['BIOREP_ID'].nunique()}")
    print(f"  unique design_well values: {exp['design_well'].nunique()}")

    out_path = "./in_vitro_master_metadata.csv"
    df.to_csv(out_path, index=False)
    print(f"\nWrote {out_path}")
    print(df.head(10).to_string())