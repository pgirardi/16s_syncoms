"""
Primer pair / sequencing-run lookup, keyed strictly off PHYSICAL WELL POSITION
on the collection/sequencing plate -- never off biological trial.

This is the fix for the original 00_build_master_metadata.R bug: primer_run()
used to key off `trial`, which silently misrouted 66/96 CollPlate4 wells.
"""

import re

# (Plate, Primer Pair, Sequencing Run) -- transcribed directly from Patricia,
# cross-checked against growth_coll_plate_mapping.xlsx 'sequencing' sheet.
# NOTE: plate 4b primer pair is PROVISIONAL (341F2_799F2) pending confirmation
# from read output -- flagged downstream in the QC report.
IN_PLANTA_PRIMER_TABLE = {
    "plate 1":  {"primer_pair": "341F3_799F2", "seq_run": 2},
    "plate 2":  {"primer_pair": "341F6_799F2", "seq_run": 2},
    "plate 3":  {"primer_pair": "341F3_799F1", "seq_run": 3},
    "plate 4a": {"primer_pair": "341F2_799F2", "seq_run": 3},
    "plate 4b": {"primer_pair": "341F2_799F2", "seq_run": 3},  # PROVISIONAL
    "plate 5":  {"primer_pair": "341F2_799F1", "seq_run": 3},
}

IN_VITRO_PRIMER_TABLE = {
    "plate 1": {"primer_pair": "341F2_799F1", "seq_run": 2},
    "plate 2": {"primer_pair": "341F3_799F1", "seq_run": 2},
    "plate 3": {"primer_pair": "341F2_799F2", "seq_run": 2},
}

_WELL_RE = re.compile(r"^([A-H])(\d{1,2})$")


def parse_well(well: str):
    """'D03' -> ('D', 3). Raises if the well string doesn't match A01-H12 format."""
    m = _WELL_RE.match(well.strip().upper())
    if not m:
        raise ValueError(f"Could not parse well position: {well!r}")
    row, col = m.group(1), int(m.group(2))
    if not (1 <= col <= 12):
        raise ValueError(f"Well column out of range 1-12: {well!r}")
    return row, col


def resolve_in_planta_plate_label(collection_plate: str, collection_well: str) -> str:
    """
    Maps a (collection_plate, collection_well) to the plate label used in
    IN_PLANTA_PRIMER_TABLE, handling the CollPlate4 -> plate 4a/4b split
    at well D03 (4a = A01-D02, 4b = D03 onward). This split is WELL-POSITION
    based, per Patricia's explicit confirmation, not trial-based.
    """
    m = re.match(r"^CollPlate(\d+)$", collection_plate.strip())
    if not m:
        raise ValueError(f"Unrecognized collection_plate format: {collection_plate!r}")
    plate_num = int(m.group(1))

    if plate_num != 4:
        return f"plate {plate_num}"

    row, col = parse_well(collection_well)
    # 4a = rows A-C (any col), or row D cols 1-2. 4b = row D cols 3-12, or rows E-H.
    if row in ("A", "B", "C") or (row == "D" and col <= 2):
        return "plate 4a"
    else:
        return "plate 4b"


def lookup_in_planta_primer(collection_plate: str, collection_well: str):
    """Returns dict with primer_pair, seq_run, plate_label. Raises on any lookup failure --
    this is the hard-stop replacement for the missing stopifnot() in the old R pipeline."""
    plate_label = resolve_in_planta_plate_label(collection_plate, collection_well)
    entry = IN_PLANTA_PRIMER_TABLE.get(plate_label)
    if entry is None:
        raise KeyError(
            f"No primer table entry for plate_label={plate_label!r} "
            f"(from collection_plate={collection_plate!r}, well={collection_well!r})"
        )
    return {"plate_label": plate_label, **entry}


def fastq_paths(primer_pair: str, seq_run: int, well: str, base_dir_template: str):
    """
    base_dir_template e.g. '/uufs/chpc.utah.edu/.../demulti_data_trial{seq_run}'
    Filename pattern confirmed from CHPC listing: {primer_pair}_R{1|2}_{well}.fastq
    e.g. 341F2_799F1_R1_A01.fastq
    """
    base_dir = base_dir_template.format(seq_run=seq_run)
    r1 = f"{base_dir}/{primer_pair}_R1_{well}.fastq"
    r2 = f"{base_dir}/{primer_pair}_R2_{well}.fastq"
    return r1, r2