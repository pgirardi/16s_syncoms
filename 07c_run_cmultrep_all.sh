#!/bin/bash
# Run 07c_clr_cmultrepl.R over every 07 CLR output (merged + single-basis, all
# three trials). Output files keep the same names in a parallel directory, so
# 08 / 09a / 09b can be pointed at clr_out_cmult/ instead of clr_out/ unchanged.
#
# Usage:  bash run_07c_all.sh [clr_out] [clr_out_cmult] [community_strain_abundance.csv]
set -euo pipefail

IN_DIR=${1:-clr_out}
OUT_DIR=${2:-clr_out_cmult}
COMMUNITY=${3:-community_strain_abundance.csv}

module load R 2>/dev/null || true
mkdir -p "$OUT_DIR"

shopt -s nullglob
files=("$IN_DIR"/clr_*_long.csv)
if [ ${#files[@]} -eq 0 ]; then
  echo "No clr_*_long.csv files found in $IN_DIR" >&2; exit 1
fi

for f in "${files[@]}"; do
  echo "=== $(basename "$f")"
  Rscript 07b_clr_cmultrep.R \
    --template "$f" \
    --community "$COMMUNITY" \
    --output "$OUT_DIR/$(basename "$f")"
  echo
done
echo "Done: ${#files[@]} files -> $OUT_DIR/"
