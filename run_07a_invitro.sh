#!/bin/bash
# Run the three in vitro PERMANOVA analyses, dropping low-confidence taxa.
#
#   bash run_07a_invitro.sh                  # uses clr_out/
#   bash run_07a_invitro.sh clr_out_cmult    # same runs on the cmultRepl CLR files
#
# Outputs carry a _dropLC suffix, so they won't overwrite runs made without the flag.
# To run WITHOUT dropping, set DROP_FLAG="" below.
set -euo pipefail

CLR_DIR=${1:-clr_out}
DROP_FLAG="--drop_low_confidence_taxa"

module load R 2>/dev/null || true   # same R that has vegan installed

run() {
  echo "=================================================================="
  echo ">>> Rscript 07a_permanova_invitro.R $* $DROP_FLAG"
  echo "=================================================================="
  Rscript 07a_permanova_invitro.R "$@" $DROP_FLAG
  echo
}

# 1. full split-plot analysis, all timepoints (full CLR basis)
run --clr_file "$CLR_DIR/clr_invitro_standard_long.csv" --out_dir permanova_invitro

# 2. 24H only (free permutation)
run --clr_file "$CLR_DIR/clr_invitro_standard_long.csv" --out_dir permanova_invitro --timepoint_subset 24H

# 3. indirect-effect analysis, p25C2 excluded from the CLR basis
run --clr_file "$CLR_DIR/clr_invitro_p25C2_excluded_long.csv" --out_dir permanova_invitro_excl_p25C2

echo "All three runs finished."
