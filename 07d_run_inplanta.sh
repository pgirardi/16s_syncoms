#!/bin/bash
# Run the three in planta PERMANOVA analyses, dropping low-confidence taxa.
#
#   bash 07d_run_inplanta.sh                  # uses clr_out/
#   bash 07d_run_inplanta.sh clr_out_cmult    # same runs on the cmultRepl CLR files
#
# Outputs carry a _dropLC suffix, so they won't overwrite runs made without the flag.
# To run WITHOUT dropping, set DROP_FLAG="" below.
set -euo pipefail

CLR_DIR=${1:-clr_out}
DROP_FLAG="--drop_low_confidence_taxa"

module load R 2>/dev/null || true   # same R that has vegan installed

run() {
  echo "=================================================================="
  echo ">>> Rscript 07b_permanova_inplanta.R $* $DROP_FLAG"
  echo "=================================================================="
  Rscript 07b_permanova_inplanta.R "$@" $DROP_FLAG
  echo
}

# 1. trial 2, full treatment x dpi model
run --clr_file "$CLR_DIR/clr_inplanta_trial2_standard_long.csv" --out_dir permanova_planta_trial2

# 2. trial 2, dpi 5 only
run --clr_file "$CLR_DIR/clr_inplanta_trial2_standard_long.csv" --out_dir permanova_planta_trial2 --dpi_subset 5

# 3. trial 1, full treatment x dpi model
run --clr_file "$CLR_DIR/clr_inplanta_trial1_standard_long.csv" --out_dir permanova_planta_trial1

echo "All three runs finished."
