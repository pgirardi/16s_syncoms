#!/usr/bin/env Rscript
#
# export_seqtab_to_csv.R
#
# Format-bridge only, no analysis logic: converts a DADA2 seqtab_nochim .rds
# (samples x ASV-sequence count matrix, ASV sequences as column names) into a
# long-format CSV (sample_uid, sequence, count) that's easy to work with in
# Python/pandas. All actual reference-building logic lives in
# build_monoculture_reference_fasta.py -- this script does nothing else.
#
# Usage:
#   Rscript export_seqtab_to_csv.R <seqtab.rds> <output.csv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript export_seqtab_to_csv.R <seqtab.rds> <output.csv>")
}
rds_path <- args[1]
out_path <- args[2]

seqtab <- readRDS(rds_path)

long_df <- data.frame(
  sample_uid = rep(rownames(seqtab), times = ncol(seqtab)),
  sequence   = rep(colnames(seqtab), each = nrow(seqtab)),
  count      = as.vector(seqtab),
  stringsAsFactors = FALSE
)
long_df <- long_df[long_df$count > 0, ]  # drop zero entries, matrix is sparse

write.csv(long_df, out_path, row.names = FALSE)
cat(sprintf("Wrote %d rows (%d samples x %d ASVs, non-zero only) to %s\n",
            nrow(long_df), nrow(seqtab), ncol(seqtab), out_path))
