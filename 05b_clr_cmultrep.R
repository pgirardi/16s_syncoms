#!/usr/bin/env Rscript
# =============================================================================
# 05b_clr_cmultrepl.R
#
# Recomputes CLR values with Bayesian-multiplicative zero replacement
# (zCompositions::cmultRepl) instead of 05a_clr_transform.py's flat pseudocount
# (half the smallest non-zero raw count in the pass matrix), so the two can be
# compared head-to-head.
#
# DESIGN: template-driven, not a re-implementation of 05a.
#   Reads an EXISTING 05a output and keeps every row and every column exactly
#   as-is (samples, taxa, clr_basis, flags, metadata, rel_abundance_pct) --
#   only clr_value is swapped. The old value is kept as clr_value_pseudocount.
#   Output therefore drops straight into 06 / 07a / 07b with no code changes,
#   so the statistical design (incl. in vitro restricted permutation) is
#   identical between the two runs.
#
# CLR PASSES -- mirrors 05a exactly. Every pass runs over ALL samples in the
#   file (05a computes one matrix per trial, not per background):
#     * merged file (has clr_basis column, e.g. clr_invitro_long.csv):
#         - "standard" pass: parts = every taxon in the file, p25C2 included,
#           for every sample. Supplies rows with clr_basis 'standard' (e.g.
#           all of S1) and 'standard_p25C2_only' (p25C2 in S1+C2 / SynCom1C).
#         - "p25C2_excluded" pass: parts = every taxon EXCEPT p25C2, for every
#           sample. Supplies rows with clr_basis 'p25C2_excluded'.
#     * single-basis file (no clr_basis column, e.g.
#       clr_invitro_standard_long.csv or clr_invitro_p25C2_excluded_long.csv
#       -- the ones PERMANOVA should read): one pass, parts = every taxon in
#       the file.
#   Parts are always the taxa present in the template (05a already dropped
#   Unassigned, the ambiguous 59+p25C2 call, and never-observed taxa).
#
# BUILT-IN CHECK: before computing anything new, the script recomputes CLR
#   with 05a's OLD rule from the raw counts and compares it with the template's
#   clr_value. If they don't match, it stops -- a comparison on a different
#   basis would not be like-for-like.
#
# Zero replacement: cmultRepl(method = "GBM" by default, output = "p-counts",
#   z.delete = FALSE, z.warning = 1). z.delete = FALSE matters: cmultRepl's
#   default silently DROPS any part that is >80% zeros, which would remove
#   Leaf391 (~90% zeros) from the basis. If GBM can't be estimated (a part
#   non-zero in <2 samples), the pass falls back to CZM and says so; the
#   method actually used is recorded per row in `zero_replacement`.
#
# Usage (CHPC):
#   module load R
#   Rscript 05b_clr_cmultrepl.R \
#       --template clr_out/clr_invitro_long.csv \
#       --community community_strain_abundance.csv \
#       --output clr_out_cmult/clr_invitro_long.csv \
#       [--method GBM]
#   Run once per 05 output file used downstream (merged + single-basis,
#   for each of in vitro / in planta trial 1 / trial 2)
# =============================================================================

suppressPackageStartupMessages(library(zCompositions))

# ---- args -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) {
    if (is.null(default)) stop("Missing required argument ", flag)
    return(default)
  }
  args[i + 1]
}
template_path  <- get_arg("--template")
community_path <- get_arg("--community")
output_path    <- get_arg("--output")
method         <- get_arg("--method", "GBM")
match_tol      <- as.numeric(get_arg("--match-tol", "1e-6"))
P25            <- get_arg("--dual-clr-taxon", "p25C2")
if (!method %in% c("GBM", "SQ", "BL", "CZM")) stop("--method must be one of GBM, SQ, BL, CZM")

# ---- load -------------------------------------------------------------------
tmpl <- read.csv(template_path, stringsAsFactors = FALSE, check.names = FALSE)
comm <- read.csv(community_path, stringsAsFactors = FALSE,
                 colClasses = c(strain_or_group = "character"))
for (col in c("sample_uid", "taxon", "clr_value")) {
  if (!col %in% names(tmpl)) stop("Template is missing column '", col, "'")
}

# same relabel() as 05a_clr_transform.py
relabel <- function(s) {
  s <- trimws(as.character(s))
  ifelse(tolower(s) == "unassigned", "Unassigned",
  ifelse(s == "p25C2", "p25C2",
  ifelse(s == "59+p25C2", "Leaf59+p25C2 (ambiguous)",
         paste0("Leaf", sub("\\.0$", "", s)))))
}
comm$taxon <- relabel(comm$strain_or_group)

tmpl_taxa <- sort(unique(tmpl$taxon))
missing_taxa <- setdiff(tmpl_taxa, comm$taxon)
if (length(missing_taxa) > 0) {
  stop("Template taxa with no matching label in the community table: ",
       paste(missing_taxa, collapse = ", "))
}

samples <- unique(tmpl$sample_uid)
counts <- tapply(comm$R2_count, list(comm$sample_uid, comm$taxon), sum)
counts[is.na(counts)] <- 0
absent <- setdiff(samples, rownames(counts))
if (length(absent) > 0) {
  stop(length(absent), " template samples have no rows in the community table, e.g. ",
       paste(head(absent, 5), collapse = ", "))
}
counts <- counts[samples, , drop = FALSE]

# ---- define passes (see header) ----------------------------------------------
has_basis <- "clr_basis" %in% names(tmpl)
if (has_basis) {
  is_excl_row <- grepl("excluded", tmpl$clr_basis)
  row_pass <- ifelse(is_excl_row, "p25C2_excluded", "standard")
  passes <- list(
    standard       = tmpl_taxa,
    p25C2_excluded = setdiff(tmpl_taxa, P25)
  )
  passes <- passes[unique(row_pass)]
} else {
  row_pass <- rep("single", nrow(tmpl))
  passes <- list(single = tmpl_taxa)
}

cat(sprintf("Template: %s\n  %d samples, %d taxa, %s\n", template_path, length(samples),
            length(tmpl_taxa), if (has_basis) "merged dual-CLR file" else "single-basis file"))
for (nm in names(passes)) {
  cat(sprintf("  - pass %-15s n=%3d  rows kept=%4d  parts=%s\n", nm, length(samples),
              sum(row_pass == nm), paste(passes[[nm]], collapse = ",")))
}

# ---- zero-replacement rules -------------------------------------------------
clr <- function(m) { l <- log(m); sweep(l, 1, rowMeans(l)) }

old_rule <- function(x) {             # 05a: half the smallest non-zero raw count
  pc <- min(x[x > 0]) / 2
  x[x == 0] <- pc
  list(values = clr(x), method = "pseudocount_halfmin")
}

cmult_rule <- function(x) {
  if (all(x > 0)) return(list(values = clr(x), method = "none_needed"))
  allzero <- colnames(x)[colSums(x) == 0]
  if (length(allzero) > 0) {
    stop("Part(s) zero in every sample, cmultRepl cannot impute: ",
         paste(allzero, collapse = ", "))
  }
  used <- method
  n_pos <- colSums(x > 0)
  if (method == "GBM" && any(n_pos < 2)) {
    used <- "CZM"
    message(sprintf("  NOTE: GBM not estimable (non-zero in <2 samples: %s); using CZM for this pass",
                    paste(names(n_pos)[n_pos < 2], collapse = ", ")))
  }
  rep <- cmultRepl(x, label = 0, method = used, output = "p-counts",
                   z.delete = FALSE, z.warning = 1, suppress.print = TRUE)
  list(values = clr(as.matrix(rep)), method = used)
}

compute <- function(rule) {
  vals <- rep(NA_real_, nrow(tmpl)); meth <- rep(NA_character_, nrow(tmpl))
  for (nm in names(passes)) {
    res <- rule(counts[, passes[[nm]], drop = FALSE])
    idx <- which(row_pass == nm)
    vals[idx] <- res$values[cbind(tmpl$sample_uid[idx], tmpl$taxon[idx])]
    meth[idx] <- res$method
  }
  list(values = vals, method = meth)
}

# ---- 1. verify the reconstruction against 05a's own clr_value -----------------
old <- compute(old_rule)
err <- max(abs(old$values - tmpl$clr_value), na.rm = TRUE)
cat(sprintf("Reconstruction check (05a rule, half-min raw-count pseudocount): max |diff| = %.3g\n", err))
if (!is.finite(err) || err > match_tol) {
  worst <- order(-abs(old$values - tmpl$clr_value))[1:min(5, nrow(tmpl))]
  print(data.frame(tmpl[worst, intersect(c("sample_uid", "taxon", "clr_basis"), names(tmpl))],
                   template = tmpl$clr_value[worst], recomputed = old$values[worst]))
  stop("Could not reproduce 05a's clr_value (max |diff| = ", signif(err, 3),
       "). Comparison would not be like-for-like -- send this output so the ",
       "reconstruction can be fixed.")
}
cat("  -> reproduced 05a exactly.\n")

# ---- 2. cmultRepl CLR --------------------------------------------------------
new <- compute(cmult_rule)

out <- tmpl
out$clr_value_pseudocount <- tmpl$clr_value
out$clr_value <- new$values
out$zero_replacement <- paste0("cmultRepl_", new$method)
out$raw_count <- counts[cbind(tmpl$sample_uid, tmpl$taxon)]
out$is_zero <- out$raw_count == 0

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(out, output_path, row.names = FALSE)
cat("Wrote ", output_path, "\n", sep = "")

# ---- 3. diagnostics ----------------------------------------------------------
diag <- do.call(rbind, lapply(split(out, out$taxon), function(d) {
  delta <- d$clr_value - d$clr_value_pseudocount
  data.frame(
    taxon = d$taxon[1],
    frac_zero = round(mean(d$is_zero), 2),
    mean_shift_zeros = if (any(d$is_zero)) round(mean(delta[d$is_zero]), 2) else NA,
    mean_shift_nonzeros = round(mean(delta[!d$is_zero]), 2),
    cor_nonzeros = round(suppressWarnings(
      cor(d$clr_value[!d$is_zero], d$clr_value_pseudocount[!d$is_zero])), 3))
}))
cat("\nPer-taxon change (cmultRepl - pseudocount):\n")
print(diag, row.names = FALSE)
diag_path <- sub("\\.csv$", "_zero_handling_diagnostics.csv", output_path)
write.csv(diag, diag_path, row.names = FALSE)
cat("Wrote ", diag_path, "\n", sep = "")
