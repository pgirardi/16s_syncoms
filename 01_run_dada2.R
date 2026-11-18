#!/usr/bin/env Rscript
#
# 01_run_dada2.R
#
# Runs DADA2 denoising on R1 and R2 INDEPENDENTLY (own error model each, no
# read merging -- confirmed appropriate for the ~450bp amplicon in prior
# sessions: --rev_also, not --paired). Samples are grouped by
# (primer_pair, seq_run) -- a "demux group" -- since each such group was
# sequenced together and needs its own learned error model; groups are never
# pooled across each other.
#
# Usage:
#   Rscript 01_run_dada2.R <master_metadata.csv> <output_dir>
#
# Input CSV must have columns: sample_uid, fastq_R1, fastq_R2, primer_pair, seq_run


suppressMessages(library(dada2))
suppressMessages(library(dplyr))

n_threads <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
if (is.na(n_threads) || n_threads < 1) {
  n_threads <- suppressWarnings(as.integer(Sys.getenv("SLURM_NTASKS", "")))
}
if (is.na(n_threads) || n_threads < 1) {
  n_threads <- 1
  cat("WARNING: could not detect SLURM core allocation from environment variables ",
      "(SLURM_CPUS_PER_TASK / SLURM_NTASKS both unset or invalid). Defaulting to ",
      "multithread = 1. If you know your salloc allocation size, consider setting ",
      "n_threads manually at the top of this script.\n", sep = "")
}
cat(sprintf("Using multithread = %d (detected from SLURM allocation)\n", n_threads))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 01_run_dada2.R <master_metadata.csv> <output_dir>")
}
metadata_path <- args[1]
output_dir <- args[2]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

TRUNC_LEN <- list(
  # "341F3_799F2_run2" = list(R1 = 220, R2 = 200),
)

get_trunc_len <- function(demux_group, direction) {
  if (!is.null(TRUNC_LEN[[demux_group]]) && !is.null(TRUNC_LEN[[demux_group]][[direction]])) {
    return(TRUNC_LEN[[demux_group]][[direction]])
  }
  cat(sprintf(
    "    WARNING: no truncLen set for %s/%s -- using 0 (no truncation). Set this after inspecting quality profiles.\n",
    demux_group, direction
  ))
  return(0)
}

meta <- read.csv(metadata_path, stringsAsFactors = FALSE)

required_cols <- c("sample_uid", "fastq_R1", "fastq_R2", "primer_pair", "seq_run")
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("Metadata missing required columns: ", paste(missing_cols, collapse = ", "))
}

if ("sample_type" %in% colnames(meta)) {
  n_before <- nrow(meta)
  meta <- meta %>% filter(sample_type != "empty")
  cat(sprintf("Excluded %d rows with sample_type == 'empty' (%d remain)\n",
              n_before - nrow(meta), nrow(meta)))
}

meta$demux_group <- paste0(meta$primer_pair, "_run", meta$seq_run)

read_track_all <- list()

run_one_direction <- function(grp, grp_meta, grp_dir, direction) {
  cat(sprintf("  -- direction: %s\n", direction))

  out_rds <- file.path(grp_dir, sprintf("seqtab_nochim_%s.rds", direction))
  if (file.exists(out_rds)) {
    cat(sprintf("    %s already exists, skipping (delete it to force a rerun).\n", out_rds))
    return(NULL)
  }

  fastq_col <- paste0("fastq_", direction)
  raw_paths <- grp_meta[[fastq_col]]
  sample_uids <- grp_meta$sample_uid

  exists_mask <- file.exists(raw_paths)
  if (!all(exists_mask)) {
    missing <- raw_paths[!exists_mask]
    cat(sprintf("    WARNING: %d/%d fastq files not found, skipping those samples:\n",
                sum(!exists_mask), length(raw_paths)))
    cat(paste("     ", missing, collapse = "\n"), "\n")
  }
  raw_paths <- raw_paths[exists_mask]
  uids <- sample_uids[exists_mask]

  if (length(raw_paths) == 0) {
    cat("    No files found for this group/direction, skipping.\n")
    return(NULL)
  }

  filt_dir <- file.path(grp_dir, paste0("filtered_", direction))
  dir.create(filt_dir, showWarnings = FALSE, recursive = TRUE)
  filt_paths <- file.path(filt_dir, paste0(uids, "_filt.fastq.gz"))
  names(filt_paths) <- uids
  names(raw_paths) <- uids

  filt_out <- filterAndTrim(
    raw_paths, filt_paths,
    truncLen = get_trunc_len(grp, direction),
    maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE,
    compress = TRUE, multithread = n_threads
  )

  # filterAndTrim() names its output rows by basename(raw_paths), NOT by
  # names(raw_paths) -- confirmed the hard way (subscript-out-of-bounds when
  # indexing by sample_uid). Row order IS preserved relative to the input
  # vector though, so overwrite rownames explicitly rather than rely on
  # name-matching against a naming scheme filterAndTrim doesn't actually use.
  if (nrow(filt_out) != length(uids)) {
    stop(sprintf(
      "filterAndTrim() returned %d rows but expected %d (one per input file) -- row order assumption is broken, do not proceed blindly.",
      nrow(filt_out), length(uids)
    ))
  }
  rownames(filt_out) <- uids

  filt_exists <- file.exists(filt_paths)
  if (!all(filt_exists)) {
    cat(sprintf("    %d samples had 0 reads survive filtering, excluding from error learning:\n",
                sum(!filt_exists)))
    cat(paste("     ", uids[!filt_exists], collapse = "\n"), "\n")
  }
  filt_paths_ok <- filt_paths[filt_exists]

  if (length(filt_paths_ok) == 0) {
    cat("    No samples survived filtering for this group/direction, skipping.\n")
    return(NULL)
  }

  err <- learnErrors(filt_paths_ok, multithread = n_threads)
  dada_out <- dada(filt_paths_ok, err = err, multithread = n_threads)
  seqtab <- makeSequenceTable(dada_out)
  seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                       multithread = n_threads, verbose = TRUE)

  saveRDS(seqtab_nochim, out_rds)
  cat(sprintf("    Wrote %s (%d samples x %d ASVs)\n",
              out_rds, nrow(seqtab_nochim), ncol(seqtab_nochim)))

  getN <- function(x) sum(getUniques(x))
  dada_reads <- if (length(filt_paths_ok) == 1) {
    setNames(getN(dada_out), names(filt_paths_ok))
  } else {
    sapply(dada_out, getN)
  }

  track <- data.frame(
    sample_uid = uids,
    demux_group = grp,
    direction = direction,
    input = filt_out[uids, "reads.in"],
    filtered = filt_out[uids, "reads.out"],
    stringsAsFactors = FALSE
  )
  track$denoised <- NA
  track$denoised[match(names(dada_reads), track$sample_uid)] <- dada_reads
  track
}

for (grp in unique(meta$demux_group)) {
  cat(sprintf("\n=== Processing demux group: %s ===\n", grp))
  grp_meta <- meta %>% filter(demux_group == grp)
  grp_dir <- file.path(output_dir, grp)
  dir.create(grp_dir, showWarnings = FALSE, recursive = TRUE)

  for (direction in c("R1", "R2")) {
    track <- run_one_direction(grp, grp_meta, grp_dir, direction)
    if (!is.null(track)) {
      read_track_all[[paste(grp, direction)]] <- track
    }
  }
}

if (length(read_track_all) > 0) {
  read_track_df <- bind_rows(read_track_all)
  out_path <- file.path(output_dir, "read_tracking.csv")
  write.csv(read_track_df, out_path, row.names = FALSE)
  cat(sprintf("\nWrote read tracking table: %s\n", out_path))
} else {
  cat("\nNo samples were successfully processed -- read_tracking.csv not written.\n")
}

inspect_quality_profiles <- function(metadata_path, output_dir, n_samples = 6) {
  meta <- read.csv(metadata_path, stringsAsFactors = FALSE)
  if ("sample_type" %in% colnames(meta)) {
    meta <- meta %>% filter(sample_type != "empty")
  }
  meta$demux_group <- paste0(meta$primer_pair, "_run", meta$seq_run)
  qc_dir <- file.path(output_dir, "quality_profiles")
  dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

  for (grp in unique(meta$demux_group)) {
    grp_meta <- meta %>% filter(demux_group == grp)
    for (direction in c("R1", "R2")) {
      fastq_col <- paste0("fastq_", direction)
      paths <- grp_meta[[fastq_col]]
      paths <- paths[file.exists(paths)]
      if (length(paths) == 0) next
      sample_paths <- head(paths, n_samples)
      pdf(file.path(qc_dir, sprintf("%s_%s_quality.pdf", grp, direction)))
      print(plotQualityProfile(sample_paths))
      dev.off()
    }
  }
  cat(sprintf("Quality profile PDFs written to %s -- inspect these before setting TRUNC_LEN.\n", qc_dir))
}
