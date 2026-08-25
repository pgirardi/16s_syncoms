#!/usr/bin/env Rscript
# =============================================================================
# 09b_permanova_inplanta.R
#
# PERMANOVA + betadisper + Aitchison PCA for the IN PLANTA trials.
#
# Ported from 05b_permanova_inplanta.R to consume the current pipeline's CLR
# output (07_clr_transform.py's long-format CSVs) instead of the old
# clr_wide_<method>_....rds files. Statistical design is UNCHANGED from 05b:
# each seedling is destructively sampled once, so ordinary free-permutation
# adonis2() is correct here (no split-plot blocking needed, unlike 09a).
# See "WHAT CHANGED FROM 05b" below for the input-format specifics.
#
# Runs BOTH in planta trials by default (pass one long-format --clr_file at a
# time; call the script twice, once per trial's file). Per project
# convention, trial 2 is the primary trial for the manuscript and trial 1
# goes in the supplement -- the script itself doesn't enforce this, just
# label your output directories accordingly when you call it.
#
# ---------------------------------------------------------------------------
# GROWTH_PLATE AND SEQ_PLATE AS COVARIATES (unchanged reasoning from 05b)
# ---------------------------------------------------------------------------
# growth_plate: expected to be crossed with dpi in both trials (every growth
# plate has samples at both dpi3 and dpi5) -- checked at runtime, not
# assumed; if it comes back fully aliased for your data, that CONTRADICTS
# the expected design and is worth investigating rather than just reporting
# as a limitation.
# seq_plate (== collection_plate): expected NOT clean -- trial 1's
# CollPlate1/CollPlate2 were perfectly aliased with dpi in the old pipeline
# (2 plates, 2 timepoints, one-to-one); trial 2 had CollPlate4 as the only
# plate with real within-dpi crossing. Re-verify with
# table(meta$collection_plate, meta$dpi) per trial if plate assignments
# change -- checked at runtime below, not hardcoded.
#
# ---------------------------------------------------------------------------
# WHY BUFFER / t0 / MONOCULTURE WELLS AREN'T HERE BY DEFAULT
# ---------------------------------------------------------------------------
# 07_clr_transform.py's in-planta CLR pass restricts to
# sample_type == 'experimental' AND truly_inoculated == TRUE -- t0, dropped,
# and monoculture wells never appear in the CLR long CSVs this script reads,
# so there's no equivalent of 05b's --include_t0 / --exclude_dropped flags to
# carry over; that filtering already happened upstream in Python.
# Buffer IS present in the CLR file (community_background == 'Buffer') since
# it's a real experimental sample_type, but is EXCLUDED from the default
# --communities list below (matching 05b, which only ever analyzed
# SynCom1/SynCom1C -- Buffer wells have no inoculated community to test a
# compositional shift within). Pass --communities Buffer,SynCom1,SynCom1C
# explicitly if you want Buffer analyzed too, e.g. as a negative-control
# check that treatment doesn't spuriously shift the background/contaminant
# signal in uninoculated wells.
#
# ---------------------------------------------------------------------------
# WHAT CHANGED FROM 05b
# ---------------------------------------------------------------------------
# - Input: --clr_file points at ONE trial's long-format CSV (e.g.
#   clr_inplanta_trial2_standard_long.csv) instead of --clr_dir + --method +
#   looping both trials internally -- call this script once per trial. Use
#   the _standard_long.csv file for the full-basis (old Fig 3) run, or the
#   _p25C2_excluded_long.csv file for the indirect-effect (old Fig 4) run --
#   same basis-consistency reasoning as 09a, see 07_clr_transform.py's
#   docstring; do NOT use the dual-basis clr_inplanta_trial{N}_long.csv here.
# - Communities inferred from community_background values in the file
#   (default: SynCom1, SynCom1C; override with --communities).
# - No --metadata / --exclude_dropped / --include_t0: see above.
# - "collection_plate" is the seq_plate equivalent (same physical plate,
#   renamed for clarity coming out of 07_clr_transform.py) -- column is
#   called seq_plate in this script's `samp` frame either way.
# - stat_ellipse() REMOVED from the PCA plot per request -- points only.
#   betadisper()/permutest() are UNCHANGED and still run and reported.
# - NEW: --dpi_subset lets you rerun restricted to one dpi (e.g.
#   --dpi_subset 5) for a "which single dpi should the figure use"
#   comparison -- same free-permutation treatment-only model as the full run,
#   just on a smaller sample set, with its own betadisper and PCA. The
#   default (no --dpi_subset) still produces the full treatment*dpi analysis
#   across both dpi values -- both are useful, run both.
#
# Inputs : --clr_file <path to a clr_inplanta_trial{1,2}_standard_long.csv or
#           _p25C2_excluded_long.csv, from 07_clr_transform.py>
# Outputs: <out_dir>/permanova_planta_<community>[_dpi<N>].csv
#          <out_dir>/betadisper_planta_<community>[_dpi<N>]_by_<factor>.csv
#          <out_dir>/betadisper_summary_planta_<community>[_dpi<N>].csv  (F/p per factor, NEW)
#          <out_dir>/pca_scores_planta_<community>[_dpi<N>].csv
#          <out_dir>/pca_loadings_planta_<community>[_dpi<N>].csv
#          <out_dir>/pca_variance_explained_planta_<community>[_dpi<N>].csv
#          <out_dir>/plots/pca_planta_<community>[_dpi<N>].png
#          <out_dir>/plots/scree_planta_<community>[_dpi<N>].png
#
# Usage:
#   Rscript 09b_permanova_inplanta.R --clr_file clr_out/clr_inplanta_trial2_standard_long.csv --out_dir permanova_planta_trial2
#   Rscript 09b_permanova_inplanta.R --clr_file clr_out/clr_inplanta_trial2_standard_long.csv --out_dir permanova_planta_trial2 --dpi_subset 5
#   Rscript 09b_permanova_inplanta.R --clr_file clr_out/clr_inplanta_trial1_standard_long.csv --out_dir permanova_planta_trial1
# =============================================================================

suppressPackageStartupMessages({
  library(vegan); library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(optparse)
})

opt_list <- list(
  make_option("--clr_file",  type="character",
              help="Path to a clr_inplanta_trial{1,2}_standard_long.csv or _p25C2_excluded_long.csv (from 07_clr_transform.py)"),
  make_option("--out_dir",   type="character", default="permanova_planta"),
  make_option("--n_perm",    type="integer", default=999),
  make_option("--dpi_subset", type="double", default=NA_real_,
              help="If set (e.g. 5), restrict to that single dpi. Default: NA (full analysis across all dpi)."),
  make_option("--communities", type="character", default="SynCom1,SynCom1C",
              help="Comma-separated community_background values to analyze. Default: 'SynCom1,SynCom1C' (excludes Buffer, matching 05b -- see header)."),
  make_option("--seed", type="integer", default=1)
)
opt <- parse_args(OptionParser(option_list=opt_list))
if (is.null(opt$clr_file)) stop("--clr_file is required (see header for which file to point it at).")
opt$out_dir <- sub("/+$", "", opt$out_dir)
plot_dir <- file.path(opt$out_dir, "plots")
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)
set.seed(opt$seed)

check_plate_identifiability <- function(samp, plate_col, group_col) {
  tab <- table(samp[[plate_col]], samp[[group_col]])
  informative <- rowSums(tab > 0) > 1
  list(tab = tab, fully_aliased = !any(informative),
       n_informative_plates = sum(informative), n_plates = nrow(tab))
}

clr_long <- read.csv(opt$clr_file, stringsAsFactors=FALSE)
required_cols <- c("sample_uid","taxon","clr_value","community_background","treatment","dpi")
missing_cols <- setdiff(required_cols, names(clr_long))
if (length(missing_cols)) stop("--clr_file is missing required column(s): ", paste(missing_cols, collapse=", "))

communities <- strsplit(opt$communities, ",")[[1]]
cat("Communities to analyze:", paste(communities, collapse=", "), "\n")
if (!is.na(opt$dpi_subset))
  cat("NOTE: --dpi_subset =", opt$dpi_subset, "-- SINGLE-DPI mode.\n")
cat("\n")

for (community in communities) {
  suffix <- if (!is.na(opt$dpi_subset)) paste0("_dpi", opt$dpi_subset) else ""
  cat("=== in planta |", community, if (!is.na(opt$dpi_subset)) paste0("| dpi=", opt$dpi_subset) else "", "===\n")

  sub_long <- clr_long %>% filter(community_background == community)
  if (!is.na(opt$dpi_subset)) sub_long <- sub_long %>% filter(dpi == opt$dpi_subset)
  if (nrow(sub_long) == 0) { cat("  SKIP -- no rows for this community/dpi\n\n"); next }

  clr_wide <- sub_long %>%
    select(sample_uid, taxon, clr_value) %>%
    pivot_wider(names_from = taxon, values_from = clr_value) %>%
    as.data.frame()
  rownames(clr_wide) <- clr_wide$sample_uid
  clr_mat <- as.matrix(clr_wide[, setdiff(colnames(clr_wide), "sample_uid"), drop=FALSE])

  samp <- sub_long %>%
    distinct(sample_uid, treatment, dpi, collection_well, growth_plate, seq_plate) %>%
    rename(timepoint = dpi)  # kept as `timepoint` internally to match 05b's variable names/plot code
  samp <- samp[match(rownames(clr_mat), samp$sample_uid), ]
  stopifnot(identical(rownames(clr_mat), samp$sample_uid))

  samp$treatment <- factor(samp$treatment)
  samp$timepoint <- factor(samp$timepoint)
  samp$seq_plate <- factor(samp$seq_plate)
  samp$growth_plate <- factor(samp$growth_plate)

  if (nrow(clr_mat) < 4 || nlevels(samp$treatment) < 2) {
    cat("  SKIP -- insufficient wells/factor levels after filtering\n\n"); next
  }
  d <- vegdist(clr_mat, method="euclidean")

  if (!is.na(opt$dpi_subset)) {
    mod_treatment <- adonis2(d ~ treatment, data=samp, permutations=opt$n_perm, by="margin")
    cat("  [treatment, free permutation -- single dpi]\n"); print(mod_treatment)

    combined <- bind_rows(
      as.data.frame(mod_treatment) %>% rownames_to_column("term") %>% filter(term == "treatment"),
      as.data.frame(mod_treatment) %>% rownames_to_column("term") %>% filter(term == "Residual") %>% mutate(term="Residual (treatment model)")
    )
    disp_facs <- c("treatment")
  } else {
    if (nlevels(samp$timepoint) < 2) {
      cat("  SKIP -- fewer than 2 dpi levels present\n\n"); next
    }
    mod_treatment <- adonis2(d ~ treatment, data=samp, permutations=opt$n_perm, by="margin")
    mod_timepoint <- adonis2(d ~ timepoint, data=samp, permutations=opt$n_perm, by="margin")
    mod <- adonis2(d ~ treatment * timepoint, data=samp, permutations=opt$n_perm, by="margin")
    cat("  [treatment]\n"); print(mod_treatment)
    cat("  [dpi]\n"); print(mod_timepoint)
    cat("  [interaction]\n"); print(mod)

    mod_growthplate <- NULL
    if (nlevels(samp$growth_plate) < 2) {
      cat("  growth_plate: SKIPPED -- fewer than 2 plates present\n")
    } else {
      gp_check <- check_plate_identifiability(samp, "growth_plate", "timepoint")
      if (gp_check$fully_aliased) {
        cat("  growth_plate: SKIPPED -- fully aliased with dpi for this subset (this contradicts\n",
            "  the expected crossed design -- worth investigating, see header).\n")
      } else {
        mod_growthplate <- adonis2(d ~ growth_plate, data=samp, permutations=opt$n_perm, by="margin")
        cat("  [growth_plate]\n"); print(mod_growthplate)
      }
    }

    mod_seqplate <- NULL
    mod_seqplate_sub <- NULL
    if (nlevels(samp$seq_plate) < 2) {
      cat("  seq_plate: SKIPPED -- only one plate present\n")
    } else {
      sp_check <- check_plate_identifiability(samp, "seq_plate", "timepoint")
      cat("  [seq_plate x dpi crosstab]\n"); print(sp_check$tab)
      if (sp_check$fully_aliased) {
        cat("  seq_plate: SKIPPED -- fully aliased with dpi; not estimable (expected for trial 1).\n")
      } else {
        if (sp_check$n_informative_plates == 1)
          cat("  seq_plate: fitting, but independent plate signal comes almost entirely from one\n",
              "  plate's internal split -- see sensitivity subset below.\n")
        mod_seqplate <- adonis2(d ~ seq_plate, data=samp, permutations=opt$n_perm, by="margin")
        cat("  [seq_plate]\n"); print(mod_seqplate)

        tp_by_plate <- sp_check$tab
        informative_plates <- rownames(tp_by_plate)[rowSums(tp_by_plate > 0) > 1]
        if (length(informative_plates) >= 1) {
          plate_sub <- informative_plates[1]
          samp_sub <- samp[samp$seq_plate == plate_sub, ]
          samp_sub$timepoint <- droplevels(samp_sub$timepoint)
          if (nlevels(samp_sub$timepoint) >= 2 && nrow(samp_sub) >= 4) {
            d_sub <- vegdist(clr_mat[rownames(clr_mat) %in% samp_sub$sample_uid, , drop=FALSE], method="euclidean")
            mod_seqplate_sub <- adonis2(d_sub ~ timepoint, data=samp_sub, permutations=opt$n_perm, by="margin")
            cat("  [SENSITIVITY: within seq_plate==", plate_sub, " only]\n", sep="")
            print(mod_seqplate_sub)
          }
        }
      }
    }

    saveRDS(list(treatment=mod_treatment, timepoint=mod_timepoint, interaction=mod,
                 growth_plate=mod_growthplate, seq_plate=mod_seqplate, seq_plate_sensitivity=mod_seqplate_sub),
            file.path(opt$out_dir, paste0("permanova_model_planta_", community, ".rds")))

    combined <- bind_rows(
      as.data.frame(mod_treatment) %>% rownames_to_column("term") %>% filter(term == "treatment"),
      as.data.frame(mod_timepoint) %>% rownames_to_column("term") %>% filter(term == "timepoint"),
      as.data.frame(mod) %>% rownames_to_column("term") %>% filter(term == "treatment:timepoint"),
      as.data.frame(mod_treatment) %>% rownames_to_column("term") %>% filter(term == "Residual") %>% mutate(term="Residual (treatment model)"),
      as.data.frame(mod_timepoint) %>% rownames_to_column("term") %>% filter(term == "Residual") %>% mutate(term="Residual (timepoint model)"),
      as.data.frame(mod) %>% rownames_to_column("term") %>% filter(term == "Residual") %>% mutate(term="Residual (interaction model)")
    )
    if (!is.null(mod_growthplate)) {
      combined <- bind_rows(combined,
        as.data.frame(mod_growthplate) %>% rownames_to_column("term") %>% filter(term == "growth_plate"),
        as.data.frame(mod_growthplate) %>% rownames_to_column("term") %>% filter(term == "Residual") %>% mutate(term="Residual (growth_plate model)"))
    }
    if (!is.null(mod_seqplate)) {
      combined <- bind_rows(combined,
        as.data.frame(mod_seqplate) %>% rownames_to_column("term") %>% filter(term == "seq_plate"),
        as.data.frame(mod_seqplate) %>% rownames_to_column("term") %>% filter(term == "Residual") %>% mutate(term="Residual (seq_plate model)"))
    }
    disp_facs <- c("treatment", "timepoint")
    if (!is.null(mod_growthplate)) disp_facs <- c(disp_facs, "growth_plate")
    if (!is.null(mod_seqplate)) disp_facs <- c(disp_facs, "seq_plate")
  }
  write.csv(combined, file.path(opt$out_dir, paste0("permanova_planta_", community, suffix, ".csv")), row.names=FALSE)

  # ── betadisper: location vs. spread. Kept exactly as in 05b -- a
  # statistical test, unaffected by the ellipse being removed from the plot
  # below. Ordinary free permutation throughout (destructive design, every
  # row independent, no stratum restriction needed here unlike 09a).
  # NOTE: F/p from permutest() now also collected into a
  # betadisper_summary_*.csv (previously console-only, see 09a for why). ────
  betadisper_summary_rows <- list()
  for (fac in disp_facs) {
    bd <- betadisper(d, samp[[fac]], type="centroid")
    pt <- permutest(bd, permutations=opt$n_perm)
    bd_out <- tibble(sample_uid=samp$sample_uid, group=samp[[fac]], dist_to_centroid=bd$distances)
    write.csv(bd_out, file.path(opt$out_dir, paste0("betadisper_planta_", community, suffix, "_by_", fac, ".csv")), row.names=FALSE)
    cat("  betadisper by", fac, "-- F =", round(pt$tab$F[1], 3), ", p =", round(pt$tab$`Pr(>F)`[1], 4), "\n")

    group_means <- bd_out %>% group_by(group) %>% summarise(mean_dist=mean(dist_to_centroid), .groups="drop")
    betadisper_summary_rows[[fac]] <- tibble(
      factor = fac,
      df_group = pt$tab$Df[1], df_residual = pt$tab$Df[2],
      F = pt$tab$F[1], n_perm = pt$tab$`N.Perm`[1], p_value = pt$tab$`Pr(>F)`[1],
      min_group_mean_dist = min(group_means$mean_dist), max_group_mean_dist = max(group_means$mean_dist),
      group_with_min_dist = group_means$group[which.min(group_means$mean_dist)],
      group_with_max_dist = group_means$group[which.max(group_means$mean_dist)]
    )
  }
  betadisper_summary <- bind_rows(betadisper_summary_rows)
  write.csv(betadisper_summary, file.path(opt$out_dir, paste0("betadisper_summary_planta_", community, suffix, ".csv")), row.names=FALSE)
  cat("  wrote", file.path(opt$out_dir, paste0("betadisper_summary_planta_", community, suffix, ".csv")), "\n")

  # ── Aitchison PCA ──────────────────────────────────────────────────────
  pca <- prcomp(clr_mat, center=TRUE, scale.=FALSE)
  var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  scores <- as.data.frame(pca$x) %>% rownames_to_column("sample_uid") %>% left_join(samp, by="sample_uid")
  write.csv(scores, file.path(opt$out_dir, paste0("pca_scores_planta_", community, suffix, ".csv")), row.names=FALSE)
  write.csv(as.data.frame(pca$rotation), file.path(opt$out_dir, paste0("pca_loadings_planta_", community, suffix, ".csv")))

  # NOTE: stat_ellipse() layer REMOVED per request -- points only.
  p <- ggplot(scores, aes(x=PC1, y=PC2, colour=treatment, shape=timepoint)) +
    geom_point(size=2.2, alpha=0.85) +
    theme_minimal(base_size=10) +
    labs(title=paste0("Aitchison PCA -- in planta | ", community, if (!is.na(opt$dpi_subset)) paste0(" | dpi ", opt$dpi_subset) else ""),
         subtitle=paste0("PC1 ", round(var_exp[1],1), "%, PC2 ", round(var_exp[2],1), "% of variance."),
         x=paste0("PC1 (", round(var_exp[1],1), "%)"),
         y=paste0("PC2 (", round(var_exp[2],1), "%)"))
  f_png <- file.path(plot_dir, paste0("pca_planta_", community, suffix, ".png"))
  ggsave(f_png, p, width=8, height=6)
  cat("  wrote", f_png, "\n")

  scree_df <- tibble(PC = factor(paste0("PC", seq_along(var_exp)), levels=paste0("PC", seq_along(var_exp))),
                      pct_variance = var_exp, cum_pct_variance = cumsum(var_exp))
  write.csv(scree_df, file.path(opt$out_dir, paste0("pca_variance_explained_planta_", community, suffix, ".csv")), row.names=FALSE)

  p_scree <- ggplot(scree_df, aes(x=PC)) +
    geom_col(aes(y=pct_variance), fill="grey70") +
    geom_line(aes(y=cum_pct_variance, group=1), colour="firebrick", linewidth=0.7) +
    geom_point(aes(y=cum_pct_variance), colour="firebrick", size=1.6) +
    scale_y_continuous(limits=c(0, 100), expand=expansion(mult=c(0, 0.05))) +
    theme_minimal(base_size=10) +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7)) +
    labs(title=paste0("Scree plot -- in planta | ", community, if (!is.na(opt$dpi_subset)) paste0(" | dpi ", opt$dpi_subset) else ""),
         x="Principal component", y="% variance explained")
  f_scree <- file.path(plot_dir, paste0("scree_planta_", community, suffix, ".png"))
  ggsave(f_scree, p_scree, width=8, height=5)
  cat("  wrote", f_scree, "\n\n")
}

cat("Done.\n")
