#!/usr/bin/env Rscript
# =============================================================================
# 07b_permanova_inplanta.R
#
# PERMANOVA + betadisper + Aitchison PCA for the IN PLANTA trials.
#
# each seedling is destructively sampled once, so ordinary free-permutation
# adonis2() is correct here (no split-plot blocking needed, unlike 07a).
#
# Runs BOTH in planta trials by default (pass one long-format --clr_file at a
# time; call the script twice, once per trial's file). Per project
# convention, trial 2 is the primary trial for the manuscript and trial 1
# goes in the supplement -- the script itself doesn't enforce this, just
# label your output directories accordingly when you call it.
#
# ---------------------------------------------------------------------------
# GROWTH_PLATE AND SEQ_PLATE AS COVARIATES
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
# Inputs : --clr_file <path to a clr_inplanta_trial{1,2}_standard_long.csv or
#           _p25C2_excluded_long.csv, from 05_clr_transform.py>
# Outputs: <out_dir>/permanova_planta_<community>[_dpi<N>].csv
#          <out_dir>/betadisper_planta_<community>[_dpi<N>]_by_<factor>.csv
#          <out_dir>/betadisper_summary_planta_<community>[_dpi<N>].csv  (F/p per factor, NEW)
#          <out_dir>/pca_scores_planta_<community>[_dpi<N>].csv
#          <out_dir>/pca_loadings_planta_<community>[_dpi<N>].csv
#          <out_dir>/pca_variance_explained_planta_<community>[_dpi<N>].csv
#          <out_dir>/plots/pca_planta_<community>[_dpi<N>].png
#          <out_dir>/plots/pca_planta_<community>[_dpi<N>].pdf   (NEW, vector)
#          <out_dir>/plots/scree_planta_<community>[_dpi<N>].png
#          <out_dir>/plots/scree_planta_<community>[_dpi<N>].pdf (NEW, vector)
#
# Usage:
#   Rscript 07b_permanova_inplanta.R --clr_file clr_out/clr_inplanta_trial2_standard_long.csv --out_dir permanova_planta_trial2
#   Rscript 07b_permanova_inplanta.R --clr_file clr_out/clr_inplanta_trial2_standard_long.csv --out_dir permanova_planta_trial2 --dpi_subset 5
#   Rscript 07b_permanova_inplanta.R --clr_file clr_out/clr_inplanta_trial1_standard_long.csv --out_dir permanova_planta_trial1
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
  make_option("--drop_low_confidence_taxa", action="store_true", default=FALSE,
              help="Drop taxa flagged low_confidence_flag (from 07) and re-center the CLR on the remaining taxa before PERMANOVA/betadisper/PCA. Sensitivity analysis; output files get a _dropLC suffix. Default: keep all taxa."),
  make_option("--seed", type="integer", default=1)
)
opt <- parse_args(OptionParser(option_list=opt_list))
if (is.null(opt$clr_file)) stop("--clr_file is required (see header for which file to point it at).")
opt$out_dir <- sub("/+$", "", opt$out_dir)
plot_dir <- file.path(opt$out_dir, "plots")
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)
set.seed(opt$seed)

# Shared treatment palette/label standard -- IDENTICAL in 07a and 07b so
# Tailocin/delTFA/Buffer render as the same color regardless of which system
# or script produced the plot. Colors from NatParksPalettes "DeathValley"
# (the project palette; flagged colorblind-safe in the package), chosen to
# differ in both hue and lightness: slate blue / gold / dark rust.
# drop=FALSE in scale_colour_manual (used below) keeps the mapping fixed even
# if a given subset is missing one treatment level.
TREATMENT_COLORS <- c(Buffer="#435F90", delTFA="#FEB359", Tailocin="#8C2B0E")
# Scree plot colors, same palette: bars (individual PC) + line (cumulative).
SCREE_BAR_COLOR  <- "#B47E83"
SCREE_LINE_COLOR <- "#132F5B"

FONT <- list(base = 14, axis_text = 12, scree_pc_labels = 10, title = 16, subtitle = 13)
PLOT_THEME <- theme_minimal(base_size = FONT$base) +
  theme(axis.text     = element_text(size = FONT$axis_text),
        legend.text   = element_text(size = FONT$axis_text),
        legend.title  = element_text(size = FONT$base),
        plot.title    = element_text(size = FONT$title, face = "bold"),
        plot.subtitle = element_text(size = FONT$subtitle))
TREATMENT_LEVELS <- names(TREATMENT_COLORS)
# in planta uses shorthand treatment codes (Buf/TFA/Tail) -- standardized to
# the same full names used in vitro so labels/colors line up across systems.
TREATMENT_RELABEL <- c(Buf="Buffer", TFA="delTFA", Tail="Tailocin")

check_plate_identifiability <- function(samp, plate_col, group_col) {
  tab <- table(samp[[plate_col]], samp[[group_col]])
  informative <- rowSums(tab > 0) > 1
  list(tab = tab, fully_aliased = !any(informative),
       n_informative_plates = sum(informative), n_plates = nrow(tab))
}

# ── optional: drop low-confidence taxa (--drop_low_confidence_taxa) ─────────
# Dropping columns from an already-computed CLR matrix is NOT enough on its
# own: every CLR value was centered on the geometric mean of ALL taxa,
# including the ones being dropped, so the leftover columns are not a valid
# CLR of the reduced community and their Euclidean distances would be wrong.
# Re-centering each row on the mean of the RETAINED columns gives exactly the
# CLR of the subcomposition (clr_S(x)_i = clr(x)_i - mean_{j in S} clr(x)_j),
# i.e. the same result as recomputing CLR from scratch on only the retained
# taxa, using the same zero-replaced values.
drop_low_conf_taxa <- function(clr_mat, flagged) {
  keep <- setdiff(colnames(clr_mat), flagged)
  if (length(keep) < 2) stop("--drop_low_confidence_taxa would leave fewer than 2 taxa")
  sub <- clr_mat[, keep, drop = FALSE]
  sweep(sub, 1, rowMeans(sub))
}
P25C2_INOCULATED_BACKGROUNDS <- c("S1+C2", "SynCom1C")
SYSTEM_TAG <- "planta"

clr_long <- read.csv(opt$clr_file, stringsAsFactors=FALSE)

as_flag <- function(x) toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1")
if ("low_confidence_flag" %in% names(clr_long)) {
  clr_long$low_confidence_flag <- as_flag(clr_long$low_confidence_flag)
} else {
  if (opt$drop_low_confidence_taxa) stop("--drop_low_confidence_taxa set, but --clr_file has no low_confidence_flag column")
  clr_long$low_confidence_flag <- FALSE
}
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
  if (opt$drop_low_confidence_taxa) suffix <- paste0(suffix, "_dropLC")
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

  if (!(community %in% P25C2_INOCULATED_BACKGROUNDS) && "p25C2" %in% colnames(clr_mat)) {
    clr_mat <- drop_low_conf_taxa(clr_mat, "p25C2")
    cat("  p25C2 removed from parts (never inoculated into", community, ") -- CLR re-centered on",
        ncol(clr_mat), "taxa\n")
  }

  # Informational by default; --drop_low_confidence_taxa removes these and
  # re-centers the CLR on the remaining taxa (see drop_low_conf_taxa()).
  taxon_confidence <- sub_long %>% distinct(taxon, low_confidence_flag)
  flagged_taxa <- taxon_confidence %>% filter(low_confidence_flag == TRUE) %>% pull(taxon)
  flagged_taxa <- intersect(flagged_taxa, colnames(clr_mat))
  taxon_confidence <- taxon_confidence %>% filter(taxon %in% colnames(clr_mat))
  if (length(flagged_taxa) > 0) {
    cat("  NOTE: taxa flagged low_confidence in this basis (blank-derived QC, from 07):",
        paste(flagged_taxa, collapse=", "), "\n")
  }
  if (opt$drop_low_confidence_taxa && length(flagged_taxa) > 0) {
    clr_mat <- drop_low_conf_taxa(clr_mat, flagged_taxa)
    taxon_confidence <- taxon_confidence %>% filter(!taxon %in% flagged_taxa)
    cat("  --drop_low_confidence_taxa: REMOVED", paste(flagged_taxa, collapse=", "),
        "-- CLR re-centered on the remaining", ncol(clr_mat), "taxa:",
        paste(colnames(clr_mat), collapse=", "), "\n")
    write.csv(data.frame(taxon = flagged_taxa),
              file.path(opt$out_dir, paste0("dropped_low_confidence_taxa_", SYSTEM_TAG, "_", community, suffix, ".csv")),
              row.names = FALSE)
  } else if (opt$drop_low_confidence_taxa) {
    cat("  --drop_low_confidence_taxa: no flagged taxa in this basis -- nothing dropped\n")
  }

  samp <- sub_long %>%
    distinct(sample_uid, treatment, dpi, collection_well, growth_plate, seq_plate) %>%
    rename(timepoint = dpi)  # kept as `timepoint` internally to match 05b's variable names/plot code
  samp <- samp[match(rownames(clr_mat), samp$sample_uid), ]
  stopifnot(identical(rownames(clr_mat), samp$sample_uid))

  samp$treatment <- factor(recode(as.character(samp$treatment), !!!TREATMENT_RELABEL), levels=TREATMENT_LEVELS)
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
            file.path(opt$out_dir, paste0("permanova_model_planta_", community, suffix, ".rds")))

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

  # ── betadisper: location vs. spread.  --
  # Ordinary free permutation throughout (destructive design, every
  # row independent, no stratum restriction needed here unlike 07a).
  # NOTE: F/p from permutest() also collected into a betadisper_summary_*.csv
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
  loadings_df <- as.data.frame(pca$rotation) %>% rownames_to_column("taxon") %>%
    left_join(taxon_confidence, by="taxon")
  write.csv(loadings_df, file.path(opt$out_dir, paste0("pca_loadings_planta_", community, suffix, ".csv")), row.names=FALSE)

  p <- ggplot(scores, aes(x=PC1, y=PC2, colour=treatment, shape=timepoint)) +
    geom_point(size=3.2, alpha=0.9) +
    scale_colour_manual(values=TREATMENT_COLORS, drop=FALSE) +
    PLOT_THEME +
    labs(title=paste0("Aitchison PCA -- in planta | ", community, if (!is.na(opt$dpi_subset)) paste0(" | dpi ", opt$dpi_subset) else ""),
         subtitle=paste0("PC1 ", round(var_exp[1],1), "%, PC2 ", round(var_exp[2],1), "% of variance."),
         x=paste0("PC1 (", round(var_exp[1],1), "%)"),
         y=paste0("PC2 (", round(var_exp[2],1), "%)"))
  f_png <- file.path(plot_dir, paste0("pca_planta_", community, suffix, ".png"))
  ggsave(f_png, p, width=8, height=6)
  cat("  wrote", f_png, "\n")
  # Vector PDF companion -- same convention as 07a and 06_tailocin_effect_plot.py.
  f_pdf <- file.path(plot_dir, paste0("pca_planta_", community, suffix, ".pdf"))
  ggsave(f_pdf, p, width=8, height=6, device=cairo_pdf)
  cat("  wrote", f_pdf, "\n")

  scree_df <- tibble(PC = factor(paste0("PC", seq_along(var_exp)), levels=paste0("PC", seq_along(var_exp))),
                      pct_variance = var_exp, cum_pct_variance = cumsum(var_exp))
  write.csv(scree_df, file.path(opt$out_dir, paste0("pca_variance_explained_planta_", community, suffix, ".csv")), row.names=FALSE)

  p_scree <- ggplot(scree_df, aes(x=PC)) +
    geom_col(aes(y=pct_variance), fill=SCREE_BAR_COLOR) +
    geom_line(aes(y=cum_pct_variance, group=1), colour=SCREE_LINE_COLOR, linewidth=0.9) +
    geom_point(aes(y=cum_pct_variance), colour=SCREE_LINE_COLOR, size=2.2) +
    scale_y_continuous(limits=c(0, 100), expand=expansion(mult=c(0, 0.05))) +
    PLOT_THEME +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=FONT$scree_pc_labels)) +
    labs(title=paste0("Scree plot -- in planta | ", community, if (!is.na(opt$dpi_subset)) paste0(" | dpi ", opt$dpi_subset) else ""),
         x="Principal component", y="% variance explained")
  f_scree <- file.path(plot_dir, paste0("scree_planta_", community, suffix, ".png"))
  ggsave(f_scree, p_scree, width=8, height=5)
  cat("  wrote", f_scree, "\n")
  f_scree_pdf <- file.path(plot_dir, paste0("scree_planta_", community, suffix, ".pdf"))
  ggsave(f_scree_pdf, p_scree, width=8, height=5, device=cairo_pdf)
  cat("  wrote", f_scree_pdf, "\n\n")
}

cat("Done.\n")
