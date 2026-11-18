#!/usr/bin/env Rscript
# =============================================================================
# 07a_permanova_invitro.R
#
# PERMANOVA + betadisper + Aitchison PCA for the IN VITRO trial.
#
# ---------------------------------------------------------------------------
# SPLIT-PLOT DESIGN
# ---------------------------------------------------------------------------
# The repeatedly-sampled unit is NOT biorep_id alone -- each biorep_id is
# split into three independently treated sub-cultures (Buffer/delTFA/
# Tailocin), and it's each of those sub-cultures (biorep_id x treatment) that
# gets non-destructively resampled across the 4 timepoints. subculture_id =
# interaction(biorep_id, treatment) is the correct plot unit. Two separately-
# scoped permutation schemes are required, not one (confirmed empirically in
# 05a: a single combined free/free scheme silently dropped the
# treatment/timepoint rows and only returned the interaction):
#   * whole-plot stratum (tests treatment): swap ENTIRE sub-cultures against
#     each other; never permute timepoint order within one.
#   * sub-plot stratum (tests timepoint AND treatment:timepoint): hold each
#     sub-culture's treatment fixed; freely permute timepoint order within it.
# This is the standard Anderson & ter Braak restricted-permutation framework
# for split-plot designs.
#
# Restricted permutation via Plots(strata=...) requires every stratum to have
# the SAME sub-structure (same set/count of timepoints). If a well was
# dropped upstream (e.g. by 05's coverage/prevalence filtering), a sub-culture
# can end up with an incomplete timepoint set. This script detects that
# automatically and drops any incomplete sub-culture, reporting exactly which
# ones and why, rather than crashing -- see "enforce balance" below.
#
# Inputs : --clr_file <path to clr_invitro_standard_long.csv or
#           clr_invitro_p25C2_excluded_long.csv, from 07_clr_transform.py>
# Outputs: <out_dir>/permanova_invitro_<community>[_<timepoint>].csv
#          <out_dir>/betadisper_invitro_<community>[_<timepoint>]_by_<factor>.csv
#          <out_dir>/betadisper_summary_invitro_<community>[_<timepoint>].csv  (F/p per factor, NEW)
#          <out_dir>/pca_scores_invitro_<community>[_<timepoint>].csv
#          <out_dir>/pca_loadings_invitro_<community>[_<timepoint>].csv
#          <out_dir>/pca_variance_explained_invitro_<community>[_<timepoint>].csv
#          <out_dir>/plots/pca_invitro_<community>[_<timepoint>].png
#          <out_dir>/plots/pca_invitro_<community>[_<timepoint>].pdf   (NEW, vector)
#          <out_dir>/plots/scree_invitro_<community>[_<timepoint>].png
#          <out_dir>/plots/scree_invitro_<community>[_<timepoint>].pdf (NEW, vector)
#
# Usage:
#   Rscript 07a_permanova_invitro.R --clr_file clr_out/clr_invitro_standard_long.csv --out_dir permanova_invitro
#   Rscript 07a_permanova_invitro.R --clr_file clr_out/clr_invitro_standard_long.csv --out_dir permanova_invitro --timepoint_subset 24H
#   Rscript 07a_permanova_invitro.R --clr_file clr_out/clr_invitro_p25C2_excluded_long.csv --out_dir permanova_invitro_excl_p25C2
# =============================================================================

suppressPackageStartupMessages({
  library(vegan); library(permute); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(optparse)
})

opt_list <- list(
  make_option("--clr_file",  type="character",
              help="Path to clr_invitro_standard_long.csv or clr_invitro_p25C2_excluded_long.csv (from 07_clr_transform.py)"),
  make_option("--out_dir",   type="character", default="permanova_invitro"),
  make_option("--n_perm",    type="integer", default=999),
  make_option("--timepoint_subset", type="character", default=NA_character_,
              help="If set (e.g. '24H'), restrict to that single timepoint and use ordinary free permutation instead of the split-plot scheme. Default: NA (full split-plot analysis across all timepoints)."),
  make_option("--communities", type="character", default=NA_character_,
              help="Comma-separated community_background values to analyze, e.g. 'S1,S1+C2'. Default: all values present in --clr_file."),
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

# ── plate identifiability check ──────────────────────────────────────
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
SYSTEM_TAG <- "invitro"

# ══════════════════════════════════════════════════════════════════════════
# 1. Load the long CLR table and pivot to one wide matrix per community
# ══════════════════════════════════════════════════════════════════════════
clr_long <- read.csv(opt$clr_file, stringsAsFactors=FALSE)

as_flag <- function(x) toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1")
if ("low_confidence_flag" %in% names(clr_long)) {
  clr_long$low_confidence_flag <- as_flag(clr_long$low_confidence_flag)
} else {
  if (opt$drop_low_confidence_taxa) stop("--drop_low_confidence_taxa set, but --clr_file has no low_confidence_flag column")
  clr_long$low_confidence_flag <- FALSE
}
required_cols <- c("sample_uid","taxon","clr_value","community_background","treatment","timepoint")
missing_cols <- setdiff(required_cols, names(clr_long))
if (length(missing_cols)) stop("--clr_file is missing required column(s): ", paste(missing_cols, collapse=", "))

communities <- if (!is.na(opt$communities)) strsplit(opt$communities, ",")[[1]] else unique(clr_long$community_background)
cat("Communities to analyze:", paste(communities, collapse=", "), "\n\n")

if (!is.na(opt$timepoint_subset))
  cat("NOTE: --timepoint_subset =", opt$timepoint_subset,
      "-- running SINGLE-TIMEPOINT mode (ordinary free permutation, no split-plot restriction needed).\n\n")

# ══════════════════════════════════════════════════════════════════════════
# 2. One PERMANOVA + betadisper + PCA per community
# ══════════════════════════════════════════════════════════════════════════
for (community in communities) {
  suffix <- if (!is.na(opt$timepoint_subset)) paste0("_", opt$timepoint_subset) else ""
  if (opt$drop_low_confidence_taxa) suffix <- paste0(suffix, "_dropLC")
  cat("=== in vitro |", community, if (!is.na(opt$timepoint_subset)) paste0("| ", opt$timepoint_subset) else "", "===\n")

  sub_long <- clr_long %>% filter(community_background == community)
  if (!is.na(opt$timepoint_subset)) sub_long <- sub_long %>% filter(timepoint == opt$timepoint_subset)

  if (nrow(sub_long) == 0) { cat("  SKIP -- no rows for this community/timepoint\n\n"); next }

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

  # Report which taxa in this basis are flagged low_confidence by
  # 05_clr_transform.py (blank-contamination-derived QC flag -- see that
  # script's docstring). By default this is informational only and every
  # taxon stays in the analysis. With --drop_low_confidence_taxa the flagged
  # taxa are removed and the CLR is re-centered on the remaining taxa (a
  # valid subcomposition -- see drop_low_conf_taxa() above), as a
  # sensitivity analysis.
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
    distinct(sample_uid, treatment, timepoint, BIOREP_ID, well, sequencing_plate) %>%
    rename(seq_plate = sequencing_plate, biorep_id = BIOREP_ID)
  samp <- samp[match(rownames(clr_mat), samp$sample_uid), ]
  stopifnot(identical(rownames(clr_mat), samp$sample_uid))

  samp$treatment  <- factor(samp$treatment, levels=TREATMENT_LEVELS)
  # numeric order (6H, 12H, 18H, 24H), not alphabetical -- affects legend order only, not any test
  tp_lv <- unique(as.character(samp$timepoint))
  samp$timepoint  <- factor(samp$timepoint, levels=tp_lv[order(as.numeric(sub("[^0-9.].*$", "", tp_lv)))])
  samp$biorep_id  <- factor(samp$biorep_id)
  samp$seq_plate  <- factor(samp$seq_plate)
  samp$subculture_id <- interaction(samp$biorep_id, samp$treatment, drop=TRUE)

  if (nrow(clr_mat) < 4 || nlevels(samp$treatment) < 2) {
    cat("  SKIP -- insufficient wells/factor levels after filtering\n\n"); next
  }
  d <- vegdist(clr_mat, method="euclidean")  # Aitchison distance = Euclidean on CLR

  # ────────────────────────────────────────────────────────────────────────
  # SINGLE-TIMEPOINT MODE: ordinary free permutation, treatment only.
  # ────────────────────────────────────────────────────────────────────────
  if (!is.na(opt$timepoint_subset)) {
    mod_treatment <- adonis2(d ~ treatment, data=samp, permutations=opt$n_perm, by="margin")
    cat("  [treatment, free permutation -- single timepoint, sub-cultures independent here]\n")
    print(mod_treatment)

    combined <- bind_rows(
      as.data.frame(mod_treatment) %>% rownames_to_column("term") %>% filter(term == "treatment") %>% mutate(stratum="treatment (free permutation, single timepoint)"),
      as.data.frame(mod_treatment) %>% rownames_to_column("term") %>% filter(term == "Residual") %>% mutate(stratum="treatment (free permutation, single timepoint)", term="Residual")
    )
    write.csv(combined, file.path(opt$out_dir, paste0("permanova_invitro_", community, suffix, ".csv")), row.names=FALSE)
    saveRDS(mod_treatment, file.path(opt$out_dir, paste0("permanova_model_invitro_", community, suffix, ".rds")))

    disp_facs <- c("treatment")
    disp_h <- opt$n_perm
  } else {
    # ────────────────────────────────────────────────────────────────────
    # FULL MODE: split-plot restricted permutation across all timepoints
    # ────────────────────────────────────────────────────────────────────
    full_tp_set <- levels(samp$timepoint)
    tp_by_sub <- split(as.character(samp$timepoint), samp$subculture_id)
    complete <- vapply(tp_by_sub, function(x) setequal(x, full_tp_set) && length(x) == length(full_tp_set), logical(1))
    incomplete_subs <- names(complete)[!complete]
    if (length(incomplete_subs)) {
      cat("  Dropping", length(incomplete_subs), "of", length(complete),
          "sub-cultures with incomplete timepoint coverage:", paste(incomplete_subs, collapse=", "), "\n")
      for (s in incomplete_subs)
        cat("    ", s, "has timepoints:", paste(sort(tp_by_sub[[s]]), collapse=","),
            "-- expected:", paste(full_tp_set, collapse=","), "\n")
      samp <- samp[samp$subculture_id %in% names(complete)[complete], ]
      samp$subculture_id <- droplevels(samp$subculture_id)
      samp$treatment <- droplevels(samp$treatment)
      samp$timepoint <- droplevels(samp$timepoint)
      clr_mat <- clr_mat[rownames(clr_mat) %in% samp$sample_uid, , drop=FALSE]
      samp <- samp[match(rownames(clr_mat), samp$sample_uid), ]
      d <- vegdist(clr_mat, method="euclidean")
    }

    if (nrow(clr_mat) < 4 || nlevels(samp$treatment) < 2 || nlevels(samp$timepoint) < 2) {
      cat("  SKIP -- insufficient wells/factor levels after balance enforcement\n\n"); next
    }

    h_wholeplot <- how(within = Within(type="none"),
                        plots  = Plots(strata=samp$subculture_id, type="free"),
                        nperm  = opt$n_perm)
    h_subplot <- how(within = Within(type="free"),
                      plots  = Plots(strata=samp$subculture_id, type="none"),
                      nperm  = opt$n_perm)

    mod_wholeplot <- adonis2(d ~ treatment, data=samp, permutations=h_wholeplot, by="margin")
    mod_timepoint <- adonis2(d ~ timepoint, data=samp, permutations=h_subplot, by="margin")
    mod_subplot   <- adonis2(d ~ treatment * timepoint, data=samp, permutations=h_subplot, by="margin")

    cat("  [whole-plot stratum -- treatment, tested by swapping whole sub-cultures]\n"); print(mod_wholeplot)
    cat("  [sub-plot stratum -- timepoint, tested within sub-cultures]\n"); print(mod_timepoint)
    cat("  [sub-plot stratum -- interaction, tested within sub-cultures]\n"); print(mod_subplot)

    plate_check <- check_plate_identifiability(samp, "seq_plate", "timepoint")
    cat("  [seq_plate x timepoint crosstab]\n"); print(plate_check$tab)
    mod_seqplate <- NULL
    if (nlevels(samp$seq_plate) < 2) {
      cat("  seq_plate: SKIPPED -- only one plate present after filtering\n")
    } else if (plate_check$fully_aliased) {
      cat("  seq_plate: SKIPPED -- fully aliased with timepoint; not estimable.\n")
    } else {
      if (plate_check$n_informative_plates == 1)
        cat("  seq_plate: fitting, but only", plate_check$n_informative_plates, "of",
            plate_check$n_plates, "plate(s) has samples at >1 timepoint -- result reflects\n",
            "  mostly that plate's internal split, not a broad cross-plate comparison.\n")
      mod_seqplate <- adonis2(d ~ seq_plate, data=samp, permutations=h_subplot, by="margin")
      cat("  [sub-plot stratum -- seq_plate, tested within sub-cultures]\n"); print(mod_seqplate)

      tp_counts <- colSums(plate_check$tab > 0)
      informative_tp <- names(tp_counts)[tp_counts > 1]
      if (length(informative_tp) >= 1) {
        tp_sub <- informative_tp[1]
        samp_sub <- samp[samp$timepoint == tp_sub, ]
        samp_sub$seq_plate <- droplevels(samp_sub$seq_plate)
        if (nlevels(samp_sub$seq_plate) >= 2 && nrow(samp_sub) >= 4) {
          d_sub <- vegdist(clr_mat[rownames(clr_mat) %in% samp_sub$sample_uid, , drop=FALSE], method="euclidean")
          mod_seqplate_sub <- adonis2(d_sub ~ seq_plate, data=samp_sub, permutations=opt$n_perm, by="margin")
          cat("  [SENSITIVITY: seq_plate within timepoint==", tp_sub, " only, ordinary free permutation]\n", sep="")
          print(mod_seqplate_sub)
          saveRDS(mod_seqplate_sub, file.path(opt$out_dir,
                  paste0("permanova_model_invitro_seqplate_sensitivity_", tp_sub, "_", community, suffix, ".rds")))
        }
      }
    }

    saveRDS(list(wholeplot=mod_wholeplot, timepoint=mod_timepoint, subplot=mod_subplot, seqplate=mod_seqplate),
            file.path(opt$out_dir, paste0("permanova_model_invitro_", community, suffix, ".rds")))

    combined <- bind_rows(
      as.data.frame(mod_wholeplot) %>% rownames_to_column("term") %>% filter(term == "treatment") %>% mutate(stratum="whole-plot (treatment)"),
      as.data.frame(mod_timepoint) %>% rownames_to_column("term") %>% filter(term == "timepoint") %>% mutate(stratum="sub-plot (timepoint)"),
      as.data.frame(mod_subplot) %>% rownames_to_column("term") %>% filter(term == "treatment:timepoint") %>% mutate(stratum="sub-plot (interaction)"),
      as.data.frame(mod_wholeplot) %>% rownames_to_column("term") %>% filter(term == "Residual") %>% mutate(stratum="whole-plot (treatment)", term="Residual (whole-plot)"),
      as.data.frame(mod_subplot) %>% rownames_to_column("term") %>% filter(term == "Residual") %>% mutate(stratum="sub-plot (timepoint, interaction)", term="Residual (sub-plot)")
    )
    if (!is.null(mod_seqplate)) {
      combined <- bind_rows(combined,
        as.data.frame(mod_seqplate) %>% rownames_to_column("term") %>% filter(term == "seq_plate") %>% mutate(stratum="sub-plot (seq_plate)"),
        as.data.frame(mod_seqplate) %>% rownames_to_column("term") %>% filter(term == "Residual") %>% mutate(stratum="sub-plot (seq_plate)", term="Residual (seq_plate model)"))
    }
    write.csv(combined, file.path(opt$out_dir, paste0("permanova_invitro_", community, suffix, ".csv")), row.names=FALSE)

    disp_facs <- c("treatment", "timepoint")
    if (!is.null(mod_seqplate)) disp_facs <- c(disp_facs, "seq_plate")
    disp_h <- list(treatment=h_wholeplot, timepoint=h_subplot, seq_plate=h_subplot)
  }

  # ── betadisper: is a significant PERMANOVA a location or spread effect? ──
  betadisper_summary_rows <- list()
  for (fac in disp_facs) {
    bd <- betadisper(d, samp[[fac]], type="centroid")
    h_this <- if (is.list(disp_h)) disp_h[[fac]] else disp_h
    pt <- permutest(bd, permutations=h_this)
    bd_out <- tibble(sample_uid=samp$sample_uid, group=samp[[fac]], dist_to_centroid=bd$distances)
    write.csv(bd_out, file.path(opt$out_dir, paste0("betadisper_invitro_", community, suffix, "_by_", fac, ".csv")), row.names=FALSE)
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
  write.csv(betadisper_summary, file.path(opt$out_dir, paste0("betadisper_summary_invitro_", community, suffix, ".csv")), row.names=FALSE)
  cat("  wrote", file.path(opt$out_dir, paste0("betadisper_summary_invitro_", community, suffix, ".csv")), "\n")

  # ── Aitchison PCA (ordinary PCA on CLR coordinates -- no extra scaling) ──
  pca <- prcomp(clr_mat, center=TRUE, scale.=FALSE)
  var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  scores <- as.data.frame(pca$x) %>% rownames_to_column("sample_uid") %>% left_join(samp, by="sample_uid")
  write.csv(scores, file.path(opt$out_dir, paste0("pca_scores_invitro_", community, suffix, ".csv")), row.names=FALSE)
  loadings_df <- as.data.frame(pca$rotation) %>% rownames_to_column("taxon") %>%
    left_join(taxon_confidence, by="taxon")
  write.csv(loadings_df, file.path(opt$out_dir, paste0("pca_loadings_invitro_", community, suffix, ".csv")), row.names=FALSE)

  p <- ggplot(scores, aes(x=PC1, y=PC2, colour=treatment, shape=timepoint)) +
    geom_point(size=3.2, alpha=0.9) +
    scale_colour_manual(values=TREATMENT_COLORS, drop=FALSE) +
    PLOT_THEME +
    labs(title=paste0("Aitchison PCA -- in vitro | ", community, if (!is.na(opt$timepoint_subset)) paste0(" | ", opt$timepoint_subset) else ""),
         subtitle=paste0("PC1 ", round(var_exp[1],1), "%, PC2 ", round(var_exp[2],1), "% of variance."),
         x=paste0("PC1 (", round(var_exp[1],1), "%)"),
         y=paste0("PC2 (", round(var_exp[2],1), "%)"))
  f_png <- file.path(plot_dir, paste0("pca_invitro_", community, suffix, ".png"))
  ggsave(f_png, p, width=8, height=6)
  cat("  wrote", f_png, "\n")

  f_pdf <- file.path(plot_dir, paste0("pca_invitro_", community, suffix, ".pdf"))
  ggsave(f_pdf, p, width=8, height=6, device=cairo_pdf)
  cat("  wrote", f_pdf, "\n")

  # ── scree plot ────────────────────────────────────────────────────────
  scree_df <- tibble(PC = factor(paste0("PC", seq_along(var_exp)), levels=paste0("PC", seq_along(var_exp))),
                      pct_variance = var_exp, cum_pct_variance = cumsum(var_exp))
  write.csv(scree_df, file.path(opt$out_dir, paste0("pca_variance_explained_invitro_", community, suffix, ".csv")), row.names=FALSE)

  p_scree <- ggplot(scree_df, aes(x=PC)) +
    geom_col(aes(y=pct_variance), fill=SCREE_BAR_COLOR) +
    geom_line(aes(y=cum_pct_variance, group=1), colour=SCREE_LINE_COLOR, linewidth=0.9) +
    geom_point(aes(y=cum_pct_variance), colour=SCREE_LINE_COLOR, size=2.2) +
    scale_y_continuous(limits=c(0, 100), expand=expansion(mult=c(0, 0.05))) +
    PLOT_THEME +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=FONT$scree_pc_labels)) +
    labs(title=paste0("Scree plot -- in vitro | ", community, if (!is.na(opt$timepoint_subset)) paste0(" | ", opt$timepoint_subset) else ""),
         x="Principal component", y="% variance explained")
  f_scree <- file.path(plot_dir, paste0("scree_invitro_", community, suffix, ".png"))
  ggsave(f_scree, p_scree, width=8, height=5)
  cat("  wrote", f_scree, "\n")
  f_scree_pdf <- file.path(plot_dir, paste0("scree_invitro_", community, suffix, ".pdf"))
  ggsave(f_scree_pdf, p_scree, width=8, height=5, device=cairo_pdf)
  cat("  wrote", f_scree_pdf, "\n\n")
}

cat("Done.\n")
