#!/usr/bin/env Rscript
# =============================================================================
# 09a_permanova_invitro.R
#
# PERMANOVA + betadisper + Aitchison PCA for the IN VITRO trial.
#
# Ported from 05a_permanova_invitro.R to consume the current pipeline's CLR
# output (07_clr_transform.py's long-format CSVs) instead of the old
# clr_wide_<method>_....rds files. All of the statistical design from 05a is
# UNCHANGED -- same split-plot restricted permutation, same plate-
# identifiability checks, same betadisper-per-stratum logic. Only the input
# format and a few metadata column names changed. See "WHAT CHANGED FROM 05a"
# below for the specifics.
#
# ---------------------------------------------------------------------------
# SPLIT-PLOT DESIGN (unchanged from 05a -- read this if you haven't seen it)
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
# dropped upstream (e.g. by 07's coverage/prevalence filtering), a sub-culture
# can end up with an incomplete timepoint set. This script detects that
# automatically and drops any incomplete sub-culture, reporting exactly which
# ones and why, rather than crashing -- see "enforce balance" below.
#
# ---------------------------------------------------------------------------
# SEQ_PLATE AS A COVARIATE -- read before trusting the plate p-value
# ---------------------------------------------------------------------------
# Re-verify this against your actual data with
# table(meta$sequencing_plate, meta$timepoint) if plate assignments change --
# it's checked at runtime below (check_plate_identifiability()), not
# hardcoded, but the STRUCTURE of the confound (which timepoints share a
# plate) is a property of your specific plate layout and can shift between
# reruns of 00b_build_master_metadata_invitro.py if the layout sheet changes.
#
# ---------------------------------------------------------------------------
# WHY t0 / CONTROL / MONOCULTURE WELLS AREN'T HERE
# ---------------------------------------------------------------------------
# 07_clr_transform.py's in-vitro CLR pass already restricts to
# sample_type == 'experimental' -- t0, monoculture, and blank wells never
# appear in the CLR long CSVs this script reads, so there's no equivalent of
# 05a's --include_t0 flag to carry over. If you want those wells included in
# a PERMANOVA for some reason, that filtering needs to change in
# 07_clr_transform.py, not here.
#
# ---------------------------------------------------------------------------
# WHAT CHANGED FROM 05a
# ---------------------------------------------------------------------------
# - Input: --clr_file points directly at ONE long-format CSV (sample_uid,
#   taxon, clr_value, community_background, treatment, timepoint, BIOREP_ID,
#   well, sequencing_plate) instead of --clr_dir + --method + a wide .rds
#   per community. Use clr_invitro_standard_long.csv for the full-basis (old
#   Fig 3) run, or clr_invitro_p25C2_excluded_long.csv for the indirect-
#   effect (old Fig 4, --exclude_suffix _excl-p25C2) run -- these are NOT
#   interchangeable with the dual-basis clr_invitro_long.csv used by the
#   effect-size plot script (08_tailocin_effect_plot.py); that file mixes two
#   different CLR bases in one table (p25C2 on one basis, everything else on
#   another) and is NOT valid input for a distance-based test where every
#   taxon needs to be on the same basis. See 07_clr_transform.py's docstring.
# - Communities are inferred from community_background values present in the
#   file (S1, S1+C2) rather than from filenames -- no --exclude_suffix
#   needed, that's now just "which CLR file did you point --clr_file at."
# - No --metadata / --exclude_dropped: dropped_flag filtering doesn't apply
#   to in vitro (no such concept in this design) and isn't needed here.
# - No --include_t0: see above, t0 was never in the input to begin with.
# - biorep_id -> BIOREP_ID, treatment -> treatment (unchanged), seq_plate ->
#   sequencing_plate (column rename only, same meaning).
# - stat_ellipse() REMOVED from the PCA plot per request -- points only, no
#   drawn confidence ellipse. betadisper()/permutest() are UNCHANGED and
#   still run and reported; removing the ellipse is a plotting choice, not a
#   statistical one, and doesn't touch the dispersion test at all.
# - NEW: --timepoint_subset lets you rerun the analysis restricted to one
#   timepoint (e.g. --timepoint_subset 24H) for a "which single timepoint
#   should the actual figure use" comparison. At a single timepoint each
#   sub-culture contributes exactly one well, so the split-plot restriction
#   is moot -- this mode uses ordinary free-permutation adonis2(d~treatment),
#   still with betadisper (no ellipse) and its own PCA. The default (no
#   --timepoint_subset) still produces the full split-plot analysis across
#   all four timepoints, matching 05a's original scope -- both are useful,
#   run both.
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
#          <out_dir>/plots/scree_invitro_<community>[_<timepoint>].png
#
# Usage:
#   Rscript 09a_permanova_invitro.R --clr_file clr_out/clr_invitro_standard_long.csv --out_dir permanova_invitro
#   Rscript 09a_permanova_invitro.R --clr_file clr_out/clr_invitro_standard_long.csv --out_dir permanova_invitro --timepoint_subset 24H
#   Rscript 09a_permanova_invitro.R --clr_file clr_out/clr_invitro_p25C2_excluded_long.csv --out_dir permanova_invitro_excl_p25C2
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
  make_option("--seed", type="integer", default=1)
)
opt <- parse_args(OptionParser(option_list=opt_list))
if (is.null(opt$clr_file)) stop("--clr_file is required (see header for which file to point it at).")
opt$out_dir <- sub("/+$", "", opt$out_dir)
plot_dir <- file.path(opt$out_dir, "plots")
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)
set.seed(opt$seed)

# ── plate identifiability check (unchanged from 05a) ────────────────────────
check_plate_identifiability <- function(samp, plate_col, group_col) {
  tab <- table(samp[[plate_col]], samp[[group_col]])
  informative <- rowSums(tab > 0) > 1
  list(tab = tab, fully_aliased = !any(informative),
       n_informative_plates = sum(informative), n_plates = nrow(tab))
}

# ══════════════════════════════════════════════════════════════════════════
# 1. Load the long CLR table and pivot to one wide matrix per community
# ══════════════════════════════════════════════════════════════════════════
clr_long <- read.csv(opt$clr_file, stringsAsFactors=FALSE)
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
#    (S1 and S1+C2 are on different compositional bases whenever --clr_file
#    is the p25C2-excluded table for S1+C2 wells -- see 07_clr_transform.py's
#    dual-CLR handling. Always analyzed separately, same as 05a.)
# ══════════════════════════════════════════════════════════════════════════
for (community in communities) {
  suffix <- if (!is.na(opt$timepoint_subset)) paste0("_", opt$timepoint_subset) else ""
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

  samp <- sub_long %>%
    distinct(sample_uid, treatment, timepoint, BIOREP_ID, well, sequencing_plate) %>%
    rename(seq_plate = sequencing_plate, biorep_id = BIOREP_ID)
  samp <- samp[match(rownames(clr_mat), samp$sample_uid), ]
  stopifnot(identical(rownames(clr_mat), samp$sample_uid))

  samp$treatment  <- factor(samp$treatment)
  samp$timepoint  <- factor(samp$timepoint)
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
    # FULL MODE: split-plot restricted permutation across all timepoints,
    # unchanged from 05a.
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
                  paste0("permanova_model_invitro_seqplate_sensitivity_", tp_sub, "_", community, ".rds")))
        }
      }
    }

    saveRDS(list(wholeplot=mod_wholeplot, timepoint=mod_timepoint, subplot=mod_subplot, seqplate=mod_seqplate),
            file.path(opt$out_dir, paste0("permanova_model_invitro_", community, ".rds")))

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
    write.csv(combined, file.path(opt$out_dir, paste0("permanova_invitro_", community, ".csv")), row.names=FALSE)

    disp_facs <- c("treatment", "timepoint")
    if (!is.null(mod_seqplate)) disp_facs <- c(disp_facs, "seq_plate")
    disp_h <- list(treatment=h_wholeplot, timepoint=h_subplot, seq_plate=h_subplot)
  }

  # ── betadisper: is a significant PERMANOVA a location or spread effect? ──
  # Kept exactly as in 05a -- this is a statistical test, not a plot element,
  # and runs regardless of the ellipse being removed from the PCA plot below.
  # NOTE: previously only per-sample distances were written to file and the
  # actual F/p from permutest() only went to the console (cat()) -- that
  # meant the significance test itself wasn't reproducible from the output
  # files alone. Now also collected into a betadisper_summary_*.csv per
  # community/suffix, one row per factor.
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
  write.csv(as.data.frame(pca$rotation), file.path(opt$out_dir, paste0("pca_loadings_invitro_", community, suffix, ".csv")))

  # NOTE: stat_ellipse() layer REMOVED per request -- points only. betadisper
  # above is untouched; this is a plotting-only change.
  p <- ggplot(scores, aes(x=PC1, y=PC2, colour=treatment, shape=timepoint)) +
    geom_point(size=2.2, alpha=0.85) +
    theme_minimal(base_size=10) +
    labs(title=paste0("Aitchison PCA -- in vitro | ", community, if (!is.na(opt$timepoint_subset)) paste0(" | ", opt$timepoint_subset) else ""),
         subtitle=paste0("PC1 ", round(var_exp[1],1), "%, PC2 ", round(var_exp[2],1), "% of variance."),
         x=paste0("PC1 (", round(var_exp[1],1), "%)"),
         y=paste0("PC2 (", round(var_exp[2],1), "%)"))
  f_png <- file.path(plot_dir, paste0("pca_invitro_", community, suffix, ".png"))
  ggsave(f_png, p, width=8, height=6)
  cat("  wrote", f_png, "\n")

  # ── scree plot ────────────────────────────────────────────────────────
  scree_df <- tibble(PC = factor(paste0("PC", seq_along(var_exp)), levels=paste0("PC", seq_along(var_exp))),
                      pct_variance = var_exp, cum_pct_variance = cumsum(var_exp))
  write.csv(scree_df, file.path(opt$out_dir, paste0("pca_variance_explained_invitro_", community, suffix, ".csv")), row.names=FALSE)

  p_scree <- ggplot(scree_df, aes(x=PC)) +
    geom_col(aes(y=pct_variance), fill="grey70") +
    geom_line(aes(y=cum_pct_variance, group=1), colour="firebrick", linewidth=0.7) +
    geom_point(aes(y=cum_pct_variance), colour="firebrick", size=1.6) +
    scale_y_continuous(limits=c(0, 100), expand=expansion(mult=c(0, 0.05))) +
    theme_minimal(base_size=10) +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7)) +
    labs(title=paste0("Scree plot -- in vitro | ", community, if (!is.na(opt$timepoint_subset)) paste0(" | ", opt$timepoint_subset) else ""),
         x="Principal component", y="% variance explained")
  f_scree <- file.path(plot_dir, paste0("scree_invitro_", community, suffix, ".png"))
  ggsave(f_scree, p_scree, width=8, height=5)
  cat("  wrote", f_scree, "\n\n")
}

cat("Done.\n")
