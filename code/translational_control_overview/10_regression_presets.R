# Shared definitions for the REGRESSION fork (10_ / 10b_ / 10c_ / 11_ / 11d_ / 11f_ / 12_).
#
# Mirrors the role 02_sweep_presets.R plays for the classification stack, and exists for the
# same reason: a sweep is only interpretable if the arms that produced it are still on disk.
# 02_sweep_presets.R is SOURCED by the regression stack, never edited - 02b_/02c_/03_/03f_ all
# depend on it, and sweep_suffix() is how 12_ resolves the classifier run it compares against.
#
# The replication axis differs from the classification stack and that is the whole point.
# Notebook 02_ keeps every positive and draws a size-matched negative class, so neg_seed drives
# the noise (test AUC sd 0.0165). A regression on all 10,172 genes has NO draw, so neg_seed has
# nothing to seed. split_seed replaces it: it seeds ONLY createDataPartition, while the forest
# and the CV stay pinned at forest_seed, so within a seed every arm sees identical rows and
# paired differences isolate the arm alone.
#
# CAVEAT that must travel with every regression floor: two 80/20 splits share ~81% of their
# training rows, against the negative draw's ~52% gene overlap. So the regression's floor is
# tighter than the classifier's, and shap_noise_floor.csv is NOT a valid yardstick for
# regression shares. 11f_ measures shapreg_noise_floor.csv for that purpose.

REG_SPLIT_SEEDS <- 1:20

# --- output namespace -------------------------------------------------------------------
# The classification stack owns these prefixes. A regression run that wrote one of them would
# silently destroy a frozen baseline whose two gates are green, so this is an ASSERTED
# invariant rather than a naming convention. Every write path in 10_ and 11_ goes through it,
# and the output_prefix_collide mutation in 10z_ proves it fires.
CLF_PREFIXES <- c("rf_classifier_", "rf_importance_", "rf_performance_",
                  "rf_model_gene_sets_", "rf_model_data_", "shap_")

# The guard rests entirely on "shapreg_..." not matching the "shap_" prefix (the underscore
# sits in a different position). Assert that here rather than trusting it.
stopifnot(
  "the namespace guard is broken: shapreg_ matches the shap_ prefix" =
    !startsWith("shapreg_feature_x.csv", "shap_"),
  "the namespace guard is broken: rfreg_ matches a classification prefix" =
    !any(vapply(CLF_PREFIXES, function(p) startsWith("rfreg_importance_x.csv", p), logical(1)))
)

assert_no_collision <- function(paths) {
  bn  <- basename(paths)
  bad <- bn[vapply(bn, function(b) any(vapply(CLF_PREFIXES,
                    function(p) startsWith(b, p), logical(1))), logical(1))]
  stopifnot(
    "a regression output would land in the classification namespace" = length(bad) == 0,
    "regression must never write to output/predictive_modeling" =
      !any(grepl("output/predictive_modeling", paths, fixed = TRUE))
  )
  invisible(paths)
}

# Write helper that cannot bypass the guard.
reg_write_csv <- function(x, path) { assert_no_collision(path); readr::write_csv(x, path); path }
reg_save_rds  <- function(x, path) { assert_no_collision(path); saveRDS(x, path);          path }

# --- suffix -----------------------------------------------------------------------------
# Deliberately string-alignable with sweep_suffix() so 12_ can pair a regression run with the
# classifier run it is being compared against. Everything up to _clip{mode} is identical; the
# regression-only tokens are appended after it, and each appears ONLY when it is not the
# default, so adding a param never renames an existing run's outputs.
#
# The "eif33d" double-3 is inherited from 02_rf_model.Rmd's paste0("eif3", gene, ...) with
# gene = "3d". It is reproduced deliberately, not a typo: changing it here would break the
# string alignment with every classification output on disk.
REG_SUFFIX_DEFAULTS <- list(
  gene = "3d", direction = "promotes", readout = "TE", condition = "hypoxia",
  timepoint = "1hr", lfc = "0.5", feature_set_label = "tcoreg",
  variant_set = "corrected", initiation_score = "both", clip_mode = "aggregate",
  target_transform = "raw", winsorize_y = 0, min_node_size = 5, mtry = 0,
  num_trees = 1000, splitrule = "variance", sample_fraction = 0, weight_mode = "none",
  restrict_to_classifier_genes = FALSE, forest_seed = 9
)

reg_suffix <- function(arm = list(), split_seed = 9L) {
  p <- utils::modifyList(REG_SUFFIX_DEFAULTS, arm)
  paste0(
    "eif3", p$gene, "_", p$direction, "_",
    if (p$readout != "TE") paste0(tolower(p$readout), "_") else "",
    p$condition, "_", p$timepoint, "_lfc", p$lfc,
    "_", p$feature_set_label,
    if (p$variant_set != "current") paste0("_", p$variant_set) else "",
    if (p$initiation_score != "pwm") paste0("_", p$initiation_score) else "",
    "_clip", p$clip_mode,
    "_y", gsub("_", "", p$target_transform),
    if (p$winsorize_y > 0) paste0("_w", p$winsorize_y) else "",
    if (p$min_node_size != 5) paste0("_mns", p$min_node_size) else "",
    if (p$mtry != 0) paste0("_mtry", p$mtry) else "",
    if (p$num_trees != 1000) paste0("_t", p$num_trees) else "",
    if (p$splitrule != "variance") paste0("_sr", p$splitrule) else "",
    # sample_fraction and forest_seed are live params of 10_rf_regression.Rmd and both change the
    # forest, so a run at a non-default value must not land on the default run's filenames.
    # Before they were tokenised, sample_fraction = 0.632 silently overwrote rfreg_model_,
    # rfreg_importance_, rfreg_performance_ and every downstream shapreg_* of the default run.
    # 0 means "let ranger choose" (1 with replacement for regression), which is the default.
    if (p$sample_fraction != 0) paste0("_sf", p$sample_fraction) else "",
    if (p$weight_mode != "none") paste0("_wt", p$weight_mode) else "",
    if (isTRUE(p$restrict_to_classifier_genes)) "_clfgenes" else "",
    if (as.integer(p$forest_seed) != 9L) paste0("_fseed", as.integer(p$forest_seed)) else "",
    if (as.integer(split_seed) != 9L) paste0("_split", as.integer(split_seed)) else ""
  )
}

# A family sweep is always the same three arms, and the shuffled arm is the point: impurity
# importance is non-negative, so an uninformative family never scores ZERO. Whatever the
# permuted copy still absorbs is the floor the real family has to clear. Permuting preserves
# every marginal, the NA pattern and the within-family correlation structure.
#
# Not reused from 02_sweep_presets.R's family_arms() because the regression arms need their own
# stem convention and that file is sourced read-only by four other notebooks.
family_arms_reg <- function(families, stem) {
  list(
    none = list(feature_set_label = paste0(stem, "none"), exclude_families = families),
    real = list(feature_set_label = paste0(stem, "real")),
    shuf = list(feature_set_label = paste0(stem, "shuf"), shuffle_families = families)
  )
}

# --- presets ----------------------------------------------------------------------------
# Each arm carries its OWN feature_set_label. Two arms sharing one overwrite each other's
# outputs, which is the bug 02_sweep_presets.R's header records.
REG_PRESETS <- list(
  # The question that motivated the fork. raw / winsor3 / rank test the SAME hypothesis - tail
  # robustness - at three strengths, which is what makes the result decisive rather than
  # merely a preference: if winsor3 recovers most of rank's Spearman gain, the gain was the 8
  # extreme genes and the cheaper interpretable arm wins; if it recovers little, rank's gain is
  # genuine non-linearity across the whole distribution.
  transform = list(
    raw    = list(feature_set_label = "tcoreg_tf", target_transform = "raw"),
    winsor = list(feature_set_label = "tcoreg_tf", target_transform = "raw", winsorize_y = 3),
    rank   = list(feature_set_label = "tcoreg_tf", target_transform = "rank_ecdf")
  ),
  # min.node.size is the ONLY ranger default that differs between probability and regression
  # forests on this matrix (5 vs 10); mtry is 8 = floor(sqrt(69)) for BOTH. So this preset is
  # small on purpose. mtry 23 = floor(p/3) is randomForest's regression default, included to
  # show it is not better rather than because ranger would ever pick it.
  hyperparams = list(
    mns5   = list(feature_set_label = "tcoreg_hp", min_node_size = 5),
    mns10  = list(feature_set_label = "tcoreg_hp", min_node_size = 10),
    mtryp3 = list(feature_set_label = "tcoreg_hp", mtry = 23),
    t500   = list(feature_set_label = "tcoreg_hp", num_trees = 500)
  ),
  # Separates "different model type" from "different gene population" in the SHAP comparison.
  # Without this arm, agreement between the classifier and the regression cannot be told apart
  # from two effects cancelling.
  population = list(
    all     = list(feature_set_label = "tcoreg_pop"),
    clfgene = list(feature_set_label = "tcoreg_pop", restrict_to_classifier_genes = TRUE)
  ),
  # Family nulls transfer unchanged: shuffle_families still works on a regression, and
  # impurity importance is still non-negative, so the shuffled arm is still the floor.
  g4      = family_arms_reg("g4", "tcoreg_g4"),
  peptide = family_arms_reg("nascent_peptide,polya_track", "tcoreg_pep")
)

PRESET_FAMILIES_REG <- list(
  transform   = NA_character_,
  hyperparams = NA_character_,
  population  = NA_character_,
  g4          = "g4",
  peptide     = c("nascent_peptide", "polya_track")
)

# --- one representative per correlated cluster -------------------------------------------
# Builds the exclude_features list from 14_feature_correlation.Rmd's cluster table instead of
# having it typed by hand, and guarantees the rule that matters for a cross-experiment heatmap:
# EVERY cluster keeps exactly one survivor, never zero.
#
# Why never zero: a feature that looks like noise in one experiment may carry signal in another,
# and if it is absent here its row is missing from this column while present in that one, so the
# rows stop lining up. That is the whole thing the heatmap depends on. Dropping a whole cluster
# also discards the concept, not just the redundancy.
#
# Representative choice, in order: largest SHAP share, then most central within the cluster
# (highest mean |rho| to the other members), then fewest missing values. SHAP leads because it is
# what the model actually used; centrality breaks ties, which for a 2-member cluster is always a
# tie since both members share the single correlation.
reduce_to_representatives <- function(clusters_csv, shap_csv, threshold = 0.8,
                                      centrality = NULL) {
  stopifnot("cluster table not found - run 14_feature_correlation.Rmd first" =
              file.exists(clusters_csv),
            "SHAP table not found - run 11_shap_regression.Rmd first" = file.exists(shap_csv))
  cl <- readr::read_csv(clusters_csv, show_col_types = FALSE)
  cl <- cl[abs(cl$threshold - threshold) < 1e-9, , drop = FALSE]
  stopifnot("no clusters at that threshold" = nrow(cl) > 0)
  sh <- readr::read_csv(shap_csv, show_col_types = FALSE)
  cl$share <- sh$share_pct[match(cl$feature, sh$feature)]
  if (!is.null(centrality)) cl$centrality <- centrality[cl$feature] else cl$centrality <- 0

  keep <- vapply(split(cl, cl$cluster), function(d) {
    o <- order(-d$share, -d$centrality, d$feature)
    d$feature[o[1]]
  }, character(1))
  drop <- setdiff(cl$feature, keep)

  # The invariant, asserted rather than trusted.
  per_cluster_kept <- vapply(split(cl, cl$cluster),
                             function(d) sum(!d$feature %in% drop), integer(1))
  stopifnot(
    "a cluster would lose every member - never drop a whole cluster" = all(per_cluster_kept == 1L),
    "a kept feature is also in the drop list" = !any(keep %in% drop)
  )
  list(keep = unname(keep), drop = drop,
       n_clusters = length(keep), n_dropped = length(drop))
}
