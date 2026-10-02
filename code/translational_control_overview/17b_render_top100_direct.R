# 17b_render_top100_direct.R
# Renders the top-100 mean-feature-value heatmap (17_) straight from the outcome files, without
# waiting for any regression run: each column's genes and LFC come from its outcome table joined
# to its own feature matrix (17_'s `outcomes` param). Same columns and labels as the combined
# intrinsic SHAP heatmap (15b_), so the two figures can be read side by side.
#   Rscript code/translational_control_overview/17b_render_top100_direct.R [matrix tags]
# matrix tags limits the figure to matrices that are ready, e.g. "231" or "231,mcf7tx,helatx".
# Default: every column whose matrix and outcome file exist.
suppressMessages({ library(here) })
O <- here("output", "translational_control_overview")
a <- commandArgs(TRUE)
tags  <- if (length(a) >= 1 && a[1] != "all") strsplit(a[1], ",")[[1]] else NULL
top_n <- if (length(a) >= 2) as.integer(a[2]) else 100L   # 100 / 200 / 500 (Kate, 2026-10-01)
# A third argument "rho" renders the Spearman-rho heatmap (16_) for the same columns instead.
which_nb <- if (length(a) >= 3 && a[3] == "rho") "16" else "17"
tc <- function(g, c_, t) paste0("translation_categories_si", g, "_vs_sictrl_", c_, "_", t, ".csv")
cols <- list(
  # label, matrix tag (as in 15b's suffixes), outcome file, LFC column, padj column
  list("MDA-MB-231 3d hypoxia 1hr",    "231",     tc("3d", "hypoxia", "1hr"),  "te_lfc", "te_padj"),
  list("MDA-MB-231 3d hypoxia 4hr",    "231",     tc("3d", "hypoxia", "4hr"),  "te_lfc", "te_padj"),
  list("MDA-MB-231 3d hypoxia 24hr",   "231",     tc("3d", "hypoxia", "24hr"), "te_lfc", "te_padj"),
  list("MDA-MB-231 3d normoxia 1+4hr", "231",     tc("3d", "normoxia", "1and4hr"), "te_lfc", "te_padj"),
  list("MDA-MB-231 3e hypoxia 1hr",    "231",     tc("3e", "hypoxia", "1hr"),  "te_lfc", "te_padj"),
  list("MDA-MB-231 3e hypoxia 4hr",    "231",     tc("3e", "hypoxia", "4hr"),  "te_lfc", "te_padj"),
  list("MDA-MB-231 3e hypoxia 24hr",   "231",     tc("3e", "hypoxia", "24hr"), "te_lfc", "te_padj"),
  list("MDA-MB-231 3e normoxia 1+4hr", "231",     "te_si3e_vs_sictrl_normoxia_1and4hr.csv", "te_lfc", "te_padj"),
  list("MCF7-SIX1 3d hypoxia 1hr",     "mcf7tx",  tc("3d", "mcf7six1_hypoxia", "1hr"),  "te_lfc", "te_padj"),
  list("MCF7-SIX1 3d normoxia 1hr",    "mcf7tx",  tc("3d", "mcf7six1_normoxia", "1hr"), "te_lfc", "te_padj"),
  list("MCF7-SIX1 3e hypoxia 1hr",     "mcf7tx",  tc("3e", "mcf7six1_hypoxia", "1hr"),  "te_lfc", "te_padj"),
  list("MCF7-SIX1 3e normoxia 1hr",    "mcf7tx",  tc("3e", "mcf7six1_normoxia", "1hr"), "te_lfc", "te_padj"),
  list("HeLa 3d no stress (Herrmannova)", "helatx", tc("3d", "hela_normoxia", "steadystate"), "te_lfc", "te_padj"),
  list("HeLa 3e no stress (Herrmannova)", "helatx", tc("3e", "hela_normoxia", "steadystate"), "te_lfc", "te_padj"),
  list("HeLa 4E TE level, E.V. (Teleman)",   "helatx", "teleman_absolute_te_ev.csv", "te_level", "te_level_padj"),
  list("HeLa 4E TE level, 4E-BP1 (Teleman)", "helatx", "teleman_absolute_te_bp.csv", "te_level", "te_level_padj"),
  list("HEK293T DHX29 IP enrichment (Hia)",  "hektx",  "dhx29_ip_hek293t_enrichment.csv", "enrich_lfc", "enrich_padj"),
  list("HEK293T DAP5 no stress (Weber)*",    "manetx", "translation_categories_sidap5_vs_sictrl_hek293t_normoxia_steadystate.csv", "te_lfc", "te_padj"))
mat_of <- c(`231` = "feature_matrix_tco.rds", mcf7tx = "feature_matrix_tco_mcf7six1.rds",
            helatx = "feature_matrix_tco_hela.rds", hektx = "feature_matrix_tco_dhx29.rds",
            manetx = "feature_matrix_tco_mane.rds")
keep <- Filter(function(z) (is.null(tags) || z[[2]] %in% tags) &&
                 file.exists(here("output", z[[3]])) && file.exists(file.path(O, mat_of[[z[[2]]]])), cols)
skipped <- setdiff(vapply(cols, `[[`, "", 1), vapply(keep, `[[`, "", 1))
if (length(skipped)) cat("not included:", paste(skipped, collapse = "; "), "\n")
# 17_ reads BOTH the matrix (matrix_map) and the annotation bars - cell line, initiation factor,
# condition - from the suffix, so it must be a real run name in 15b_'s form, not a placeholder.
# (A placeholder suffix rendered every column as MDA-MB-231 / eIF3d / "other".)
stem <- c("MDA-MB-231 3d hypoxia 1hr" = "eif33d_promotes_hypoxia_1hr", "MDA-MB-231 3d hypoxia 4hr" = "eif33d_promotes_hypoxia_4hr",
          "MDA-MB-231 3d hypoxia 24hr" = "eif33d_promotes_hypoxia_24hr", "MDA-MB-231 3d normoxia 1+4hr" = "eif33d_promotes_normoxia_1and4hr",
          "MDA-MB-231 3e hypoxia 1hr" = "eif33e_promotes_hypoxia_1hr", "MDA-MB-231 3e hypoxia 4hr" = "eif33e_promotes_hypoxia_4hr",
          "MDA-MB-231 3e hypoxia 24hr" = "eif33e_promotes_hypoxia_24hr", "MDA-MB-231 3e normoxia 1+4hr" = "eif33e_promotes_normoxia_1and4hr",
          "MCF7-SIX1 3d hypoxia 1hr" = "eif33d_promotes_mcf7six1_hypoxia_1hr", "MCF7-SIX1 3d normoxia 1hr" = "eif33d_promotes_mcf7six1_normoxia_1hr",
          "MCF7-SIX1 3e hypoxia 1hr" = "eif33e_promotes_mcf7six1_hypoxia_1hr", "MCF7-SIX1 3e normoxia 1hr" = "eif33e_promotes_mcf7six1_normoxia_1hr",
          "HeLa 3d no stress (Herrmannova)" = "eif33d_promotes_hela_normoxia_steadystate",
          "HeLa 3e no stress (Herrmannova)" = "eif33e_promotes_hela_normoxia_steadystate",
          "HeLa 4E TE level, E.V. (Teleman)" = "eif34e_promotes_hela_teleman_ev_level",
          "HeLa 4E TE level, 4E-BP1 (Teleman)" = "eif34e_promotes_hela_teleman_bp_level",
          "HEK293T DHX29 IP enrichment (Hia)" = "eif3dhx29_promotes_hek293t_dhx29ip_ribo",
          "HEK293T DAP5 no stress (Weber)*" = "eif3dap5_promotes_hek293t_normoxia_steadystate")
stopifnot("a column has no run-name stem" = all(vapply(keep, `[[`, "", 1) %in% names(stem)))
sfx <- vapply(keep, function(z) paste0(stem[[z[[1]]]], "_lfc0.5_tcoreg",
                                        if (z[[2]] == "231") "" else paste0("_", z[[2]]), "_intrinsic_direct"), "")
exps <- paste(paste0(vapply(keep, `[[`, "", 1), "=", sfx), collapse = ";")
outs <- paste(vapply(keep, function(z) paste0(z[[1]], "=", z[[3]], "|", z[[4]], "|", z[[5]]), ""), collapse = ";")
star <- any(grepl("*", vapply(keep, `[[`, "", 1), fixed = TRUE))
FOOT <- paste("* features computed on the MANE Select transcript; all other columns use each",
              "dataset's own most abundant or representative transcript.")
tag  <- paste0("direct_", if (is.null(tags)) "all" else paste(tags, collapse = "_"))
if (which_nb == "16") {
  rmarkdown::render(here("code/translational_control_overview/16_feature_correlation_heatmap.Rmd"),
    params = list(experiments = exps, outcomes = outs, out_tag = tag, footnote = if (star) FOOT else ""),
    output_file = paste0("16_feature_correlation_heatmap_", tag, ".html"), intermediates_dir = tempfile(), quiet = TRUE, envir = new.env())
  cat("rendered", length(keep), "columns -> plots/16_feature_lfc_correlation_heatmap_", tag, ".pdf\n", sep = "")
  quit(save = "no")
}
rmarkdown::render(here("code/translational_control_overview/17_top100_feature_means_heatmap.Rmd"),
  params = list(experiments = exps, outcomes = outs, out_tag = tag, top_n = top_n,
                footnote = if (star) FOOT else ""),
  output_file = paste0("17_top", top_n, "_feature_means_heatmap_", tag, ".html"),
  intermediates_dir = tempfile(), quiet = TRUE, envir = new.env())
cat("rendered", length(keep), "columns ->", paste0("plots/17_top", top_n, "_feature_means_heatmap_", tag, ".pdf"), "\n")
