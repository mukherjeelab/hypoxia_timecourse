# Renders 15_shap_heatmap.Rmd with every MDA-MB-231 and MCF7-SIX1 regression run already on disk
# (8 MDA-MB-231 + 4 MCF7-SIX1 on its own isoforms + 2 HeLa from Herrmannova et al. 2024). Gene populations differ between cell lines,
# hence allow_mixed_genes = TRUE; the notebook reports the overlap.
#   Rscript code/translational_control_overview/15b_render_all_celllines.R full|intrinsic
T <- "_corrected_both_clipaggregate_yraw"
mk <- function(fs) {
  a <- c("3d hypoxia 1hr"="eif33d_promotes_hypoxia_1hr","3d hypoxia 4hr"="eif33d_promotes_hypoxia_4hr",
         "3d hypoxia 24hr"="eif33d_promotes_hypoxia_24hr","3d normoxia 1+4hr"="eif33d_promotes_normoxia_1and4hr",
         "3e hypoxia 1hr"="eif33e_promotes_hypoxia_1hr","3e hypoxia 4hr"="eif33e_promotes_hypoxia_4hr",
         "3e hypoxia 24hr"="eif33e_promotes_hypoxia_24hr","3e normoxia 1+4hr"="eif33e_promotes_normoxia_1and4hr")
  m <- c("3d hypoxia 1hr"="eif33d_promotes_mcf7six1_hypoxia_1hr","3d normoxia 1hr"="eif33d_promotes_mcf7six1_normoxia_1hr",
         "3e hypoxia 1hr"="eif33e_promotes_mcf7six1_hypoxia_1hr","3e normoxia 1hr"="eif33e_promotes_mcf7six1_normoxia_1hr")
  h <- c("3d no stress (Herrmannova)"="eif33d_promotes_hela_normoxia_steadystate",
         "3e no stress (Herrmannova)"="eif33e_promotes_hela_normoxia_steadystate")
  fx <- function(stem, lab) paste0(stem, "_lfc0.5_tcoreg", lab, T)
  e <- c(setNames(fx(a, fs), paste("MDA-MB-231", names(a))),
         # MCF7-SIX1 on its OWN dominant isoforms only (the _mcf7tx runs). The earlier fits on the
         # 231 isoform matrix differ from these by little and are not shown.
         setNames(fx(m, paste0("_mcf7tx", fs)), paste("MCF7-SIX1", names(m))),
         # HeLa (Herrmannova et al. 2024, eLife): gene-level TE, tested genes only, features on
         # HeLa's own dominant isoforms (the _helatx runs; Roiuk control RNA). Outcome files and
         # sign checks are in 10h_herrmannova_outcome_prep.Rmd.
         setNames(fx(h, paste0("_helatx", fs)), paste("HeLa", names(h))))
  have <- function(x) file.exists(here::here("output", "translational_control_overview",
                                             paste0("rfreg_importance_", x, ".csv")))
  # Intrinsic-only columns, each on its own dataset's transcripts; added once their run exists.
  if (fs == "_intrinsic") {
    extra <- c("HeLa 4E TE level, E.V. (Teleman)"     = fx("eif34e_promotes_hela_teleman_ev_level", "_helatx_intrinsic"),
               "HeLa 4E TE level, 4E-BP1 (Teleman)"   = fx("eif34e_promotes_hela_teleman_bp_level", "_helatx_intrinsic"),
               "HEK293T DHX29 IP enrichment (Hia)"    = fx("eif3dhx29_promotes_hek293t_dhx29ip_ribo", "_hektx_intrinsic"))
    e <- c(e, extra[vapply(extra, have, logical(1))])
  }
  # DAP5 knockout (Weber 2022, HEK293T; 10w_weber_dap5_outcome_prep.Rmd), the one column whose
  # features are on MANE Select rather than the dataset's own most abundant transcript - hence the
  # star and the footnote. Added only once its run exists, so the other columns still render.
  d <- fx("eif3dap5_promotes_hek293t_normoxia_steadystate", paste0("_manetx", fs))
  # Intrinsic heatmap only: the full DAP5 run excludes its three DAP5-derived features (circular
  # against a DAP5 outcome), so it has 66 features where every other full column has 69, and 15_
  # requires one shared feature set. Its full-model SHAP is in its own 11_ report.
  if (fs == "_intrinsic" && have(d)) e <- c(e, "HEK293T DAP5 no stress (Weber)*" = d)
  if (fs == "_intrinsic" && exists("with_capbind") && with_capbind) {
    cb <- c("HEK293T cap binding, complete media*" = "cm", "HEK293T cap binding, glucose deprivation*" = "glu",
            "HEK293T cap binding, thapsigargin*" = "tg", "HEK293T cap binding, glucose / complete*" = "gluvscm",
            "HEK293T cap binding, thapsigargin / DMSO*" = "tgvsdmso")
    cbs <- setNames(fx(paste0("eif33d_promotes_hek293t_capbind_", cb, "_rpkm"), "_manetx_intrinsic"), names(cb))
    stopifnot("a cap-binding run is missing" = all(vapply(cbs, have, logical(1))))
    e <- c(e, cbs)
  }
  paste(paste0(names(e), "=", e), collapse = ";")
}
FOOT <- paste("* features computed on the MANE Select transcript; all other columns use each",
              "dataset's own most abundant or representative transcript.")
args <- commandArgs(TRUE); fs <- if (args[1] %in% c("intrinsic", "intrinsic_capbind")) "_intrinsic" else ""
# "intrinsic_capbind" adds the five subunit-seq cap-binding columns (HEK293T, MANE Select, starred)
# to the intrinsic figure. They are kept out of the default figures: provisional, RPKM-based.
with_capbind <- args[1] == "intrinsic_capbind"
tag <- if (fs == "") "all_celllines" else if (with_capbind) "all_celllines_intrinsic_capbind" else "all_celllines_intrinsic"
rmarkdown::render(here::here("code/translational_control_overview/15_shap_heatmap.Rmd"),
  params = list(experiments = mk(fs), allow_mixed_genes = TRUE, reduced_suffix = "",
                # Correlated-feature clusters computed on THIS model's feature set: the intrinsic
                # heatmap uses 14_ run on the 57 intrinsic features, so no group counts an
                # external feature the model never saw (memberships verified identical otherwise).
                clusters_suffix = paste0("eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg", fs, T),
                out_tag = tag,
                footnote = if (grepl("*", mk(fs), fixed = TRUE)) FOOT else ""),
  output_file = paste0("15_shap_heatmap_", tag, ".html"),
  intermediates_dir = tempfile(), quiet = TRUE, envir = new.env())
# The model-free companions (Spearman rho per feature; top-100 mean values) use the same columns,
# intrinsic only by design. Each column's matrix comes from 16_/17_'s matrix_map by suffix tag.
if (fs == "_intrinsic") for (nb in c("16_feature_correlation_heatmap", "17_top100_feature_means_heatmap"))
  rmarkdown::render(here::here("code/translational_control_overview", paste0(nb, ".Rmd")),
    params = list(experiments = mk(fs), out_tag = tag,
                  footnote = if (grepl("*", mk(fs), fixed = TRUE)) FOOT else ""),
    output_file = paste0(nb, "_", tag, ".html"), intermediates_dir = tempfile(), quiet = TRUE, envir = new.env())
