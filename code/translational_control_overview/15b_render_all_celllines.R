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
         # HeLa (Herrmannova et al. 2024, eLife): gene-level TE on MDA-MB-231's isoforms, tested
         # genes only. Outcome files and sign checks are in 10h_herrmannova_outcome_prep.Rmd.
         setNames(fx(h, fs), paste("HeLa", names(h))))
  # DHX29 selective ribosome IP (Hia et al., HEK293T; 10x_dhx29_ip_outcome_prep.R) joins the
  # intrinsic heatmap once it has been run on its own HEK293T-isoform matrix. Its earlier run on
  # MDA-MB-231 isoforms was stopped deliberately, so there is no SHAP output to show yet:
  #   if (fs == "_intrinsic") e <- c(e, "HEK293T DHX29 IP enrichment (Hia)" =
  #     fx("eifdhx29_promotes_hek293t_dhx29ip_ribo", fs))
  paste(paste0(names(e), "=", e), collapse = ";")
}
args <- commandArgs(TRUE); fs <- if (args[1] == "intrinsic") "_intrinsic" else ""
tag <- if (fs == "") "all_celllines" else "all_celllines_intrinsic"
rmarkdown::render(here::here("code/translational_control_overview/15_shap_heatmap.Rmd"),
  params = list(experiments = mk(fs), allow_mixed_genes = TRUE, reduced_suffix = "",
                # Correlated-feature clusters computed on THIS model's feature set: the intrinsic
                # heatmap uses 14_ run on the 57 intrinsic features, so no group counts an
                # external feature the model never saw (memberships verified identical otherwise).
                clusters_suffix = paste0("eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg", fs, T),
                out_tag = tag),
  output_file = paste0("15_shap_heatmap_", tag, ".html"),
  intermediates_dir = tempfile(), quiet = TRUE)
