# Renders 15_shap_heatmap.Rmd for the eIF3d cap-binding (subunit-seq) columns ONLY, kept apart
# from the combined TE heatmap on purpose. All columns are intrinsic-only. The MDA-MB-231 si3d
# hypoxia 1hr TE column is included as a reference.
#
# PROVISIONAL: the binding outcomes are computed from deposited RPKM, not DESeq2 on counts, and
# must be rerun from the raw reads. Raw enrichment is length-confounded (rho 0.58-0.78) and
# neither study has a background IP. Everything is in CAP_BINDING_DATA.md.
#   Rscript code/translational_control_overview/15c_render_subunit_seq.R
T <- "_lfc0.5_tcoreg_intrinsic_corrected_both_clipaggregate_yraw"
e <- c("MDA-MB-231 3d hypoxia 1hr (TE reference)" = "eif33d_promotes_hypoxia_1hr",
       "HEK293T cap binding, complete media"      = "eif33d_promotes_hek293t_capbind_cm_rpkm",
       "HEK293T cap binding, glucose deprivation" = "eif33d_promotes_hek293t_capbind_glu_rpkm",
       "HEK293T cap binding, thapsigargin"        = "eif33d_promotes_hek293t_capbind_tg_rpkm",
       "HEK293T cap binding, glucose / complete"  = "eif33d_promotes_hek293t_capbind_gluvscm_rpkm",
       "HEK293T cap binding, thapsigargin / DMSO" = "eif33d_promotes_hek293t_capbind_tgvsdmso_rpkm")
exps <- paste(paste0(names(e), "=", e, T), collapse = ";")
rmarkdown::render(here::here("code/translational_control_overview/15_shap_heatmap.Rmd"),
  params = list(experiments = exps, allow_mixed_genes = TRUE, reduced_suffix = "",
                clusters_suffix = "eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg_intrinsic_corrected_both_clipaggregate_yraw",
                out_tag = "subunit_seq_rpkm"),
  output_file = "15_shap_heatmap_subunit_seq_rpkm.html",
  intermediates_dir = tempfile(), quiet = TRUE)
