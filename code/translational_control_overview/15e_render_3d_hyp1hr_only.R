# 15e_render_3d_hyp1hr_only.R
# Single-column SHAP heatmaps for MDA-MB-231 si3d hypoxia 1hr alone (full or intrinsic model). Same
# notebook and settings as 15d_; with one column the cells carry their numbers.
#   Rscript code/translational_control_overview/15e_render_3d_hyp1hr_only.R full|intrinsic
fs <- if (commandArgs(TRUE)[1] == "intrinsic") "_intrinsic" else ""
sfx <- paste0("eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg", fs, "_corrected_both_clipaggregate_yraw")
tag <- paste0("mdamb231_3d_hyp1hr", fs)
rmarkdown::render(here::here("code/translational_control_overview/15_shap_heatmap.Rmd"),
  params = list(experiments = paste0("MDA-MB-231 3d hypoxia 1hr", if (nzchar(fs)) ", intrinsic model" else ", full model", "=", sfx), allow_mixed_genes = TRUE,
                reduced_suffix = "", clusters_suffix = sfx, out_tag = tag),
  output_file = paste0("15_shap_heatmap_", tag, ".html"), intermediates_dir = tempfile(), quiet = TRUE, envir = new.env())
