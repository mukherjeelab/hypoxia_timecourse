# 15d_render_231_only.R
# MDA-MB-231-only SHAP heatmaps (8 columns, full or intrinsic), for reading the MDA-MB-231 results
# while the other cell lines' reruns are still going. Same notebook and settings as 15b_.
#   Rscript code/translational_control_overview/15d_render_231_only.R full|intrinsic
T <- "_corrected_both_clipaggregate_yraw"
a <- c("3d hypoxia 1hr"="eif33d_promotes_hypoxia_1hr","3d hypoxia 4hr"="eif33d_promotes_hypoxia_4hr",
       "3d hypoxia 24hr"="eif33d_promotes_hypoxia_24hr","3d normoxia 1+4hr"="eif33d_promotes_normoxia_1and4hr",
       "3e hypoxia 1hr"="eif33e_promotes_hypoxia_1hr","3e hypoxia 4hr"="eif33e_promotes_hypoxia_4hr",
       "3e hypoxia 24hr"="eif33e_promotes_hypoxia_24hr","3e normoxia 1+4hr"="eif33e_promotes_normoxia_1and4hr")
fs <- if (commandArgs(TRUE)[1] == "intrinsic") "_intrinsic" else ""
e <- setNames(paste0(a, "_lfc0.5_tcoreg", fs, T), paste("MDA-MB-231", names(a)))
tag <- paste0("mdamb231_only", fs)
rmarkdown::render(here::here("code/translational_control_overview/15_shap_heatmap.Rmd"),
  params = list(experiments = paste(paste0(names(e), "=", e), collapse = ";"), allow_mixed_genes = TRUE,
                reduced_suffix = "", clusters_suffix = paste0("eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg", fs, T),
                out_tag = tag),
  output_file = paste0("15_shap_heatmap_", tag, ".html"), intermediates_dir = tempfile(), quiet = TRUE, envir = new.env())
