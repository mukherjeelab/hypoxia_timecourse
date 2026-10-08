# 15f_render_231_mcf7six1.R
# SHAP heatmaps for our own data only: 8 MDA-MB-231 + 4 MCF7-SIX1 columns (MCF7-SIX1 on its own
# dominant isoforms, the _mcf7tx runs), full or intrinsic. Same notebook and settings as 15b_.
#   Rscript code/translational_control_overview/15f_render_231_mcf7six1.R full|intrinsic
T <- "_corrected_both_clipaggregate_yraw"
a <- c("3d hypoxia 1hr"="eif33d_promotes_hypoxia_1hr","3d hypoxia 4hr"="eif33d_promotes_hypoxia_4hr",
       "3d hypoxia 24hr"="eif33d_promotes_hypoxia_24hr","3d normoxia 1+4hr"="eif33d_promotes_normoxia_1and4hr",
       "3e hypoxia 1hr"="eif33e_promotes_hypoxia_1hr","3e hypoxia 4hr"="eif33e_promotes_hypoxia_4hr",
       "3e hypoxia 24hr"="eif33e_promotes_hypoxia_24hr","3e normoxia 1+4hr"="eif33e_promotes_normoxia_1and4hr")
m <- c("3d hypoxia 1hr"="eif33d_promotes_mcf7six1_hypoxia_1hr","3d normoxia 1hr"="eif33d_promotes_mcf7six1_normoxia_1hr",
       "3e hypoxia 1hr"="eif33e_promotes_mcf7six1_hypoxia_1hr","3e normoxia 1hr"="eif33e_promotes_mcf7six1_normoxia_1hr")
fs <- if (commandArgs(TRUE)[1] == "intrinsic") "_intrinsic" else ""
fx <- function(stem, lab) paste0(stem, "_lfc0.5_tcoreg", lab, T)
e <- c(setNames(fx(a, fs), paste("MDA-MB-231", names(a))),
       setNames(fx(m, paste0("_mcf7tx", fs)), paste("MCF7-SIX1", names(m))))
tag <- paste0("mdamb231_mcf7six1", fs)
rmarkdown::render(here::here("code/translational_control_overview/15_shap_heatmap.Rmd"),
  params = list(experiments = paste(paste0(names(e), "=", e), collapse = ";"), allow_mixed_genes = TRUE,
                reduced_suffix = "", clusters_suffix = paste0("eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg", fs, T),
                out_tag = tag),
  output_file = paste0("15_shap_heatmap_", tag, ".html"), intermediates_dir = tempfile(), quiet = TRUE, envir = new.env())
