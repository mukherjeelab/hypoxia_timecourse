# Render 39_gc3_length_models.Rmd once per TE contrast (sictrl vs si3d).
# Usage: Rscript code/translational_control_overview/39_render_contrasts.R [hypoxia_1hr|normoxia_1and4hr ...]
# Needs RSTUDIO_PANDOC set (see project notes).
library(here)
contrasts <- list(
  hypoxia_1hr = list(de_file = "translation_categories_si3d_vs_sictrl_hypoxia_1hr.csv",
                     label = "MDA-MB-231 si3d hypoxia 1hr", tag = "hypoxia_1hr", te_word = "hypoxia 1hr"),
  normoxia_1and4hr = list(de_file = "translation_categories_si3d_vs_sictrl_normoxia_1and4hr.csv",
                          label = "MDA-MB-231 si3d normoxia 1+4hr (time-regressed)", tag = "normoxia_1and4hr",
                          te_word = "normoxia 1+4hr, time-regressed"))
args <- commandArgs(trailingOnly = TRUE)
todo <- if (length(args)) args else names(contrasts)
stopifnot("unknown contrast" = all(todo %in% names(contrasts)))
for (nm in todo) {
  rmarkdown::render(here("code", "translational_control_overview", "39_gc3_length_models.Rmd"),
                    params = contrasts[[nm]], output_file = paste0("39_gc3_length_models_", nm, ".html"), quiet = TRUE)
  cat("rendered", nm, "\n")
}
