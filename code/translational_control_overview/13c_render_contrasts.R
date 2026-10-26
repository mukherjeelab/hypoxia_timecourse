# 13c_render_contrasts.R - knit the GC + length linear-model notebook once per contrast.
#   Rscript code/translational_control_overview/13c_render_contrasts.R            # all contrasts
#   Rscript code/translational_control_overview/13c_render_contrasts.R normoxia_1and4hr
# Each contrast writes its own plots / CSVs (tagged by `tag`) and its own HTML, so runs never
# overwrite each other. The notebook itself checks that every param belongs to its tag.
suppressPackageStartupMessages(library(here))
Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")
rf <- "_lfc0.5_tcoreg_intrinsic_corrected_both_clipaggregate_yraw"
runs <- list(
  hypoxia_1hr = list(
    tag = "hypoxia_1hr", de_file = "translation_categories_si3d_vs_sictrl_hypoxia_1hr.csv",
    contrast_kind = "knockdown", label = "MDA-MB-231 si3d hypoxia 1hr",
    effect = "eIF3d-dependent TE change", outcome_label = "TE LFC (sictrl vs si3d)",
    positive_means = "TE is lower without eIF3d, i.e. eIF3d promotes translation",
    rf_suffix = paste0("eif33d_promotes_hypoxia_1hr", rf)),
  normoxia_1and4hr = list(
    tag = "normoxia_1and4hr", de_file = "translation_categories_si3d_vs_sictrl_normoxia_1and4hr.csv",
    contrast_kind = "knockdown", label = "MDA-MB-231 si3d normoxia 1+4hr (time-regressed)",
    effect = "eIF3d-dependent TE change", outcome_label = "TE LFC (sictrl vs si3d)",
    positive_means = "TE is lower without eIF3d, i.e. eIF3d promotes translation",
    rf_suffix = paste0("eif33d_promotes_normoxia_1and4hr", rf)),
  sictrl_hypvsnor_1hr = list(
    tag = "sictrl_hypvsnor_1hr", de_file = "sictrl_hypoxia_vs_normoxia_te_1hr.csv",
    contrast_kind = "condition", label = "MDA-MB-231 siCTRL hypoxia vs normoxia 1hr",
    effect = "Hypoxia-induced TE change", outcome_label = "TE LFC (hypoxia vs normoxia)",
    positive_means = "TE is higher in hypoxia than normoxia (siCTRL cells, no knockdown)",
    rf_suffix = paste0("eif3ctrl_promotes_hypvsnor_1hr", rf)))

which <- commandArgs(trailingOnly = TRUE)
if (length(which) == 0) which <- names(runs)
stopifnot("unknown contrast" = all(which %in% names(runs)))
nb <- here("code", "translational_control_overview", "13c_linear_model_gc_length.Rmd")
for (w in which) {
  cat("\n=== ", w, " ===\n", sep = "")
  rmarkdown::render(nb, params = c(runs[[w]], save_plots = TRUE),
                    output_file = paste0("13c_linear_model_gc_length_", w, ".html"),
                    envir = new.env(), quiet = TRUE)
}
