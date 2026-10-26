# 26b_render_sets.R - knit the target-set SHAP clustering notebook once per target set.
#   Rscript code/translational_control_overview/26b_render_sets.R                  # all sets
#   Rscript code/translational_control_overview/26b_render_sets.R sictrl_hypvsnor_1hr
# Plots, CSVs and the null cache are keyed by `suffix`, and each set gets its own HTML, so runs
# never overwrite each other. The notebook checks each set against its rule and its saved CSV.
suppressPackageStartupMessages(library(here))
Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")
rf <- "_lfc0.5_tcoreg_intrinsic_corrected_both_clipaggregate_yraw"
runs <- list(
  hypoxia_1hr = list(
    suffix = paste0("eif33d_promotes_hypoxia_1hr", rf), label = "MDA-MB-231 si3d hypoxia 1hr",
    later_sets = "4hr=hypoxia_3d_promotes_TE_4hr_lfc0.5.csv,24hr=hypoxia_3d_promotes_TE_24hr_lfc0.5.csv",
    set_label = "eIF3d promoted", set_rule = "knockdown", set_word = "eIF3d targets",
    geneset_csv = "hypoxia_3d_promotes_TE_1hr_lfc0.5.csv", lfc_label = "LFC sictrl/si3d (cap 1.5)", k_report = 3,
    report_profile = "features, normalized"),
  sictrl_hypvsnor_1hr = list(
    suffix = paste0("eif3ctrl_promotes_hypvsnor_1hr", rf), label = "MDA-MB-231 siCTRL hypoxia vs normoxia 1hr",
    later_sets = "4hr=hypvsnor_ctrl_promotes_TE_4hr_lfc0.5.csv,24hr=hypvsnor_ctrl_promotes_TE_24hr_lfc0.5.csv",
    set_label = "resists hypoxic TE drop", set_rule = "relative", set_word = "hypoxia-promoted mRNAs",
    geneset_csv = "hypvsnor_ctrl_promotes_TE_1hr_lfc0.5.csv", lfc_label = "LFC hypoxia/normoxia (cap 1.5)", k_report = 3,
    report_profile = "features, raw"))

which <- commandArgs(trailingOnly = TRUE)
if (length(which) == 0) which <- names(runs)
stopifnot("unknown set" = all(which %in% names(runs)))
nb <- here("code", "translational_control_overview", "26b_target_subtypes_features.Rmd")
for (w in which) {
  cat("\n=== ", w, " ===\n", sep = "")
  rmarkdown::render(nb, params = c(runs[[w]], save_plots = TRUE),
                    output_file = paste0("26b_target_subtypes_features_", w, ".html"),
                    envir = new.env(), quiet = TRUE)
}
