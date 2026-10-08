#!/bin/bash
# 32_run_si3e_clustering.sh - out-of-fold SHAP + per-timepoint clustering (23_) for eIF3e in
# MDA-MB-231: hypoxia 1, 4, 24hr and normoxia 1+4hr, intrinsic model - the si3e counterpart of
# 25b_run_timecourse_shap.sh, so the two knockdowns' per-gene SHAP can be compared directly.
# The 10_ regression models already exist; this only runs 23_. Negative controls are the 02b si3e
# set (same rule as the si3d set 23_ uses by default). 23_ caches folds, so a rerun resumes.
# ~4-6 h per timepoint.
#   bash code/translational_control_overview/32_run_si3e_clustering.sh
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
export RSTUDIO_PANDOC="${RSTUDIO_PANDOC:-/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools}"
TCO=code/translational_control_overview
STEM="eif33e_promotes_%s_lfc0.5_tcoreg_intrinsic_corrected_both_clipaggregate_yraw"
for t in hypoxia_1hr hypoxia_4hr hypoxia_24hr normoxia_1and4hr; do
  [ -f "output/translational_control_overview/rfreg_model_data_$(printf "$STEM" "$t").rds" ] || { echo "FAILED: no 10_ model data for $t"; exit 1; }
done
run23() {  # timepoint  label  positive-geneset ("" = derive from the outcome)
  local sfx; sfx=$(printf "$STEM" "$1")
  echo "=== 23_ si3e $1  $(date +%H:%M)"
  Rscript --vanilla -e "rmarkdown::render('$TCO/23_shap_gene_clustering.Rmd',
    params = list(suffix = '$sfx', label = '$2', pos_geneset = '$3', pos_label = 'eIF3e promoted',
                  neg_geneset = 'predictive_modeling/negative_control_genes_si3e.csv', lfc_label = 'LFC sictrl/si3e'),
    output_file = '23_shap_gene_clustering_si3e_$1.html', intermediates_dir = tempfile(), quiet = TRUE)" \
    || { echo "FAILED 23_ si3e $1"; exit 1; }
}
run23 hypoxia_1hr      "MDA-MB-231 si3e hypoxia 1hr"    hypoxia_3e_promotes_TE_1hr_lfc0.5.csv
run23 hypoxia_4hr      "MDA-MB-231 si3e hypoxia 4hr"    hypoxia_3e_promotes_TE_4hr_lfc0.5.csv
run23 hypoxia_24hr     "MDA-MB-231 si3e hypoxia 24hr"   hypoxia_3e_promotes_TE_24hr_lfc0.5.csv
run23 normoxia_1and4hr "MDA-MB-231 si3e normoxia 1+4hr" ""
echo "=== ALL DONE  $(date +%H:%M)"
