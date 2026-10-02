#!/bin/bash
# 25b_run_timecourse_shap.sh - out-of-fold SHAP + per-timepoint clustering (23_) for eIF3d hypoxia
# 4hr, 24hr and normoxia 1+4hr in MDA-MB-231 (1hr is already done), then the joint time-course
# clustering (25_). Each 23_ run caches its folds, so a rerun resumes. ~3-6 h per timepoint.
#   bash code/translational_control_overview/25b_run_timecourse_shap.sh
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
export RSTUDIO_PANDOC="${RSTUDIO_PANDOC:-/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools}"
TCO=code/translational_control_overview
STEM="eif33d_promotes_%s_lfc0.5_tcoreg_intrinsic_corrected_both_clipaggregate_yraw"
run23() {  # timepoint  label  positive-geneset ("" = derive from the outcome)
  local sfx; sfx=$(printf "$STEM" "$1")
  echo "=== 23_ $1  $(date +%H:%M)"
  Rscript --vanilla -e "rmarkdown::render('$TCO/23_shap_gene_clustering.Rmd',
    params = list(suffix = '$sfx', label = '$2', pos_geneset = '$3'),
    output_file = '23_shap_gene_clustering_$1.html', intermediates_dir = tempfile(), quiet = TRUE)" || { echo "FAILED 23_ $1"; exit 1; }
}
run23 hypoxia_4hr      "MDA-MB-231 si3d hypoxia 4hr"      hypoxia_3d_promotes_TE_4hr_lfc0.5.csv
run23 hypoxia_24hr     "MDA-MB-231 si3d hypoxia 24hr"     hypoxia_3d_promotes_TE_24hr_lfc0.5.csv
run23 normoxia_1and4hr "MDA-MB-231 si3d normoxia 1+4hr"   ""
echo "=== 25_ joint  $(date +%H:%M)"
Rscript --vanilla -e "rmarkdown::render('$TCO/25_shap_timecourse_joint_clustering.Rmd', intermediates_dir = tempfile(), quiet = TRUE)" || { echo "FAILED 25_"; exit 1; }
echo "=== ALL DONE  $(date +%H:%M)"
