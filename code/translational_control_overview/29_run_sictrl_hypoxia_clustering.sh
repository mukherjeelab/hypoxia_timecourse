#!/bin/bash
# 29_run_sictrl_hypoxia_clustering.sh - SHAP gene clustering of the cells' OWN hypoxic TE response:
# MDA-MB-231 siCTRL hypoxia vs normoxia at 1, 4 and 24hr, intrinsic model, run exactly as the si3d
# time course (10_ + 11_ via the parallel runner, then out-of-fold SHAP + clustering in 23_).
# Outcomes come from 10y_sictrl_hypoxia_outcome_prep.R. ~3-6 h per 23_ run; 23_ caches folds, so a
# rerun resumes.
#   bash code/translational_control_overview/29_run_sictrl_hypoxia_clustering.sh
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
export RSTUDIO_PANDOC="${RSTUDIO_PANDOC:-/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools}"
TCO=code/translational_control_overview
STEM="eif3ctrl_promotes_hypvsnor_%s_lfc0.5_tcoreg_intrinsic_corrected_both_clipaggregate_yraw"
echo "=== outcome prep  $(date +%H:%M)"
Rscript --vanilla $TCO/10y_sictrl_hypoxia_outcome_prep.R || { echo "FAILED 10y"; exit 1; }
echo "=== 10_ + 11_  $(date +%H:%M)"
# Only the timepoints whose 10_ model data is not already on disk, so a relaunch does not refit.
TODO=$(mktemp)
while IFS= read -r line; do
  t=$(echo "$line" | awk '{print $3}')
  [ -f "output/translational_control_overview/rfreg_model_data_$(printf "$STEM" "$t").rds" ] || echo "$line" >> "$TODO"
done < $TCO/jobs_sictrl_hypvsnor.txt
[ -s "$TODO" ] && bash $TCO/10p_parallel_runner.sh "$TODO" 3
rm -f "$TODO"
for t in 1hr 4hr 24hr; do
  [ -f "output/translational_control_overview/rfreg_model_data_$(printf "$STEM" "$t").rds" ] || { echo "FAILED: no 10_ model data for $t"; exit 1; }
done
for t in 1hr 4hr 24hr; do
  sfx=$(printf "$STEM" "$t")
  echo "=== 23_ $t  $(date +%H:%M)"
  Rscript --vanilla -e "rmarkdown::render('$TCO/23_shap_gene_clustering.Rmd',
    params = list(suffix = '$sfx', label = 'MDA-MB-231 siCTRL hypoxia vs normoxia $t',
                  pos_geneset = 'hypvsnor_ctrl_promotes_TE_${t}_lfc0.5.csv', pos_label = 'resists hypoxic TE drop',
                  neg_geneset = '', lfc_label = 'TE log2(hyp/nor)'),
    output_file = '23_shap_gene_clustering_sictrl_hypvsnor_$t.html', intermediates_dir = tempfile(), quiet = TRUE)" \
    || { echo "FAILED 23_ $t"; exit 1; }
done
echo "=== ALL DONE  $(date +%H:%M)"
