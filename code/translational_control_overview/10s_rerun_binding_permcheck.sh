#!/usr/bin/env bash
# Re-run 10_ for the five subunit-seq columns with the 20-shuffle permutation check (perm_n = 20),
# then 11_ only for glucose deprivation (the one column that failed the old single-shuffle bound
# and so has no SHAP yet), then the separate subunit-seq heatmap. The real forests are unchanged
# (same seeds), so existing SHAP outputs for the other four stay valid.
set -uo pipefail
cd "$(cd "$(dirname "$0")/../.." && pwd)"
export RSTUDIO_PANDOC="${RSTUDIO_PANDOC:-/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools}"
J=code/translational_control_overview/jobs_subunit_seq_rpkm.txt
render() {  # $1 notebook, $2 job line
  local spec="${2%%|*}" extra="${2#*|}"; set -- "$1" $spec
  local W; W="$(mktemp -d)"
  Rscript --vanilla -e "
    nb <- 'code/translational_control_overview/$1.Rmd'
    p <- list(gene='$2', condition='$3', timepoint='$4', run_cv=FALSE, save_plots=FALSE)
    $extra
    p <- p[intersect(names(p), names(rmarkdown::yaml_front_matter(nb)\$params))]
    ok <- tryCatch({rmarkdown::render(nb, params=p, output_file=file.path('$W','o.html'), intermediates_dir='$W', quiet=TRUE); 'OK'},
                   error=function(e) paste('FAILED:', conditionMessage(e)))
    cat(sprintf('%-22s %-28s %s\n', '$1', '$3', ok))" 2>&1 | grep -vE "Warning|built under"
  rm -rf "$W"
}
while IFS= read -r line; do [ -n "$line" ] && render 10_rf_regression "$line"; done < "$J"
render 11_shap_regression "$(grep 'capbind_glu ' "$J")"
Rscript code/translational_control_overview/15c_render_subunit_seq.R 2>&1 | grep -vE "Warning|built under" | tail -3
echo "heatmap: $(ls -l code/translational_control_overview/15_shap_heatmap_subunit_seq_rpkm.html 2>&1)"
