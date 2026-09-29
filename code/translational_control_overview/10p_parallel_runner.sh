#!/bin/bash
# Render 10_rf_regression.Rmd and 11_shap_regression.Rmd for many contrasts concurrently.
#
#   ./10p_parallel_runner.sh <jobfile> [max_concurrent]
#
# Each line of <jobfile> is:   <gene> <condition> <timepoint> | <extra R assignments>
# e.g.
#   3d hypoxia 1hr |
#   3e hypoxia 1hr | p$auc_bridge <- FALSE
#   3e normoxia 1and4hr | p$auc_bridge <- FALSE; p$outcome_csv <- 'te_si3e_vs_sictrl_normoxia_1and4hr.csv'
#
# WHY THIS EXISTS. knitr writes <notebook>.knit.md next to the .Rmd, so two concurrent renders of
# the same notebook delete each other's intermediate and one dies in pandoc with
# "withBinaryFile: does not exist". That cost a ~2 hour SHAP render during this work. Giving each
# render its own intermediates_dir removes the collision, which is what makes concurrency safe
# here at all. Running the eight-contrast set took ~2 h at 3 concurrent against ~6.5 h serial.
#
# Concurrency default is 3, not the core count: each render holds the feature matrix and a forest
# in memory, and this machine has 8 GB. Raise it only if you have checked the headroom.
set -uo pipefail
JOBS="${1:?usage: $0 <jobfile> [max_concurrent]}"
MAXC="${2:-3}"
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"

run_pair() {
  local G="$1" C="$2" T="$3" EXTRA="$4"
  for NB in 10_rf_regression 11_shap_regression; do
    local W; W="$(mktemp -d "/tmp/rp_${G}${C}${T}_XXXX")"
    Rscript --vanilla -e "
    if (Sys.getenv('RSTUDIO_PANDOC') == '')
      Sys.setenv(RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools')
    nb <- 'code/translational_control_overview/${NB}.Rmd'
    p <- list(gene='${G}', condition='${C}', timepoint='${T}', run_cv=FALSE, save_plots=FALSE)
    ${EXTRA}
    # 11_ declares fewer params than 10_, and rmarkdown hard-errors on an undeclared one, so
    # pass only what this notebook accepts. Read the names rather than hardcoding them.
    p <- p[intersect(names(p), names(rmarkdown::yaml_front_matter(nb)\$params))]
    t0 <- Sys.time()
    ok <- tryCatch({rmarkdown::render(nb, params = p, output_file = file.path('${W}','o.html'),
          intermediates_dir = '${W}', quiet = TRUE); 'OK'},
          error = function(e) paste('FAILED:', conditionMessage(e)))
    cat(sprintf('%-3s %-9s %-8s %-20s %-9s %5.1f min\n','${G}','${C}','${T}','${NB}', ok,
        as.numeric(difftime(Sys.time(), t0, units='mins'))))" 2>&1 | grep -vE "Warning|built under"
    rm -rf "$W"
  done
}

n=0
while IFS='|' read -r spec extra; do
  [ -z "${spec// /}" ] && continue
  # shellcheck disable=SC2086
  set -- $spec
  while [ "$(jobs -rp | wc -l)" -ge "$MAXC" ]; do sleep 10; done
  run_pair "$1" "$2" "$3" "$extra" &
  n=$((n+1))
done < "$JOBS"
wait
echo "ALL $n CONTRASTS DONE"
