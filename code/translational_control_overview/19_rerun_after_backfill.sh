#!/bin/bash
# 19_rerun_after_backfill.sh - rebuild every matrix and rerun the regression heatmap pipeline
# after the 2026-10-01 MDA-MB-231 backfill (G4 + DHX29 occupancy for the 1,600 transcripts the
# older precomputed_most_abundant_tx.rds list never covered, per-transcript CNOT3 codon score,
# late structure files). Scope "A" (Kate, 2026-10-01): the regression pipeline that feeds the
# heatmaps; heatmaps BEFORE the 20-refit noise floor. Classifier sweeps are NOT rerun.
#
#   bash code/translational_control_overview/19_rerun_after_backfill.sh [from_phase]
# Phases: 1 MDA-MB-231 matrix | 2 other matrices | 3 classifier + outcome prep | 4 models |
#         5 clusters + heatmaps | 6 noise floor. Stops at the first failure.
set -uo pipefail
FROM="${1:-1}"
cd "$(dirname "$0")/../.." || exit 1
export RSTUDIO_PANDOC="${RSTUDIO_PANDOC:-/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools}"
TCO=code/translational_control_overview
LOG="${LOG_DIR:-/tmp}/rerun19"; mkdir -p "$LOG"
say()  { printf '\n===== %s =====  %s\n' "$1" "$(date +%H:%M:%S)"; }
step() { local n="$1"; shift; "$@" > "$LOG/$n.log" 2>&1; local rc=$?
         if [ $rc -ne 0 ]; then echo "FAILED at $n (exit $rc) - $LOG/$n.log"; tail -15 "$LOG/$n.log"; exit 1; fi
         echo "ok: $n"; }
render() { Rscript --vanilla -e "rmarkdown::render('$TCO/$1.Rmd', params = list($2), output_file = '$3', intermediates_dir = tempfile(), quiet = TRUE)"; }
F231="eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg_corrected_both_clipaggregate_yraw"
I231="eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg_intrinsic_corrected_both_clipaggregate_yraw"

if [ "$FROM" -le 1 ]; then
  say "1 MDA-MB-231 matrix"
  step 231_late_struct Rscript --vanilla $TCO/01l_build_celltx_baseline.R mdamb231 struct
  step 231_matrix render 01_build_feature_matrix "g4mer_tag = '_mdamb231', dhx29_occ_rds = 'dhx29_riboseq_occupancy_mdamb231.rds', late_struct_rds = 'mdamb231_struct_from_lunp.rds', check_golden = FALSE" 01_build_feature_matrix.html
  grep -E "late structure files|cnot3 codon score|GATE|values that moved" "$LOG/231_matrix.log"
  # Only the intended columns may differ from the pre-backfill matrix: values may be FILLED in the
  # backfilled families, and G4 residuals may shift (their coefficients are refit); nothing else.
  step 231_diff Rscript --vanilla $TCO/19b_backfill_diff_check.R
  cat "$LOG/231_diff.log" | grep -v -E "Warn|built under"
  step golden Rscript --vanilla $TCO/01d_golden_transcripts.R
  step dictionary Rscript --vanilla $TCO/01b_feature_dictionary.R
fi

if [ "$FROM" -le 2 ]; then
  say "2 MCF7-SIX1, HEK293T, HeLa, MANE matrices (shared builder, value parity vs MDA-MB-231)"
  for c in mcf7six1 dhx29 hela mane; do
    step "matrix_$c" bash $TCO/01m_run_celltx_build.sh "$c"
    grep -E "PARITY PASSED|identical:" "$LOG/matrix_$c.log" | tail -2
  done
fi

if [ "$FROM" -le 3 ]; then
  say "3 classifier for the feature-parity gate, Teleman outcome prep"
  step classifier_3d_hyp1 render 02_rf_model "save_plots = FALSE" 02_rf_model.html
  step teleman_prep Rscript --vanilla $TCO/10t_teleman_outcome_prep.R
fi

if [ "$FROM" -le 4 ]; then
  say "4 regression + SHAP (MDA-MB-231 first: the other cell sets assert parity against it)"
  for j in jobs_231_full jobs_intrinsic_example jobs_mcf7six1_owntx jobs_hela_owntx jobs_hek293t_owntx jobs_dap5_mane; do
    bash $TCO/10p_parallel_runner.sh $TCO/$j.txt 3 > "$LOG/models_$j.log" 2>&1
    cat "$LOG/models_$j.log"
    if grep -q "FAILED" "$LOG/models_$j.log"; then echo "FAILED in $j - stopping"; exit 1; fi
  done
fi

if [ "$FROM" -le 5 ]; then
  say "5 correlation clusters, then heatmaps"
  step clusters_full render 14_feature_correlation "rf_suffix = '$F231'" 14_feature_correlation.html
  step clusters_intr render 14_feature_correlation "rf_suffix = '$I231'" 14_feature_correlation_intrinsic.html
  step heatmap_full Rscript --vanilla $TCO/15b_render_all_celllines.R full
  step heatmap_intr Rscript --vanilla $TCO/15b_render_all_celllines.R intrinsic
  step heatmap_capbind Rscript --vanilla $TCO/15c_render_subunit_seq.R
  say "HEATMAPS DONE"
fi

if [ "$FROM" -le 6 ]; then
  say "6 noise floor: 20 split-seed refits of MDA-MB-231 si3d hypoxia 1hr"
  step floor_sweep Rscript --vanilla $TCO/11d_shap_regression_seed_sweep.R transform raw 20
  step floor_report render 11f_shap_regression_stability "" 11f_shap_regression_stability.html
fi
say "ALL DONE (logs: $LOG)"
