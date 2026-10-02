#!/usr/bin/env bash
# End-to-end MCF7-SIX1 build: finish the baseline, score G4mer, build the matrix, fit the models,
# add the columns to the heatmap. Sequential by design - this machine has 8 GB and concurrent
# renders were being reaped under memory pressure.
#
#   bash code/translational_control_overview/01j_run_mcf7six1_pipeline.sh
#
# Stops at the first failed gate. The one deliberate exception is G4mer inference: a region that
# fails is skipped and the merge falls back to the transcripts already scored, leaving the rest
# NA to be imputed like any other missing feature. That keeps the feature SET identical across
# cell lines, which is what the heatmap depends on.
set -uo pipefail
cd "$(cd "$(dirname "$0")/../.." && pwd)"
# FROM=<n> skips earlier steps whose outputs are already on disk. G4mer inference in particular
# takes ~2.5 h and must not be repeated once its merged summaries exist.
FROM=${FROM:-1}
skip() { [ "$FROM" -gt "$1" ] && { echo "   (step $1 skipped: FROM=$FROM)"; return 0; } || return 1; }
export RSTUDIO_PANDOC="${RSTUDIO_PANDOC:-/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools}"
LOG=${LOG_DIR:-/tmp}/mcf7pipe; mkdir -p "$LOG"
R() { Rscript --vanilla "$@"; }
say() { echo; echo "=============== $* ==============="; date '+%H:%M:%S'; }
BLD=code/translational_control_overview/01i_build_mcf7six1_baseline.R
PY=/opt/homebrew/Caskroom/mambaforge/base/envs/g4mer/bin/python

say "1/7 positional CSC"
skip 1 || R $BLD csc > "$LOG/1_csc.log" 2>&1 || { echo "FAILED - see $LOG/1_csc.log"; tail -15 "$LOG/1_csc.log"; exit 1; }
[ "$FROM" -le 1 ] && grep -E "gate PASSED|cached" "$LOG/1_csc.log"

say "2/7 assemble the baseline"
skip 2 || R $BLD assemble > "$LOG/2_asm.log" 2>&1 || { echo "FAILED - see $LOG/2_asm.log"; tail -15 "$LOG/2_asm.log"; exit 1; }
[ "$FROM" -le 2 ] && grep -E "carried over|wrote feature_matrix" "$LOG/2_asm.log"

say "3/7 G4mer FASTAs for unscored transcripts"
skip 3 || R g4mer/01b_generate_fasta_mcf7six1.R > "$LOG/3_fa.log" 2>&1 || { echo "FAILED"; tail -15 "$LOG/3_fa.log"; exit 1; }
[ "$FROM" -le 3 ] && grep -E "^(utr5|cds|utr3)" "$LOG/3_fa.log"

say "4/7 G4mer inference (failures are tolerated - see the header)"
OUTD=g4mer/out_mcf7six1; mkdir -p "$OUTD"
for spec in "utr5 10" "cds 20" "utr3 20"; do
  set -- $spec; RG=$1; STR=$2
  if [ "$FROM" -gt 4 ]; then echo "  $RG: skipped (FROM=$FROM)"; continue; fi
  FA="g4mer/fasta/mcf7six1_${RG}_sequences.fa"
  if [ ! -s "$FA" ] || [ "$(grep -c '^>' "$FA")" -eq 0 ]; then
    echo "  $RG: nothing to score, skipping"; continue
  fi
  echo "  $RG: $(grep -c '^>' "$FA") sequences at stride $STR"
  if "$PY" g4mer/02_run_g4mer.py --region "$RG" --fasta "$FA" --outdir "$OUTD" \
       --stride "$STR" --device cpu --batch-size 128 > "$LOG/4_g4_$RG.log" 2>&1; then
    echo "    OK: $(tail -3 "$LOG/4_g4_$RG.log" | head -1)"
  else
    echo "    FAILED (tolerated) - see $LOG/4_g4_$RG.log"; tail -4 "$LOG/4_g4_$RG.log"
  fi
done

say "5/7 merge G4mer summaries"
if [ "$FROM" -le 5 ]; then
  R g4mer/05_merge_mcf7six1_summaries.R "$OUTD" > "$LOG/5_merge.log" 2>&1; MERGE_RC=$?
else
  echo "   (step 5 skipped: FROM=$FROM)"; MERGE_RC=0
fi
[ "$FROM" -le 5 ] && grep -E "^(utr5|cds|utr3)|WARNING|floor|wrote" "$LOG/5_merge.log"
[ $MERGE_RC -ne 0 ] && { echo "MERGE GATE FAILED - stopping before the matrix"; exit 1; }

say "6/7 build feature_matrix_tco_mcf7six1.rds"
if [ "$FROM" -gt 6 ]; then echo "   (step 6 skipped: FROM=$FROM)"; else
R -e "rmarkdown::render('code/translational_control_overview/01_build_feature_matrix.Rmd',
  params=list(cell_line='MCF7-SIX1',
              baseline_rds='feature_matrix_mcf7six1_dhx29_kd_feature.rds',
              kozak_seqs_rds='precomputed_mcf7six1_kozak_sequences.rds',
              dhx29_occ_rds='dhx29_riboseq_occupancy_mcf7six1.rds',
              g4mer_tag='_mcf7six1', cache_tag='_mcf7six1',
              codon_weights_rds='codon_weights_mdamb231.rds',
              out_rds='feature_matrix_tco_mcf7six1.rds'),
  output_file='01_build_feature_matrix_mcf7six1.html',
  intermediates_dir='$LOG/km01m', quiet=TRUE)" > "$LOG/6_matrix.log" 2>&1 \
  || { echo "FAILED - see $LOG/6_matrix.log"; tail -20 "$LOG/6_matrix.log"; exit 1; }
grep -E "annotation parity|Wrote|Cell line" "$LOG/6_matrix.log" | head
fi

say "7/7 fit the 4 contrasts x 2 feature sets, then the heatmaps"
EX="is_dap5_repressed,is_dap5_promoted,dap5_polysome_lfc,clip_eif3_count,cnot3_ribo_enrichment,cnot3_weighted_codon_score,cnot3_weighted_codon_score_tx,zhu2024_ctrl_halflife,zhu2024_cnot3ko_halflife_logfc,karner2026_mda231_log2_ct,hia2026_dhx29_kd_rna_lfc,hia2026_dhx29_occupancy_tx,clip_cap_binding_count_v2"
RF=eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg_corrected_both_clipaggregate_yraw
RI=eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg_intrinsic_corrected_both_clipaggregate_yraw
M="p\$matrix_rds <- 'feature_matrix_tco_mcf7six1.rds'"
J=code/translational_control_overview/jobs_mcf7six1_owntx.txt
: > "$J"
for g in 3d 3e; do for c in hypoxia normoxia; do
  X="p\$cell_line <- 'MCF7-SIX1'; p\$xcell_ref_csv <- 'translation_categories_si${g}_vs_sictrl_${c}_1hr.csv'; $M"
  # A DISTINCT feature_set_label: the existing tcoreg runs are the same contrasts fitted on the
  # 231-isoform matrix and are the control arm for the isoform effect. Reusing the label would
  # overwrite them.
  echo "$g mcf7six1_${c} 1hr | $X; p\$feature_set_label <- 'tcoreg_mcf7tx'; p\$parity_ref_suffix <- '$RF'" >> "$J"
  echo "$g mcf7six1_${c} 1hr | $X; p\$feature_set_label <- 'tcoreg_mcf7tx_intrinsic'; p\$parity_ref_suffix <- '$RI'; p\$exclude_features <- '$EX'" >> "$J"
done; done
bash code/translational_control_overview/10p_parallel_runner.sh "$J" 1 2>&1 | tee "$LOG/7_models.log" | grep -E "OK|FAILED|DONE"

say "8/8 heatmaps on MCF7-SIX1's own dominant transcripts"
# The columns cannot simply be appended to the existing 231 figure: 15_ asserts every column
# shares one GENE set, and these are built on the genes MCF7-SIX1 can score. A matched
# cross-cell-line figure needs the 231 runs restricted to the same genes, which is a separate
# batch (jobs_xcell_231.txt).
S=_lfc0.5_tcoreg_mcf7tx; T=_corrected_both_clipaggregate_yraw
EF="MCF7-SIX1 si3d hypoxia 1hr=eif33d_promotes_mcf7six1_hypoxia_1hr${S}${T};MCF7-SIX1 si3d normoxia 1hr=eif33d_promotes_mcf7six1_normoxia_1hr${S}${T};MCF7-SIX1 si3e hypoxia 1hr=eif33e_promotes_mcf7six1_hypoxia_1hr${S}${T};MCF7-SIX1 si3e normoxia 1hr=eif33e_promotes_mcf7six1_normoxia_1hr${S}${T}"
EI=$(echo "$EF" | sed "s/_tcoreg_mcf7tx_corrected/_tcoreg_mcf7tx_intrinsic_corrected/g")
CL=eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg_corrected_both_clipaggregate_yraw
for V in full intrinsic; do
  if [ "$V" = full ]; then E="$EF"; TG=mcf7six1_owntx; else E="$EI"; TG=mcf7six1_owntx_intrinsic; fi
  R -e "rmarkdown::render('code/translational_control_overview/15_shap_heatmap.Rmd',
    params=list(experiments='$E', clusters_suffix='$CL', reduced_suffix='', out_tag='$TG'),
    output_file='15_shap_heatmap_${TG}.html', intermediates_dir='$LOG/hm_$TG', quiet=TRUE)"     > "$LOG/8_hm_$TG.log" 2>&1 && echo "  $V heatmap OK"     || { echo "  $V heatmap FAILED"; tail -12 "$LOG/8_hm_$TG.log"; }
done

say "PIPELINE COMPLETE"
echo "logs in $LOG"
