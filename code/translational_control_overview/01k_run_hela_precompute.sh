#!/usr/bin/env bash
# Per-cell-set precompute (CELL=hela default | dhx29 | mcf7six1): RNAplfold + G4mer for the transcripts that no
# existing run has covered. Mirrors steps 3-5 of 01j_run_mcf7six1_pipeline.sh plus the RNAplfold
# steps that ran separately for MCF7-SIX1. Sequential: this machine has 8 GB.
#   bash code/translational_control_overview/01k_run_hela_precompute.sh [cell]
set -uo pipefail
cd "$(cd "$(dirname "$0")/../.." && pwd)"
CELL="${1:-hela}"
PY=/opt/homebrew/Caskroom/mambaforge/base/envs/g4mer/bin/python
say() { echo; echo "=============== $* ==============="; date '+%H:%M:%S'; }
say "1/5 RNAplfold FASTAs (missing _lunp only)"
Rscript --vanilla cluster_scripts/viennarna/01d_generate_fasta_mcf7six1.R "$CELL" || exit 1
say "2/5 RNAplfold, 4 workers"
bash cluster_scripts/viennarna/02d_run_rnaplfold_mcf7six1.sh 4 "$CELL" || exit 1
say "3/5 G4mer FASTAs (unscored only)"
Rscript --vanilla g4mer/01b_generate_fasta_mcf7six1.R "$CELL" || exit 1
say "4/5 G4mer inference on CPU (failures tolerated, as for MCF7-SIX1)"
OUTD=g4mer/out_$CELL; mkdir -p "$OUTD"
for spec in "utr5 10" "cds 20" "utr3 20"; do
  set -- $spec; RG=$1; STR=$2; FA="g4mer/fasta/${CELL}_${RG}_sequences.fa"
  if [ ! -s "$FA" ] || [ "$(grep -c '^>' "$FA")" -eq 0 ]; then echo "  $RG: nothing to score"; continue; fi
  echo "  $RG: $(grep -c '^>' "$FA") sequences at stride $STR"
  "$PY" g4mer/02_run_g4mer.py --region "$RG" --fasta "$FA" --outdir "$OUTD" --stride "$STR" \
    --device cpu --batch-size 128 && echo "    OK $RG" || echo "    FAILED (tolerated) $RG"
done
say "5/5 merge G4mer summaries"
Rscript --vanilla g4mer/05_merge_mcf7six1_summaries.R "$OUTD" "$CELL"
say "$CELL precompute done"
