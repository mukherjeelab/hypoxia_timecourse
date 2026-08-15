#!/usr/bin/env bash
# RNAfold MFE for the Teleman (HeLa) dominant transcript set.
#
# Interactive session, not sbatch:
#   sinteractive --ntasks=1 --mem=16G --time=8:00:00
#   bash 03c_run_rnafold_teleman.sh
#
# --noPS suppresses per-sequence postscript plots, which would otherwise write one
# file per transcript into the working directory.

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate working_orf

WORK_DIR=/beevol/home/matlinka/timecourse/viennarna
cd "${WORK_DIR}"

for REGION in cds utr5 utr3; do
  FA="${WORK_DIR}/teleman_${REGION}_sequences.fa"
  OUT="${WORK_DIR}/teleman_${REGION}_rnafold.txt"

  if [[ ! -s "${FA}" ]]; then
    echo "MISSING: ${FA} - run 01c_generate_fasta_teleman.R first" >&2
    exit 1
  fi

  echo "[$(date +%H:%M:%S)] RNAfold on ${REGION} ($(grep -c '^>' "${FA}") sequences)"
  RNAfold --noPS < "${FA}" > "${OUT}"
  echo "[$(date +%H:%M:%S)] ${REGION} done -> $(basename "${OUT}") ($(wc -l < "${OUT}") lines)"
done

echo "All RNAfold regions complete."
