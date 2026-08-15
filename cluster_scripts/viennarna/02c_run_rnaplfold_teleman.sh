#!/usr/bin/env bash
# RNAplfold accessibility for the Teleman (HeLa) dominant transcript set.
#
# Written for an interactive session, not sbatch. Start one first, e.g.:
#   sinteractive --ntasks=1 --mem=16G --time=8:00:00
# then:
#   bash 02c_run_rnaplfold_teleman.sh
#
# Parameters match the existing CDS/3'UTR run (-W 150 -L 100) so the resulting
# accessibility values are comparable to the MDA-MB-231 matrix.

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate working_orf

WORK_DIR=/beevol/home/matlinka/timecourse/viennarna
cd "${WORK_DIR}"

for REGION in cds utr5 utr3; do
  FA="${WORK_DIR}/teleman_${REGION}_sequences.fa"
  OUT="${WORK_DIR}/teleman_${REGION}_rnaplfold_output"

  if [[ ! -s "${FA}" ]]; then
    echo "MISSING: ${FA} - run 01c_generate_fasta_teleman.R first" >&2
    exit 1
  fi

  mkdir -p "${OUT}"
  cd "${OUT}"
  echo "[$(date +%H:%M:%S)] RNAplfold on ${REGION} ($(grep -c '^>' "${FA}") sequences)"
  RNAplfold -W 150 -L 100 < "${FA}"
  echo "[$(date +%H:%M:%S)] ${REGION} done: $(ls -1 *_lunp 2>/dev/null | wc -l) _lunp files"
  cd "${WORK_DIR}"
done

echo "All RNAplfold regions complete."
