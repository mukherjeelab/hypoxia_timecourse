#!/usr/bin/env bash
#SBATCH --job-name=g4mer_cds_utr3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=08:00:00
#SBATCH --output=/beevol/home/matlinka/timecourse/g4mer/logs/g4mer_cds_utr3_%j.log
#SBATCH --error=/beevol/home/matlinka/timecourse/g4mer/logs/g4mer_cds_utr3_%j.err

# G4mer rG4 scoring for CDS and 3'UTR.
#
# The 5'UTR run is done locally (see g4mer/README.md) — it is the primary analysis
# and small enough for an M2. CDS + 3'UTR are ~31 Mb of sequence, roughly 3M
# windows at stride 10, which is what needs a GPU node.
#
# Prerequisites on the cluster, once:
#   mamba create -n g4mer python=3.10
#   mamba activate g4mer
#   pip install torch transformers
#   Rscript 01_generate_fasta.R --workdir /beevol/home/matlinka/timecourse/g4mer
#
# Submit: sbatch 03_run_g4mer_cluster.sbatch

source $(conda info --base)/etc/profile.d/conda.sh
conda activate g4mer

WORK_DIR=/beevol/home/matlinka/timecourse/g4mer
OUT_DIR=${WORK_DIR}/output
mkdir -p ${OUT_DIR} ${WORK_DIR}/logs

# Fail loudly rather than writing a truncated score table downstream.
set -euo pipefail

nvidia-smi || echo "WARNING: no GPU visible, falling back to CPU"

for REGION in cds utr3; do
  echo "=== ${REGION} ==="
  python ${WORK_DIR}/02_run_g4mer.py \
    --region ${REGION} \
    --basedir ${WORK_DIR} \
    --outdir ${OUT_DIR} \
    --device auto \
    --batch-size 512 \
    --stride 10 \
    --selftest
done

echo "All regions complete. Copy back to the repo:"
echo "  scp ${OUT_DIR}/g4mer_summary_*.tsv    <local>/output/g4mer/"
echo "  scp ${OUT_DIR}/g4mer_windows_*.tsv.gz <local>/output/g4mer/"
