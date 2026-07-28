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
# and small enough for a laptop. CDS + 3'UTR are ~31 Mb of sequence, about 2.8M
# windows at stride 10, which is what wants a GPU node.
#
# Prerequisites on the cluster, once (run these on the LOGIN node):
#   conda create -y -n g4mer python=3.10
#   conda activate g4mer
#   pip install torch transformers
#   hf auth login          # Biociphers/g4mer is a gated repo
#   hf download Biociphers/g4mer   # pre-cache: compute nodes often lack internet
#
# FASTA is generated locally and copied over — no need for R/BSgenome here:
#   scp g4mer/fasta/{cds,utr3}_sequences.fa <cluster>:.../g4mer/fasta/
#   scp g4mer/02_run_g4mer.py <cluster>:.../g4mer/
#
# Submit: sbatch 03_run_g4mer_cluster.sbatch

source $(conda info --base)/etc/profile.d/conda.sh
conda activate g4mer

WORK_DIR=/beevol/home/matlinka/timecourse/g4mer
OUT_DIR=${WORK_DIR}/output
mkdir -p ${OUT_DIR} ${WORK_DIR}/logs

# Keep the model cache on shared storage rather than a quota-limited home dir.
export HF_HOME=${WORK_DIR}/hf_cache
# Model is pre-cached on the login node, so no network call is needed here.
# Comment this out if the compute node does have outbound internet.
export HF_HUB_OFFLINE=1

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
    --stride 20 \
    --selftest
done

echo "All regions complete. Copy back to the repo:"
echo "  scp ${OUT_DIR}/g4mer_summary_*.tsv    <local>/output/g4mer/"
echo "  scp ${OUT_DIR}/g4mer_windows_*.tsv.gz <local>/output/g4mer/"
