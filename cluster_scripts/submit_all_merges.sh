#!/bin/bash
#BSUB -q rna
#BSUB -R "select[mem>30] rusage[mem=30]"
#BSUB -R "span[hosts=1]"
#BSUB -o logs/submit_all_%J.log
#BSUB -e logs/submit_all_%J.out

cd $(dirname $0)/merge_bams

bsub < merge_hypoxia_si3d_RiboSeq_1h.bsub
bsub < merge_hypoxia_si3d_RiboSeq_4h.bsub
bsub < merge_hypoxia_si3d_RiboSeq_24h.bsub
bsub < merge_hypoxia_si3e_RiboSeq_1h.bsub
bsub < merge_hypoxia_si3e_RiboSeq_4h.bsub
bsub < merge_hypoxia_si3e_RiboSeq_24h.bsub
bsub < merge_hypoxia_siCTRL_RiboSeq_1h.bsub
bsub < merge_hypoxia_siCTRL_RiboSeq_4h.bsub
bsub < merge_hypoxia_siCTRL_RiboSeq_24h.bsub
bsub < merge_normoxia_si3d_RiboSeq_1h.bsub
bsub < merge_normoxia_si3d_RiboSeq_4h.bsub
bsub < merge_normoxia_si3d_RiboSeq_24h.bsub
bsub < merge_normoxia_si3e_RiboSeq_1h.bsub
bsub < merge_normoxia_si3e_RiboSeq_4h.bsub
bsub < merge_normoxia_si3e_RiboSeq_24h.bsub
bsub < merge_normoxia_siCTRL_RiboSeq_1h.bsub
bsub < merge_normoxia_siCTRL_RiboSeq_4h.bsub
bsub < merge_normoxia_siCTRL_RiboSeq_24h.bsub
