library(devtools)
library("RiboseQC")

RiboseQC_analysis(
  annotation_file = "/beevol/home/matlinka/eif3e_hypoxia/ORFquant/orfquant_hypoxia_timecourse/riboseqc/gencode.v45.annotation.gtf_Rannot",
  bam_files = "/beevol/home/matlinka/eif3e_hypoxia/ORFquant/orfquant_hypoxia_timecourse/merge_bams/merged_bams/normoxia_si3e_RiboSeq_4h_merged.bam",
  report_file = "normoxia_si3e_RiboSeq_4h_riboseqc.html",
  write_tmp_files = FALSE
)
