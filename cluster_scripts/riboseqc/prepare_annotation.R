library(RiboseQC)

prepare_annotation_files(
  annotation_directory = "/beevol/home/matlinka/eif3e_hypoxia/ORFquant/orfquant_hypoxia_timecourse/riboseqc",
  twobit_file = "/beevol/home/matlinka/eif3e_hypoxia/ORFquant/orfquant_hypoxia_timecourse/riboseqc/GRCh38.primary_assembly.genome.2bit",
  gtf_file = "/beevol/home/bartholk/riboPipe_comps/star_index/v45/gencode.v45.annotation.gtf.gz",
  genome_seq = "/beevol/home/bartholk/riboPipe_comps/star_index/v45/GRCh38.primary_assembly.genome.fa",
  scientific_name = "Homo.sapiens",
  annotation_name = "gencode.v45",
  export_bed_tables_TxDb = FALSE,
  forge_BSgenome = FALSE,
  create_TxDb = TRUE
)
