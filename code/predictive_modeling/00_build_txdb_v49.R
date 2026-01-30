# 00_build_txdb_v49.R
# Build a TxDb object from GENCODE v49 GTF and save as SQLite
# Only needs to run once.

library(GenomicFeatures)
library(here)

gtf_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
gtf_local <- here("accessories", "human", "gencode.v49.annotation.gtf.gz")

# Download if not already present
if (!file.exists(gtf_local)) {
  download.file(gtf_url, gtf_local, mode = "wb")
}

# Build TxDb
txdb <- makeTxDbFromGFF(gtf_local, format = "gtf",
                        organism = "Homo sapiens",
                        chrominfo = NULL)

# Save
out_path <- here("accessories", "human", "txdb.gencode49.sqlite")
saveDb(txdb, file = out_path)
cat("Saved TxDb to", out_path, "\n")
