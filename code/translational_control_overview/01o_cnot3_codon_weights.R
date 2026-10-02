# 01o_cnot3_codon_weights.R
# The three P-site weights behind cnot3_weighted_codon_score, derived from the Zhu 2024 source
# file exactly as code/predictive_modeling/22_cnot3_riboseq_slamseq.Rmd does (chunks
# riboseq_codon_enrichment and compute_codon_score), so 01_ can score EACH matrix's own
# transcript instead of carrying notebook 22's per-gene value.
#
# Why (Kate, 2026-10-01): notebook 22 scores the CDS of the transcript in
# precomputed_most_abundant_tx.rds (Salmon abundance summed over all 108 libraries, ribosome
# footprints included) and joins it by GENE. The matrix uses siCTRL normoxia RNA-seq, and the two
# lists disagree on 1,600 MDA-MB-231 transcripts, which therefore carry another isoform's score or
# none; every other cell set inherits the same old-list value.
#
# Nothing is transcribed: the weights are recomputed from the xlsx, then proven by reproducing
# notebook 22's saved scores on every transcript the two lists share.
#   Rscript code/translational_control_overview/01o_cnot3_codon_weights.R
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readxl); library(here)
  library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38)
})
PD  <- here("output", "predictive_modeling")
OUT <- here("output", "translational_control_overview")

ribo_raw <- read_xlsx(here("accessories", "csc_data", "zhu_slamseq",
                           "GSE268324_CNOT3_selective_ribosome_profiling_All_footprints.xlsx"),
                      sheet = "Footprints") %>%
  mutate(across(c(A_site_codon, P_site_codon, E_site_codon), as.character))
ribo_codon <- ribo_raw %>%
  filter((`Input-1` + `Input-2`) > 0) %>%
  mutate(enrichment = (`CNOT3-IP-1` + `CNOT3-IP-2`) / (`Input-1` + `Input-2`)) %>%
  pivot_longer(cols = c(P_site_codon, A_site_codon, E_site_codon),
               names_to = "site", values_to = "codon") %>%
  mutate(site = sub("_codon$", "", site)) %>%
  group_by(site, codon) %>%
  summarise(mean_enrichment = mean(enrichment, na.rm = TRUE), .groups = "drop")
w <- ribo_codon %>% filter(site == "P_site", codon %in% c("CGG", "CGA", "AGG"))
weights <- setNames(w$mean_enrichment, w$codon)
stopifnot("did not recover exactly the three CNOT3 codons" = setequal(names(weights), c("CGG", "CGA", "AGG")),
          # notebook 22's Figure 1F replication: these are the top three P-site codons
          "CGG/CGA/AGG are not the top three P-site codons - wrong column or site" =
            setequal(names(weights), ribo_codon %>% filter(site == "P_site") %>%
                       slice_max(mean_enrichment, n = 3) %>% pull(codon)))
print(round(weights, 4))

# notebook 22's score, verbatim: CDS from cdsBy (stop codon included), frame-checked, >= 10 codons
score_cds <- function(s) {
  n <- nchar(s) / 3
  if (n < 10) return(NA_real_)
  st <- seq(1, nchar(s) - 2, by = 3)
  sum(weights[substring(s, st, st + 2)], na.rm = TRUE) / n
}

# Proof: reproduce notebook 22's saved per-gene scores on the old-list transcripts.
txdb <- loadDb(here("accessories", "human", "txdb.gencode49.sqlite"))
old  <- readRDS(file.path(PD, "precomputed_most_abundant_tx.rds"))
saved <- readRDS(file.path(PD, "feature_matrix_cnot3.rds")) %>%
  transmute(gene_id_clean = sub("[.].*", "", as.character(unlist(gene_id_clean))),
            saved = as.numeric(unlist(cnot3_weighted_codon_score))) %>%
  filter(!is.na(saved)) %>% distinct(gene_id_clean, .keep_all = TRUE)
cds <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38, cdsBy(txdb, by = "tx", use.names = TRUE))
names(cds) <- sub("[.].*", "", names(cds))
ot <- old %>% transmute(transcript_id_clean = as.character(unlist(transcript_id_clean)),
                        gene_id_clean = as.character(unlist(gene_id_clean)))
s <- cds[names(cds) %in% ot$transcript_id_clean]; s <- s[width(s) %% 3 == 0]
mine <- tibble(transcript_id_clean = names(s), mine = vapply(as.character(s), score_cds, numeric(1))) %>%
  inner_join(ot, by = "transcript_id_clean") %>% inner_join(saved, by = "gene_id_clean")
dev <- max(abs(mine$mine - mine$saved), na.rm = TRUE)
cat(sprintf("reproduced notebook 22's saved scores on %d genes: max |dev| = %.2e\n", nrow(mine), dev))
stopifnot("the recomputed weights do not reproduce notebook 22's saved scores" =
            nrow(mine) > 8000 && dev < 1e-12)

saveRDS(list(weights = weights, min_codons = 10,
             source = "GSE268324 All_footprints.xlsx, P-site mean (IP1+IP2)/(Input1+Input2)",
             verified_genes = nrow(mine), verified_max_dev = dev),
        file.path(OUT, "cnot3_psite_weights.rds"))
cat("wrote cnot3_psite_weights.rds\n")
