# Kate, 2026-11-19: take the 5% of mRNAs expressed in MDA-MB-231 with the shortest CDS and ask
# which pathways they are enriched in; then the same for the longest 5%.
#
# Expressed genes = every gene in feature_matrix_tco.rds: each gene's dominant isoform in siCTRL
# normoxia MDA-MB-231 RNA-seq, the same universe as 26e_/26f_. CDS length is measured on that
# isoform. ORA method is 26e_'s: clusterProfiler::enricher on symbols, universe = all those genes,
# FDR within each collection. Collections: Hallmark, KEGG, all of C2, GO BP.
#
# Run: Rscript code/translational_control_overview/26h_shortest_cds_pathway_ora.R [shortest|longest]

suppressPackageStartupMessages({ library(tidyverse); library(here); library(msigdbr); library(clusterProfiler) })
O <- here("output", "translational_control_overview")
FRAC <- 0.05; PADJ <- 0.05
TAIL <- commandArgs(TRUE)[1]; if (is.na(TAIL)) TAIL <- "shortest"
stopifnot(TAIL %in% c("shortest", "longest"))

fm <- readRDS(file.path(O, "feature_matrix_tco.rds"))
g <- tibble(gene_id_clean = fm$gene_id_clean, symbol = fm$symbol,
            log2_cds_length = as.numeric(unlist(fm$log2_cds_length))) %>%
  filter(!is.na(symbol), !is.na(log2_cds_length)) %>% arrange(gene_id_clean) %>%
  distinct(symbol, .keep_all = TRUE)
stopifnot(all(is.finite(g$log2_cds_length)), !anyDuplicated(g$symbol))

cut <- quantile(g$log2_cds_length, if (TAIL == "shortest") FRAC else 1 - FRAC)
short <- if (TAIL == "shortest") g %>% filter(log2_cds_length <= cut) else g %>% filter(log2_cds_length >= cut)
cat(sprintf("universe: %d expressed genes | %s %.0f%%: %d genes, CDS %s %d nt (median of set %d nt, transcriptome median %d nt)\n",
            nrow(g), TAIL, 100 * FRAC, nrow(short), if (TAIL == "shortest") "<=" else ">=", round(2^cut), round(2^median(short$log2_cds_length)),
            round(2^median(g$log2_cds_length))))

c2 <- suppressWarnings(msigdbr(species = "Homo sapiens", category = "C2"))
colls <- list(
  Hallmark = suppressWarnings(msigdbr(species = "Homo sapiens", category = "H")) %>% dplyr::select(gs_name, gene_symbol),
  KEGG     = c2 %>% filter(grepl("KEGG", gs_subcat, ignore.case = TRUE)) %>% dplyr::select(gs_name, gene_symbol),
  C2       = c2 %>% dplyr::select(gs_name, gene_symbol),
  GO_BP    = suppressWarnings(msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")) %>% dplyr::select(gs_name, gene_symbol))
stopifnot(n_distinct(colls$Hallmark$gs_name) == 50)

ora <- imap_dfr(colls, function(t2g, cn) {
  r <- enricher(gene = short$symbol, universe = g$symbol, TERM2GENE = t2g, pAdjustMethod = "fdr",
                pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)
  if (is.null(r)) return(tibble())
  as_tibble(r@result) %>% mutate(collection = cn, .before = 1)
}) %>%
  mutate(k = as.numeric(sub("/.*", "", GeneRatio)), n_q = as.numeric(sub(".*/", "", GeneRatio)),
         K = as.numeric(sub("/.*", "", BgRatio)),  N = as.numeric(sub(".*/", "", BgRatio)),
         fold_enrichment = (k / n_q) / (K / N)) %>%
  dplyr::select(collection, pathway = ID, GeneRatio, BgRatio, fold_enrichment, pvalue, p.adjust, geneID)

write_csv(short, file.path(O, paste0(TAIL, "_cds_5pct_genes.csv")))
write_csv(ora, file.path(O, paste0(TAIL, "_cds_5pct_pathway_ora.csv")))

sig <- ora %>% filter(p.adjust < PADJ)
cat("\nsignificant (FDR <", PADJ, ") per collection:\n"); print(ora %>% group_by(collection) %>%
  summarise(tested = n(), significant = sum(p.adjust < PADJ)))
options(width = 200)
for (cn in names(colls)) {
  d <- sig %>% filter(collection == cn) %>% arrange(p.adjust)
  cat("\n==", cn, "(", nrow(d), "significant ) ==\n")
  if (nrow(d)) print(d %>% transmute(pathway, GeneRatio, BgRatio, FE = round(fold_enrichment, 2),
                                     padj = signif(p.adjust, 2)) %>% head(30), n = 30)
}
