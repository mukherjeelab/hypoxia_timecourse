# 10x_dhx29_ip_outcome_prep.R
# DHX29 selective ribosome profiling (Hia et al., GSE266018, Thor-Ribo-seq, HEK293T) as a
# regression outcome: DHX29-IP over Input enrichment per gene, from
# code/dhx29_reanalysis/04_dhx29_enrichment_deseq.Rmd.
#
# THE OUTCOME (chosen by Kate, 2026-09-30): gene-level DESeq2 on the representative transcript's
# counts (one transcript per gene, rule set in dhx29_reanalysis 03), UNPAIRED design (~ assay),
# apeglm-SHRUNKEN log2FC (lfc_simple). Joined to the feature matrix by GENE, like every TE column.
#
# Changed from the unshrunk lfc_raw_simple on 2026-10-01 (Kate: "it should be shrunken"). Every
# other DESeq2 column is shrunken, and the unshrunk top ranks were low-count noise: the top 400 by
# unshrunk LFC had median baseMean ~10 against 25 for the rest, and their GC3 (~63%) inverted the
# GC3-poor enrichment notebook 05 reports on the shrunken estimate. padj is the same either way.
#
# CAVEATS carried to every result (dhx29_reanalysis 04):
#   * weak per-feature signal: 39 of 10,493 genes significant at FDR 0.1 (unpaired); the
#     transcript-level p-value histogram is depleted near 0 (conservative); IP replicates agree
#     with each other less than with Input
#   * no background IP: enrichment tracks length (+0.24), CDS GC (-0.22), expression (-0.18)
#   * gene level: the representative transcript differs from the matrix transcript for ~41% of
#     shared genes, as for HeLa / MCF7-SIX1 on the 231 matrix
#   * DHX29's own mRNA is the bait artefact (co-translational IP of the nascent protein); it is
#     dropped from the outcome, as EIF3D is dropped from the si3d RNA gene sets
#   Rscript code/translational_control_overview/10x_dhx29_ip_outcome_prep.R
suppressPackageStartupMessages({ library(dplyr); library(readr); library(here) })
src <- here("output", "dhx29_reanalysis", "dhx29_enrichment_gene.csv")
g <- read_csv(src, show_col_types = FALSE)
stopifnot("gene table is not the 10,493 genes notebook 04 asserts" = nrow(g) == 10493,
          "a gene appears twice" = !any(duplicated(g$gene_id_clean)),
          "shrunken unpaired LFC is missing or non-finite" =
            "lfc_simple" %in% names(g) && all(is.finite(g$lfc_simple)))
bait <- g %>% filter(symbol == "DHX29")
cat(sprintf("DHX29 self-enrichment (shrunken, unpaired): log2FC %+.2f, padj %.2g\n",
            bait$lfc_simple, bait$padj_simple))
stopifnot("the bait is not enriched in its own pull-down - sign is wrong" =
            nrow(bait) == 1 && bait$lfc_simple > 1 && bait$padj_simple < 0.01)
# Tested genes only (Kate, 2026-10-01), the rule already applied to Herrmannova: DESeq2's
# independent filtering left 2,241 genes with padj NA (median baseMean 8.4 vs 30.9 tested), and
# they crowded the top of the enrichment ranking with barely-expressed, GC3-rich genes.
n_untested <- sum(is.na(g$padj_simple) & g$symbol != "DHX29")
out <- g %>% filter(symbol != "DHX29", !is.na(padj_simple)) %>%
  transmute(gene_id = gene_id_clean, symbol, representative_transcript = transcript_id_clean,
            enrich_lfc = lfc_simple, enrich_padj = padj_simple, base_mean)
fm <- readRDS(here("output", "translational_control_overview", "feature_matrix_tco.rds"))
inm <- out %>% filter(gene_id %in% fm$gene_id_clean)
same_tx <- sum(inm$representative_transcript ==
               fm$transcript_id_clean[match(inm$gene_id, fm$gene_id_clean)])
cat(sprintf("dropped %d untested genes (padj NA)\n", n_untested))
cat(sprintf("outcome genes %d | in the matrix %d | same transcript as the matrix %d (%.0f%%)\n",
            nrow(out), nrow(inm), same_tx, 100 * same_tx / nrow(inm)))
write_csv(out, here("output", "dhx29_ip_hek293t_enrichment.csv"))

# Gene sets for 10_'s spot-check: enriched = padj < 0.1 & log2FC > 0 (unpaired); controls =
# padj > 0.3 & |log2FC| < 0.5. Few positives, so 10_'s held-out check will skip itself (< 20).
cond <- "hek293t_dhx29ip"
pos <- out %>% filter(!is.na(enrich_padj), enrich_padj < 0.1, enrich_lfc > 0)
neg <- out %>% filter(!is.na(enrich_padj), enrich_padj > 0.3, abs(enrich_lfc) < 0.5)
write_csv(pos %>% transmute(db_gene_symbol = symbol, ensembl_gene = gene_id, enrich_lfc,
                            gs_name = paste0(cond, "_enriched_fdr0.1"),
                            gs_description = "DHX29-IP enriched genes, unpaired DESeq2, FDR 0.1; not a TE set"),
          here("output", "genesets", paste0(cond, "_dhx29_promotes_TE_ribo_lfc0.5.csv")))
write_csv(neg %>% dplyr::select(gene_id, symbol, enrich_lfc),
          here("output", "genesets", paste0(cond, "_dhx29_negative_controls_lfc0.5.csv")))
cat(sprintf("gene sets: %d enriched (FDR 0.1), %d controls\n", nrow(pos), nrow(neg)))
stopifnot("enriched set does not score above controls - join or sign is wrong" =
            median(pos$enrich_lfc) - median(neg$enrich_lfc) > 0.3)
