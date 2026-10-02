# 01q_mane_tx_map.R
# The MANE Select transcript set (cell set "mane") shared by every HEK293T column that is fitted on
# MANE rather than on a transcript chosen from its own data (Kate, 2026-10-01):
#   - Weber 2022 DAP5 KO TE (10w_weber_dap5_outcome_prep.Rmd): no transcript-level data exists;
#   - subunit-seq eIF3d cap binding (10s_subunit_seq_outcome_prep.Rmd): MANE judged a better fit
#     than borrowing the DHX29 experiment's footprint-chosen HEK293T transcripts.
# One map for both, so each gene has one transcript and one feature row however many outcomes use
# it. MANE Select is one transcript per gene, so taking the union cannot create a conflict.
#   Rscript code/translational_control_overview/01q_mane_tx_map.R
suppressPackageStartupMessages({ library(dplyr); library(readr); library(here) })
PD <- here("output", "predictive_modeling")
mane <- read_csv(here("accessories", "human", "mane_select_v49.csv"),
                 col_names = c("gene_id_clean", "transcript_id_clean", "symbol"), show_col_types = FALSE)
stopifnot("MANE file is not one transcript per gene" = !any(duplicated(mane$gene_id_clean)),
          "MANE file is not one gene per transcript" = !any(duplicated(mane$transcript_id_clean)))

# DAP5: 10w joined by symbol (direct matches only) and wrote the MANE transcript it chose.
dap5 <- read_csv(here("output", "translation_categories_sidap5_vs_sictrl_hek293t_normoxia_steadystate.csv"),
                 show_col_types = FALSE)
stopifnot("10w's DAP5 transcripts are not MANE Select for their gene" =
            all(paste(dap5$gene_id, dap5$transcript_id_clean) %in%
                paste(mane$gene_id_clean, mane$transcript_id_clean)))
# Cap binding: gene-keyed outcome files; a gene with no MANE Select entry cannot be placed.
cb_files <- list.files(here("output"), "^subunit_seq_hek293t_capbind_.*_rpkm[.]csv$", full.names = TRUE)
stopifnot("expected the five cap-binding outcome files" = length(cb_files) == 5)
cb <- unique(unlist(lapply(cb_files, function(f) sub("[.].*", "", read_csv(f, show_col_types = FALSE)$gene_id))))

genes <- union(dap5$gene_id, cb)
map <- mane %>% filter(gene_id_clean %in% genes)
cat(sprintf("DAP5 genes %d | cap-binding genes %d, of which with a MANE Select %d (%d without) | union mapped %d\n",
            nrow(dap5), length(cb), sum(cb %in% mane$gene_id_clean), sum(!cb %in% mane$gene_id_clean), nrow(map)))
stopifnot("a DAP5 gene fell out of the map" = all(dap5$gene_id %in% map$gene_id_clean))
saveRDS(map, file.path(PD, "precomputed_mane_tx.rds"))
cat("wrote precomputed_mane_tx.rds\n")
