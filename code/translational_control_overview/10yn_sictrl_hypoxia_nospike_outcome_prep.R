# 10yn_sictrl_hypoxia_nospike_outcome_prep.R
# MDA-MB-231 siCTRL hypoxia vs normoxia TE change at 1hr WITHOUT spike-in normalization, as a
# regression outcome - the no-spike arm of the `spikein_ctrl` preset (10_regression_presets.R),
# whose spike-in arm is 10y_'s sictrl_hypoxia_vs_normoxia_te_1hr.csv (requested by Kate,
# 2026-11-24).
#
# Source: output/deseq_results_1hr_nospike.csv, written by code/timecourse_qc.Rmd (deFunction,
# mode "TE", sictrl, control_condition = "normoxia", test_condition = "hypoxia",
# spikeInMatrix = NULL, size factors from createNormMatrix(spikein = NULL)). Nothing is refitted.
# Orientation: log2(hypoxia / normoxia), positive = higher TE in hypoxia - the same as 10y_.
#
# The contrast identity is checked through the file's SPIKE-IN sibling: timecourse_qc.Rmd writes
# deseq_results_1hr_spike.csv from the same counts and samples with spikeInMatrix = spikein_1hr,
# and that file must reproduce 10y_'s outcome (which is itself pinned to the matrix te_lfc, F7).
# If it does, the no-spike file differs from the spike-in outcome by the normalization alone.
#
#   Rscript code/translational_control_overview/10yn_sictrl_hypoxia_nospike_outcome_prep.R
suppressPackageStartupMessages({ library(dplyr); library(readr); library(here) })
id_clean <- function(x) sub("[.].*", "", x)

spike_out <- read_csv(here("output", "sictrl_hypoxia_vs_normoxia_te_1hr.csv"), show_col_types = FALSE)
qc_spike  <- read_csv(here("output", "deseq_results_1hr_spike.csv"),   show_col_types = FALSE)
qc_nosp   <- read_csv(here("output", "deseq_results_1hr_nospike.csv"), show_col_types = FALSE)

stopifnot("duplicate gene_id in the no-spike file" = !anyDuplicated(qc_nosp$gene_id),
          "the two timecourse_qc files do not share one gene universe" =
            identical(qc_spike$gene_id, qc_nosp$gene_id))

# 1. The spike-in sibling reproduces the spike-in outcome -> same contrast, same orientation.
#    Tolerance, not ==: the two files are separate write_csv round trips (max dev ~1e-17).
m   <- match(spike_out$gene_id, qc_spike$gene_id)
dev <- max(abs(qc_spike$log2FoldChange[m] - spike_out$te_lfc), na.rm = TRUE)
cat(sprintf("timecourse_qc spike-in TE vs 10y_ outcome: %d/%d genes matched, max deviation %g\n",
            sum(!is.na(m)), nrow(spike_out), dev))
stopifnot("a 10y_ outcome gene is missing from timecourse_qc's spike-in file" = !anyNA(m),
          "timecourse_qc's spike-in TE does not reproduce the 10y_ outcome - not the same contrast" =
            dev < 1e-9)

out <- qc_nosp %>% filter(!is.na(log2FoldChange)) %>%
  transmute(gene_id, symbol, te_lfc = log2FoldChange, te_padj = padj)
stopifnot("non-finite te_lfc" = all(is.finite(out$te_lfc)),
          "the no-spike outcome does not cover exactly the spike-in outcome's genes" =
            setequal(id_clean(out$gene_id), id_clean(spike_out$gene_id)))

# 2. Same response, different normalization: strongly concordant but not identical. A file that
#    equalled the spike-in outcome would mean the wrong file was read.
j   <- inner_join(spike_out, out, by = "gene_id", suffix = c("_spike", "_nospike"))
rho <- cor(j$te_lfc_spike, j$te_lfc_nospike, method = "spearman")
cat(sprintf("spike vs no-spike te_lfc: Spearman %.4f, max|diff| %.3f (n=%d)\n",
            rho, max(abs(j$te_lfc_spike - j$te_lfc_nospike)), nrow(j)))
cat(sprintf("median te_lfc: spike %+.4f, no-spike %+.4f | padj<0.05: spike %d, no-spike %d\n",
            median(j$te_lfc_spike), median(j$te_lfc_nospike),
            sum(j$te_padj_spike < 0.05, na.rm = TRUE), sum(j$te_padj_nospike < 0.05, na.rm = TRUE)))
stopifnot("no-spike TE is not concordant with spike-in TE - wrong contrast?" = rho > 0.9,
          "no-spike TE is identical to spike-in TE - the no-spike file was not read" =
            max(abs(j$te_lfc_spike - j$te_lfc_nospike)) > 0.01)

write_csv(out, here("output", "sictrl_hypoxia_vs_normoxia_te_1hr_nospike.csv"))
cat("wrote output/sictrl_hypoxia_vs_normoxia_te_1hr_nospike.csv:", nrow(out), "genes\n")
