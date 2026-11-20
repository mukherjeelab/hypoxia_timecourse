# One-row-per-gene summary of the reporter controls picked so far (see POSITIVE_CONTROLS.md).
# TE / category / 5'UTR length are read from source files, never typed; only `notes` is prose.
# Usage: Rscript code/finding_positive_controls/03_controls_summary.R
# Writes code/finding_positive_controls/controls_summary.csv
#
# TE is si3d knockdown, hypoxia 1hr, positive = eIF3d promotes (sictrl / si3d) in both lines.
#  - MDA-MB-231: output/translation_categories_si3d_vs_sictrl_hypoxia_1hr.csv
#  - MCF7-SIX1 TE: output/translation_categories_si3d_vs_sictrl_mcf7six1_hypoxia_1hr.csv
#    (the saved, shrunken LFC quoted in POSITIVE_CONTROLS.md)
#  - MCF7-SIX1 category: 00_mcf7six1_hypoxia_si3d_categories.Rmd output (same padj, asserted)
# 5'UTR lengths: output/finding_positive_controls/transcript_length_gc.csv (transcript_length_gc.R)
suppressPackageStartupMessages({ library(tidyverse); library(here) })

controls <- tribble(
  ~symbol,   ~control, ~tx_mda231,        ~tx_mcf7six1,      ~tx_proposed,
  "BCL2L12", "positive", "ENST00000246784", "ENST00000246784", "ENST00000441864",
  "TPRA1",   "positive", "ENST00000965987", "ENST00000355552", "ENST00000296210",
  # MCF7-SIX1's v26 transcript (ENST00000331789) adds 108 nt of promoter; the real leader in
  # both lines is the 84-nt MANE 5'UTR (POSITIVE_CONTROLS.md, ACTB)
  "ACTB",    "negative", "ENST00000646664", "ENST00000646664", NA
)

notes <- c(
  BCL2L12 = paste(
    "Strongest candidate; Exclusive in both lines. Proposed 510-nt 5'UTR is <=3% of BCL2L12",
    "mRNA in our cells and carries 5 upstream ATGs (2 in frame, would N-extend the protein),",
    "all in the extra 367 nt; the 143-nt 5'UTR our cells express (MANE) has none."),
  TPRA1 = paste(
    "MCF7-SIX1 only (231 just misses, padj 0.097). Buffered: mRNA rises on knockdown, so",
    "protein effect likely smaller than TE suggests. Low abundance (~3 TPM). Proposed 5'UTR",
    "matches MCF7-SIX1's leader (1 in-frame uATG, 12-codon uORF) but not 231's, which uses",
    "a different first exon."),
  ACTB = paste(
    "Negative control (Roiuk 2024). TE flat in every si3d contrast EXCEPT 231 hypoxia 24hr",
    "(TE LFC -1.14, padj 0.009). Abundant (~600-800 TPM). 231 1hr is Forwarded: mRNA and",
    "footprints fall together, TE does not change.")
)

read_te <- function(f) read_csv(here("output", f), show_col_types = FALSE) %>%
  mutate(gene_id_clean = sub("\\..*", "", gene_id))
te231 <- read_te("translation_categories_si3d_vs_sictrl_hypoxia_1hr.csv")
temcf <- read_te("translation_categories_si3d_vs_sictrl_mcf7six1_hypoxia_1hr.csv")
catmcf <- read_te(file.path("finding_positive_controls",
                            "mcf7six1_hypoxia_1hr_si3d_rna_ribo_te_categories.csv"))

# direction checks (CLAUDE.md): EIF3D's own RNA falls on si3d in 231; MCF7-SIX1 file is negated
stopifnot(te231$rna_lfc[te231$symbol == "EIF3D"] > 1.5,
          isTRUE(all.equal(temcf$te_lfc, -temcf$source_log2FoldChange)),
          temcf$te_lfc[temcf$symbol == "EIF3D"] < -1)

mcf <- temcf %>% dplyr::select(gene_id_clean, te_lfc, te_padj) %>%
  inner_join(catmcf %>% dplyr::select(gene_id_clean, category, padj_chk = te_padj),
             by = "gene_id_clean")
stopifnot("MCF7-SIX1 category file does not match the TE file" =
            max(abs(mcf$te_padj - mcf$padj_chk)) < 1e-8)

pick <- function(d, prefix) controls %>%
  dplyr::select(symbol) %>%
  inner_join(d %>% dplyr::select(symbol, te_lfc, te_padj, category) %>% distinct(), by = "symbol") %>%
  rename_with(~ paste0(prefix, .x), -symbol)
mcf <- mcf %>% left_join(temcf %>% dplyr::select(gene_id_clean, symbol), by = "gene_id_clean")
te_tab <- pick(te231, "mda231_") %>% inner_join(pick(mcf, "mcf7six1_"), by = "symbol")
stopifnot("a control gene is missing from a TE file, or duplicated" =
            nrow(te_tab) == nrow(controls))

len <- read_csv(here("output", "finding_positive_controls", "transcript_length_gc.csv"),
                show_col_types = FALSE) %>%
  filter(region == "utr5") %>% dplyr::select(transcript_id_clean, length_nt)
need <- na.omit(unique(c(controls$tx_mda231, controls$tx_mcf7six1, controls$tx_proposed)))
stopifnot("run transcript_length_gc.R on these transcripts first" =
            all(need %in% len$transcript_id_clean))
utr5_len <- function(tx) len$length_nt[match(tx, len$transcript_id_clean)]

out <- controls %>%
  inner_join(te_tab, by = "symbol") %>%
  mutate(across(ends_with("te_lfc"), ~ round(.x, 2)),
         across(ends_with("te_padj"), ~ signif(.x, 2)),
         utr5_nt_mda231 = utr5_len(tx_mda231),
         utr5_nt_mcf7six1 = utr5_len(tx_mcf7six1),
         utr5_nt_proposed = utr5_len(tx_proposed),
         notes = unname(notes[symbol])) %>%
  dplyr::select(gene = symbol, control,
                mda231_te_lfc, mda231_te_padj, mda231_category,
                mcf7six1_te_lfc, mcf7six1_te_padj, mcf7six1_category,
                utr5_tx_mda231 = tx_mda231, utr5_nt_mda231,
                utr5_tx_mcf7six1 = tx_mcf7six1, utr5_nt_mcf7six1,
                utr5_tx_proposed = tx_proposed, utr5_nt_proposed, notes)
stopifnot(!anyNA(out$notes))

write_csv(out, here("code", "finding_positive_controls", "controls_summary.csv"))
print(as.data.frame(out %>% dplyr::select(-notes)))
