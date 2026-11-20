# Is each top eIF3d-promoted gene bound by eIF3 in published datasets?
# Reads every source table directly (not the feature matrix, whose clip_lee_glu_dep_bound
# column holds the complete-media set - see CAP_BINDING_DATA.md).
#
# Usage: Rscript code/finding_positive_controls/05_eif3_binding_lookup.R [n_top | SYMBOLS]
#   n_top (default 7) = how many genes to take from 04_'s ranking, or comma-separated symbols
#   (e.g. BSG,ACTB) to look up any genes; outputs then get a _SYMBOLS suffix
#   (output/finding_positive_controls/04_top_eif3d_targets_hypoxia_both_lines.csv)
# Writes output/finding_positive_controls/05_eif3_binding_lookup_{long,summary}.csv
#
# Datasets (all in /Users/katematlin/github/mpra_preliminary/data/; see
# code/translational_control_overview/EIF3_CLIP_DATA.md and CAP_BINDING_DATA.md):
#   PAR-CLIP  Lee 2015 (HEK293T, eIF3a/b/d/g)         lee_parclip.csv               Ensembl gene
#   PAR-CLIP  De Silva 2021 (Jurkat, 5 bands x +/-act)  GSE191306_...clusters.xlsx    symbol only
#   irCLIP    Cate lab 2023 (NPCs), HOMER peaks         S2_cate_2023 / S3_cate_2023   transcript
#   subunit-seq Lamper 2020 (HEK293T, CM + GLU)         lamper-s3.xlsx (enriched only) Ensembl gene
#   subunit-seq Mukhopadhyay 2023 (HEK293T, Tg)         integrated_stress_response.xlsx Ensembl gene
# Absence from a thresholded table means "not called", not "measured and unbound": detection
# follows expression (EIF3_CLIP_DATA.md), and none of these are our cell lines.
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(readxl); library(GenomicFeatures)
  library(AnnotationDbi); library(org.Hs.eg.db)
})
args <- commandArgs(trailingOnly = TRUE)
by_symbol <- length(args) > 0 && grepl("[A-Za-z]", args[1])
n_top <- if (length(args) && !by_symbol) as.integer(args[1]) else 7L
D <- "/Users/katematlin/github/mpra_preliminary/data"
strip <- function(x) sub("\\..*", "", x)

ranked <- read_csv(here("output", "finding_positive_controls",
                        "04_top_eif3d_targets_hypoxia_both_lines.csv"), show_col_types = FALSE)
cand <- ranked %>% slice_head(n = n_top) %>%
  dplyr::select(rank, gene_id_clean, symbol, lfc_231, lfc_mcf, mean_lfc)
if (by_symbol) {
  syms <- trimws(strsplit(args[1], ",")[[1]])
  cand <- tibble(rank = seq_along(syms), symbol = syms,
                 gene_id_clean = unname(mapIds(org.Hs.eg.db, syms, "ENSEMBL", "SYMBOL"))) %>%
    left_join(ranked %>% dplyr::select(gene_id_clean, lfc_231, lfc_mcf, mean_lfc), by = "gene_id_clean")
  stopifnot("a symbol did not map to an Ensembl gene" = !anyNA(cand$gene_id_clean))
}
sfx <- if (by_symbol) paste0("_", gsub(",", "_", args[1])) else ""
# spot-check genes with a known answer (removed before output)
checks <- tibble(symbol = c("JUN", "ALKBH5")) %>%
  mutate(gene_id_clean = mapIds(org.Hs.eg.db, symbol, "ENSEMBL", "SYMBOL"))
genes <- bind_rows(cand, checks)

# symbols + aliases for the symbol-only table (old names, e.g. MARCHF3 was MARCH3)
alias <- AnnotationDbi::select(org.Hs.eg.db, keys = genes$gene_id_clean, keytype = "ENSEMBL",
                               columns = c("SYMBOL", "ALIAS")) %>%
  as_tibble() %>% pivot_longer(c(SYMBOL, ALIAS), names_to = "kind", values_to = "name") %>%
  distinct(gene_id_clean = ENSEMBL, name) %>%
  bind_rows(genes %>% transmute(gene_id_clean, name = symbol)) %>% distinct() %>%
  mutate(name_u = toupper(name))

# ---- Lee 2015 PAR-CLIP -------------------------------------------------------------------
lee <- read_csv(file.path(D, "lee_parclip.csv"), show_col_types = FALSE) %>%
  transmute(gene_id_clean = strip(`Ensembl ID`), subunit = Subunit, region = Annotation, reads = Reads)
lee_hits <- lee %>% filter(gene_id_clean %in% genes$gene_id_clean) %>%
  group_by(gene_id_clean, subunit) %>%
  summarise(evidence = sprintf("%d cluster(s), %d reads [%s]", n(), sum(reads),
                               paste(sort(unique(region)), collapse = "; ")), .groups = "drop") %>%
  transmute(gene_id_clean, dataset = "PAR-CLIP Lee 2015 (HEK293T)",
            condition = paste0("eIF3", subunit), evidence)

# ---- De Silva 2021 Jurkat PAR-CLIP --------------------------------------------------------
ds_file <- file.path(D, "GSE191306_DeSilva-Jurkat-eIF3-PAR-CLIP-clusters.xlsx")
ds_key <- read_excel(ds_file, sheet = "Samples") %>% filter(!is.na(`Name of Tab`)) %>%
  transmute(sheet = `Name of Tab`, sample = `PAR-CLIP Samples`)
region_cols <- c("5'utr", "Intron", "Exon", "3'utr", "Start_codon", "Stop_codon", "Junction")
ds <- map_dfr(ds_key$sheet, ~ read_excel(ds_file, sheet = .x) %>% mutate(sheet = .x)) %>%
  left_join(ds_key, by = "sheet")
stopifnot(all(region_cols %in% names(ds)))
ds_hits <- ds %>% mutate(name_u = toupper(GeneName)) %>%
  inner_join(alias %>% dplyr::select(gene_id_clean, name_u), by = "name_u",
             relationship = "many-to-many") %>%
  rowwise() %>%
  mutate(regions = paste(region_cols[c_across(all_of(region_cols)) > 0], collapse = "; ")) %>%
  ungroup() %>%
  transmute(gene_id_clean, dataset = "PAR-CLIP De Silva 2021 (Jurkat)",
            condition = sample,
            evidence = sprintf("%s: %d site(s), %d reads [%s]", GeneName, Sites, ReadCountSum,
                               regions))

# ---- Cate lab NPC irCLIP peaks (transcript-keyed) -----------------------------------------
txdb <- loadDb(here("accessories", "human", "txdb.gencode49.sqlite"))
tx2g <- AnnotationDbi::select(txdb, keys = keys(txdb, "GENEID"), keytype = "GENEID",
                              columns = "TXNAME") %>%
  transmute(tx = strip(TXNAME), gene_id_clean = strip(GENEID)) %>%
  filter(gene_id_clean %in% genes$gene_id_clean)
cate <- map_dfr(c(S2 = "S2_cate_2023.xlsx", S3 = "S3_cate_2023.xlsx"), .id = "table", function(f)
  read_excel(file.path(D, f), col_types = "text") %>%
    transmute(PeakID, hit_tx = strip(Filter2_Hit), hit_feature = Filter2_Hit_feature_types,
              all_tx = All_ENSTs))
cate_hits <- cate %>% inner_join(tx2g, by = c("hit_tx" = "tx")) %>%
  group_by(table, gene_id_clean) %>%
  summarise(evidence = sprintf("%d peak(s) [%s]", n(), paste(sort(unique(hit_feature)), collapse = "; ")),
            .groups = "drop") %>%
  transmute(gene_id_clean, dataset = "irCLIP Cate lab 2023 (NPCs)",
            condition = paste("table", table), evidence)
# peaks touching the gene but assigned to another transcript (reported, not called a hit)
cate_any <- cate %>% mutate(tx = str_extract_all(all_tx, "ENST\\d+")) %>% unnest(tx) %>%
  inner_join(tx2g, by = "tx") %>% distinct(table, PeakID, gene_id_clean) %>%
  dplyr::count(table, gene_id_clean, name = "peaks_overlapping")

# ---- Lamper 2020 subunit-seq (enriched genes only) ----------------------------------------
lam <- read_excel(file.path(D, "lamper-s3.xlsx")) %>% mutate(gene_id_clean = strip(ensembl_gene_id))
lam_hits <- lam %>% filter(gene_id_clean %in% genes$gene_id_clean) %>%
  transmute(gene_id_clean,
            CM  = sprintf("IP/IN %.1f, log2FC %.2f, padj %.2g", `IP/IN_CM (AVG)`,
                          as.numeric(log2FoldChange_CM), as.numeric(padj_CM)),
            GLU = sprintf("IP/IN %.1f, log2FC %.2f, padj %.2g", `IP/IN_GLU (AVG)`,
                          as.numeric(log2FoldChange_Glu), as.numeric(padj_Glu))) %>%
  pivot_longer(c(CM, GLU), names_to = "condition", values_to = "evidence") %>%
  mutate(dataset = "subunit-seq Lamper 2020 (HEK293T)",
         condition = recode(condition, CM = "complete media", GLU = "glucose deprivation"))

# ---- Mukhopadhyay 2023 subunit-seq, thapsigargin -------------------------------------------
isr <- read_excel(file.path(D, "integrated_stress_response.xlsx")) %>%
  mutate(gene_id_clean = strip(ensembl_gene_id))
isr_hits <- isr %>% filter(gene_id_clean %in% genes$gene_id_clean) %>%
  transmute(gene_id_clean, dataset = "subunit-seq Mukhopadhyay 2023 (HEK293T)",
            condition = "thapsigargin",
            evidence = sprintf("log2FC %.2f, padj %.2g", as.numeric(log2FoldChange), as.numeric(padj)))

long <- bind_rows(lee_hits, ds_hits, cate_hits, lam_hits, isr_hits) %>%
  left_join(genes %>% dplyr::select(gene_id_clean, symbol), by = "gene_id_clean")

# ---- checks ---------------------------------------------------------------------------------
stopifnot("Lee 2015 should report JUN (named eIF3 target)" =
            any(long$symbol == "JUN" & grepl("Lee", long$dataset)),
          "Mukhopadhyay 2023 should report ALKBH5 (named target)" =
            any(long$symbol == "ALKBH5" & grepl("Mukhopadhyay", long$dataset)),
          "RPP14 thapsigargin log2FC should be 1.51 (Table S2)" =
            abs(as.numeric(isr$log2FoldChange[isr$gene_id_clean == "ENSG00000163684"]) - 1.514) < 1e-3)
# agree with the 231 feature matrix's bound flags where those are known to be right
fm <- readRDS(here("output", "translational_control_overview", "feature_matrix_tco.rds")) %>%
  filter(gene_id_clean %in% cand$gene_id_clean) %>%
  transmute(gene_id_clean, lee_fm = clip_lee_parclip_a_b_d_g_bound > 0,
            isr_fm = clip_lee_int_stress_bound > 0,
            cate_fm = clip_cate_undiff_npcs_bound > 0 | clip_cate_diff_npcs_bound > 0)
agree <- fm %>% mutate(lee_here  = gene_id_clean %in% lee_hits$gene_id_clean,
                       isr_here  = gene_id_clean %in% isr_hits$gene_id_clean,
                       cate_here = gene_id_clean %in% cate_hits$gene_id_clean)
print(as.data.frame(agree))
stopifnot("Lee / thapsigargin calls disagree with the feature matrix" =
            all(agree$lee_fm == agree$lee_here), all(agree$isr_fm == agree$isr_here))
if (!all(agree$cate_fm == agree$cate_here))
  message("NOTE: Cate calls differ from the feature matrix (its mapping is undocumented)")

# ---- output ---------------------------------------------------------------------------------
long <- long %>% filter(gene_id_clean %in% cand$gene_id_clean)
summary <- cand %>%
  left_join(cate_any %>% group_by(gene_id_clean) %>%
              summarise(cate_peaks_overlapping = sum(peaks_overlapping)), by = "gene_id_clean") %>%
  mutate(
    lee_parclip      = map_chr(gene_id_clean, ~ paste(long$condition[long$gene_id_clean == .x & grepl("Lee", long$dataset)], collapse = ", ")),
    desilva_eif3d_l  = map_chr(gene_id_clean, ~ paste(long$condition[long$gene_id_clean == .x & grepl("De Silva", long$dataset) & grepl("EIF3D", long$condition)], collapse = " | ")),
    desilva_other    = map_chr(gene_id_clean, ~ paste(long$condition[long$gene_id_clean == .x & grepl("De Silva", long$dataset) & !grepl("EIF3D", long$condition)], collapse = " | ")),
    cate_irclip      = map_chr(gene_id_clean, ~ paste(long$condition[long$gene_id_clean == .x & grepl("Cate", long$dataset)], collapse = ", ")),
    cap_lamper2020   = map_chr(gene_id_clean, ~ paste(long$condition[long$gene_id_clean == .x & grepl("Lamper", long$dataset)], collapse = ", ")),
    cap_thapsigargin = map_chr(gene_id_clean, ~ paste(long$evidence[long$gene_id_clean == .x & grepl("Mukhopadhyay", long$dataset)], collapse = ", ")),
    across(c(lfc_231, lfc_mcf, mean_lfc), ~ round(.x, 2))) %>%
  mutate(across(lee_parclip:cap_thapsigargin, ~ na_if(.x, "")))

out <- here("output", "finding_positive_controls")
write_csv(long %>% arrange(match(gene_id_clean, cand$gene_id_clean), dataset, condition),
          file.path(out, paste0("05_eif3_binding_lookup_long", sfx, ".csv")))
write_csv(summary, file.path(out, paste0("05_eif3_binding_lookup_summary", sfx, ".csv")))
print(as.data.frame(long %>% arrange(match(gene_id_clean, cand$gene_id_clean)) %>%
                      dplyr::select(symbol, dataset, condition, evidence)))
print(as.data.frame(summary %>% dplyr::select(rank, symbol, cate_peaks_overlapping)))
