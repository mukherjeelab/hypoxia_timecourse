# 18b_feature_table_final_regression.R
# Feature lists for the most recent si3d hypoxia 1hr REGRESSION models, split into
#   - Intrinsic features: computed from each mRNA's own sequence/structure
#   - External features:  measured data from published experiments
# with a membership column per model, so nothing is claimed to be in a model it was not in.
#
# Models (membership read from each saved rfreg_model_data_*.rds $feature_cols, never typed):
#   full        - intrinsic + external, aggregate GC          (only full-model regression run)
#   intrinsic   - intrinsic only, aggregate GC
#   posgc       - intrinsic only, GC1/GC2/GC3 with cds_gc + transcript_gc EXCLUDED by the
#                 spikein preset (10_regression_presets.R, exclude_features). The spike-in and
#                 no-spike arms are asserted to use the identical feature set.
#
# Descriptions are reused from 18_feature_table_for_slides.R (each row carries its named source);
# GC1/GC2 definitions come from feature_annotation_tco_posgc.csv (written by 01r_).
# Writes a NEW file and leaves feature_table_for_slides.xlsx untouched.
#   Rscript code/translational_control_overview/18b_feature_table_final_regression.R
suppressPackageStartupMessages({ library(dplyr); library(readr); library(here); library(openxlsx) })
o <- here("output", "translational_control_overview")
stem <- "rfreg_model_data_eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg_"
tail <- "corrected_both_clipaggregate_yraw.rds"
feats <- function(label) {
  p <- file.path(o, paste0(stem, label, tail))
  stopifnot("model data RDS not found" = file.exists(p))
  readRDS(p)$feature_cols
}
full      <- feats("")
intr      <- feats("intrinsic_")
posgc     <- feats("posgc_intrinsic_")
posgc_nsp <- feats("posgc_intrinsic_nospike_")
stopifnot(
  "spike-in and no-spike positional-GC models use different features" = setequal(posgc, posgc_nsp),
  "intrinsic model is not a subset of the full model" = all(intr %in% full),
  "positional-GC model should drop exactly cds_gc + transcript_gc and add gc1/gc2" =
    setequal(setdiff(intr, posgc), c("cds_gc", "transcript_gc")) &&
    setequal(setdiff(posgc, intr), c("gc1_internal", "gc2_internal"))
)

# Reuse the sourced description table from 18_ by evaluating only its tribble (no file writes).
src18 <- readLines(here("code", "translational_control_overview", "18_feature_table_for_slides.R"))
a <- grep("^acc5 <- ", src18); b <- grep("^stopifnot\\(\"a modelled feature has no description\"", src18)
stopifnot("could not locate the description block in 18_" = length(a) == 1, length(b) == 1, b > a)
eval(parse(text = src18[a:(b - 1)]))   # defines acc5/accC/acc3, SA/SS/SC/ST and the tribble `d`

ann_pg <- read_csv(file.path(o, "feature_annotation_tco_posgc.csv"), show_col_types = FALSE)
SP <- "feature_annotation_tco_posgc.csv (01r_positional_gc_matrix.R)"
pg_def <- function(f) ann_pg$definition[ann_pg$feature == f]
d <- bind_rows(d, tibble(
  feature = c("gc1_internal", "gc2_internal"),
  description = c("G+C at the first codon position; internal codons only (initiator and stop removed).",
                  "G+C at the second codon position; internal codons only (initiator and stop removed)."),
  units = "% (0-100)", source = SP))
# The full model uses the TRANSCRIPT-level CNOT3 score (_tx); 18_'s row describes the gene-level
# variant. Take the _tx definition from feature_annotation.csv, verbatim.
ann0 <- read_csv(file.path(o, "feature_annotation.csv"), show_col_types = FALSE)
cn <- ann0 %>% filter(feature == "cnot3_weighted_codon_score_tx")
stopifnot("no annotation row for cnot3_weighted_codon_score_tx" = nrow(cn) == 1, !is.na(cn$definition))
d <- bind_rows(d, tibble(feature = cn$feature, description = cn$definition, units = "score",
                         source = "feature_annotation.csv (01_build_feature_matrix.Rmd)"))
stopifnot("GC1/GC2 annotation definitions changed - re-check the wording above" =
            grepl("codon position 1", pg_def("gc1_internal")),
          grepl("codon position 2", pg_def("gc2_internal")))

all_feats <- union(full, posgc)
stopifnot("a modelled feature has no description" = all(all_feats %in% d$feature))

ann <- bind_rows(read_csv(file.path(o, "feature_annotation.csv"), show_col_types = FALSE) %>%
                   dplyr::select(feature, block, family),
                 ann_pg %>% filter(feature %in% c("gc1_internal", "gc2_internal")) %>%
                   dplyr::select(feature, block, family)) %>% distinct(feature, .keep_all = TRUE)

fam_label <- c(length = "Length", gc = "GC content", codon_optimality = "Codon optimality",
  initiation = "Start-codon context", structure_utr5 = "Structure: 5'UTR", structure_cds = "Structure: CDS",
  structure_utr3 = "Structure: 3'UTR", g4 = "RNA G-quadruplexes (G4)", nascent_peptide = "Nascent peptide",
  polya_track = "Poly(A) tracks", termination = "Stop codon / termination",
  clip = "eIF3 binding (CLIP / cap binding)", dap5 = "DAP5 dependence",
  stability_external = "mRNA stability & other factors")
fam_order <- names(fam_label)
gc_order  <- c("utr5_gc", "cds_gc", "utr3_gc", "transcript_gc", "gc1_internal", "gc2_internal", "gc3_internal")

gc_note <- c(
  cds_gc        = "Not in the positional-GC model: removed so the forest cannot use the aggregate in place of GC1/GC2/GC3 (CDS GC is ~ the mean of the three positions).",
  transcript_gc = "Not in the positional-GC model: removed together with CDS GC (it contains CDS GC).",
  gc1_internal  = "Only in the positional-GC model (no full-model run with GC1/GC2 exists).",
  gc2_internal  = "Only in the positional-GC model (no full-model run with GC1/GC2 exists).")

tick <- function(x) ifelse(x, "✓", "")
tab <- d %>% filter(feature %in% all_feats) %>% distinct(feature, .keep_all = TRUE) %>%
  left_join(ann, by = "feature") %>%
  mutate(family_label = fam_label[family],
         ord = ifelse(feature %in% gc_order, match(feature, gc_order), 0)) %>%
  arrange(match(family, fam_order), ord, feature) %>%
  transmute(block,
            `Feature family` = family_label, Feature = feature, `What it measures` = description,
            Units = units,
            `Positional-GC intrinsic model (57)` = tick(feature %in% posgc),
            `Intrinsic model (57)` = tick(feature %in% intr),
            `Full model (69)` = tick(feature %in% full),
            Note = ifelse(feature %in% names(gc_note), gc_note[feature], ""),
            `Description source` = source)
stopifnot("family/block missing" = !any(is.na(tab$`Feature family`) | is.na(tab$block)),
          "unexpected block" = all(tab$block %in% c("intrinsic", "external")),
          "posgc model contains an external feature" = all(posgc %in% tab$Feature[tab$block == "intrinsic"]))

intr_tab <- tab %>% filter(block == "intrinsic") %>% dplyr::select(-block)
ext_tab  <- tab %>% filter(block == "external")  %>%
  dplyr::select(-block, -`Positional-GC intrinsic model (57)`, -`Intrinsic model (57)`)
stopifnot(nrow(intr_tab) == length(union(intr, posgc)), nrow(ext_tab) == length(full) - length(intr),
          sum(intr_tab$`Positional-GC intrinsic model (57)` == "✓") == length(posgc))

summ <- tab %>% mutate(family = names(fam_label)[match(`Feature family`, fam_label)]) %>%
  group_by(block, family) %>%
  summarise(posgc = sum(Feature %in% posgc), intrinsic = sum(Feature %in% intr),
            full = sum(Feature %in% full), .groups = "drop") %>%
  arrange(desc(block == "intrinsic"), match(family, fam_order)) %>%
  transmute(Type = ifelse(block == "intrinsic", "Intrinsic", "External"),
            `Feature family` = fam_label[family],
            `Positional-GC intrinsic model` = posgc, `Intrinsic model` = intrinsic, `Full model` = full)
summ <- bind_rows(summ, tibble(Type = "TOTAL", `Feature family` = "",
  `Positional-GC intrinsic model` = length(posgc), `Intrinsic model` = length(intr), `Full model` = length(full)))
stopifnot(sum(head(summ$`Full model`, -1)) == length(full),
          sum(head(summ$`Positional-GC intrinsic model`, -1)) == length(posgc))

notes <- tibble(Note = c(
  "All models: random-forest regression of si3d vs sictrl TE log2 fold change, MDA-MB-231 hypoxia 1hr, all testable genes.",
  sprintf("Positional-GC intrinsic model (most recent, 2026-11-20): %d intrinsic features, GC split into GC1/GC2/GC3. CDS GC and transcript GC were deliberately excluded. The spike-in and no-spike-in versions use the identical feature set (asserted).", length(posgc)),
  sprintf("Intrinsic model (2026-11-18): %d intrinsic features with aggregate GC (CDS GC, transcript GC, GC3) and no GC1/GC2.", length(intr)),
  sprintf("Full model (2026-10-01): %d features = the intrinsic model's %d plus %d external features. No full model has been run with GC1/GC2, so the external sheet comes from this run only.", length(full), length(intr), length(full) - length(intr)),
  "No model contains all five of GC1, GC2, GC3, CDS GC and transcript GC together. All five are listed on the intrinsic sheet; the model columns show which run used each.",
  "Membership is read from each saved model's feature list, not typed. Descriptions carry a named source (Description source column)."))

wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#E8E8E8", border = "Bottom", wrapText = TRUE, valign = "top")
ws <- createStyle(wrapText = TRUE, valign = "top")
gcs <- createStyle(fgFill = "#FFF4D6", wrapText = TRUE, valign = "top")
for (nm in c("Summary", "Intrinsic features", "External features", "Notes")) addWorksheet(wb, nm)
writeData(wb, "Summary", summ, headerStyle = hs)
setColWidths(wb, "Summary", 1:5, c(12, 34, 18, 14, 12))
addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = nrow(summ) + 1, cols = 1:5, gridExpand = TRUE)

writeData(wb, "Intrinsic features", intr_tab, headerStyle = hs)
setColWidths(wb, "Intrinsic features", 1:ncol(intr_tab), c(24, 40, 70, 12, 16, 12, 10, 50, 45))
addStyle(wb, "Intrinsic features", ws, rows = 2:(nrow(intr_tab) + 1), cols = 1:ncol(intr_tab), gridExpand = TRUE)
gc_rows <- which(intr_tab$`Feature family` == "GC content") + 1
addStyle(wb, "Intrinsic features", gcs, rows = gc_rows, cols = 1:ncol(intr_tab), gridExpand = TRUE)
freezePane(wb, "Intrinsic features", firstRow = TRUE)

writeData(wb, "External features", ext_tab, headerStyle = hs)
setColWidths(wb, "External features", 1:ncol(ext_tab), c(30, 32, 80, 14, 10, 10, 50))
addStyle(wb, "External features", ws, rows = 2:(nrow(ext_tab) + 1), cols = 1:ncol(ext_tab), gridExpand = TRUE)
freezePane(wb, "External features", firstRow = TRUE)

writeData(wb, "Notes", notes, headerStyle = hs); setColWidths(wb, "Notes", 1, 120)
addStyle(wb, "Notes", ws, rows = 2:(nrow(notes) + 1), cols = 1)

out <- file.path(o, "feature_table_final_regression.xlsx")
saveWorkbook(wb, out, overwrite = TRUE)
cat("wrote", out, "\n"); print(as.data.frame(summ), row.names = FALSE)
print(as.data.frame(intr_tab %>% filter(`Feature family` == "GC content") %>%
  dplyr::select(Feature, 5:7)), row.names = FALSE)
