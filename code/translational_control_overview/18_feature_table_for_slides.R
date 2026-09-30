# 18_feature_table_for_slides.R
# Slide-ready table of the features the heatmap models use: 69 in the full model, 57 of them in
# the intrinsic model. Feature membership is read from the saved models (not typed); the
# descriptions are taken ONLY from a named source per row (feature_annotation.csv, the code that
# builds the feature, or the collaborator repo docs), recorded in the "Description source" column.
# Nothing is added that the source does not state.
#   Rscript code/translational_control_overview/18_feature_table_for_slides.R
suppressPackageStartupMessages({ library(dplyr); library(readr); library(here); library(openxlsx) })
o <- here("output", "translational_control_overview")
full <- read_csv(file.path(o, "rfreg_importance_eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg_corrected_both_clipaggregate_yraw.csv"), show_col_types = FALSE)$feature
intr <- read_csv(file.path(o, "rfreg_importance_eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg_intrinsic_corrected_both_clipaggregate_yraw.csv"), show_col_types = FALSE)$feature
ann  <- read_csv(file.path(o, "feature_annotation.csv"), show_col_types = FALSE)
stopifnot("intrinsic set is not a subset of the full set" = all(intr %in% full))

acc5 <- "RNAplfold -W 80 -L 40 -u 30, 5'UTR folded on its own; per-nucleotide probability of being unpaired."
accC <- "RNAplfold -W 150 -L 100 -u 30, CDS folded on its own; per-nucleotide probability of being unpaired."
acc3 <- "RNAplfold -W 150 -L 100 -u 30, 3'UTR folded on its own; per-nucleotide probability of being unpaired."
SA <- "feature_annotation.csv (01_build_feature_matrix.Rmd)"; SS <- "feature_annotation.csv; 30-nt windows from 01i_build_mcf7six1_baseline.R feats()"; SC <- "collaborator repo docs/translation_challenge_features.md (mukherjeelab/2026_TranslationResponseFeatureAnalysis)"; ST <- "collaborator repo docs/termination_features.md (mukherjeelab/2026_TranslationResponseFeatureAnalysis)"
d <- tribble(~feature, ~description, ~units, ~source,
 "log2_utr5_length",  "Length of the 5'UTR.", "log2(nt + 1)", SA,
 "log2_cds_length",   "Length of the coding sequence.", "log2(nt + 1)", SA,
 "log2_utr3_length",  "Length of the 3'UTR.", "log2(nt + 1)", SA,
 "log2_total_length", "Total length: 5'UTR + CDS + 3'UTR.", "log2(nt + 1)", "code: predictive_modeling/01_feature_extraction.Rmd (total_length, log2_total_length)",
 "utr5_gc",       "G+C content of the 5'UTR (denominator includes ambiguous bases).", "% (0-100)", SA,
 "cds_gc",        "G+C content of the CDS (denominator includes ambiguous bases).", "% (0-100)", SA,
 "utr3_gc",       "G+C content of the 3'UTR (denominator includes ambiguous bases).", "% (0-100)", SA,
 "transcript_gc", "G+C content of the transcript: length-weighted mean of the three regional GCs.", "% (0-100)", SA,
 "gc3_internal",  "G+C at the third codon position; internal codons only (initiator and stop removed).", "% (0-100)", SA,
 "cai_internal", "Codon Adaptation Index (Sharp & Li), ribosomal-protein (RPL/RPS) reference; internal codons only.", "0-1", SA,
 "fop_internal", "Frequency of optimal codons, optimal codons defined from the top-10%-TE reference; internal codons only.", "0-1", SA,
 "tai_gtrnadb",  "tRNA Adaptation Index (dos Reis-style) from GtRNAdb hg38 tRNA gene copy numbers with the collaborator wobble scheme v1.", "0-1", SA,
 "csc_internal", "iCodon predict_stability('human') on the internal-codon CDS.", "score", SA,
 "kozak_pwm_score_v2",     "Kozak 1987 Table 1 log-odds score of the start-codon context (imported file); NA unless the start is AUG.", "log-odds", SA,
 "noderer_tis_efficiency", "Noderer 2014 FACS-seq translation-initiation-site efficiency for the -6 to +5 context; NA unless the start is AUG.", "12-150", SA,
 "kozak_score",   "Kozak tier: 3 optimal (GCCACC or GCCGCC before ATG, G at +4); 2 strong (A/G at -3 and G at +4); 1 moderate (A/G at -3 without G at +4, or C/T at -3 with G at +4); 0 otherwise.", "0-3", "code: predictive_modeling/01h_feature_extraction_kozak.Rmd (match_kozak)",
 "kozak_optimal", "1 if kozak_score is 3 (optimal context), else 0.", "0/1", "code: predictive_modeling/01h_feature_extraction_kozak.Rmd",
 "tco_struct_accessibility_utr5_mean",    paste("Mean over the 5'UTR.", acc5), "0-1", SS,
 "tco_struct_accessibility_utr5_min",     paste("Minimum over the 5'UTR.", acc5), "0-1", SS,
 "tco_struct_accessibility_utr5_sd",      paste("Standard deviation over the 5'UTR.", acc5), "SD", SS,
 "tco_struct_accessibility_cap_proximal", paste("Mean over the first 30 nt of the 5'UTR.", acc5), "0-1", SS,
 "tco_struct_accessibility_aug_context",  paste("Mean over the last 30 nt of the 5'UTR.", acc5), "0-1", SS,
 "tco_struct_accessibility_cds_mean",           paste("Mean over the CDS.", accC), "0-1", SS,
 "tco_struct_accessibility_cds_min",            paste("Minimum over the CDS.", accC), "0-1", SS,
 "tco_struct_accessibility_cds_sd",             paste("Standard deviation over the CDS.", accC), "SD", SS,
 "tco_struct_accessibility_start_proximal_cds", paste("Mean over the first 30 nt of the CDS.", accC), "0-1", SS,
 "tco_struct_accessibility_stop_proximal_cds",  paste("Mean over the last 30 nt of the CDS.", accC), "0-1", SS,
 "tco_struct_accessibility_utr3_mean",          paste("Mean over the 3'UTR.", acc3), "0-1", SS,
 "tco_struct_accessibility_utr3_min",           paste("Minimum over the 3'UTR.", acc3), "0-1", SS,
 "tco_struct_accessibility_utr3_sd",            paste("Standard deviation over the 3'UTR.", acc3), "SD", SS,
 "tco_struct_accessibility_stop_proximal_utr3", paste("Mean over the first 30 nt of the 3'UTR.", acc3), "0-1", SS,
 "tco_struct_accessibility_distal_utr3",        paste("Mean over the last 30 nt of the 3'UTR.", acc3), "0-1", SS,
 "g4mer_mean_utr5",      "Mean G4mer (RNA G-quadruplex model) window score, 5'UTR; no length bias.", "score", SA,
 "g4mer_mean_cds",       "Mean G4mer window score, CDS; no length bias.", "score", SA,
 "g4mer_mean_utr3",      "Mean G4mer window score, 3'UTR; no length bias.", "score", SA,
 "g4mer_frac_above_utr5","Fraction of 5'UTR windows above the rG4 threshold.", "0-1", SA,
 "g4mer_frac_above_cds", "Fraction of CDS windows above the rG4 threshold.", "0-1", SA,
 "g4mer_frac_above_utr3","Fraction of 3'UTR windows above the rG4 threshold.", "0-1", SA,
 "g4mer_max_resid_utr5", "Maximum G4mer window score in the 5'UTR, as a residual on log2 5'UTR length + 5'UTR GC (stride 10).", "residual", SA,
 "g4mer_max_resid_cds",  "Maximum G4mer window score in the CDS, as a residual on log2 CDS length + CDS GC (stride 20).", "residual", SA,
 "g4mer_max_resid_utr3", "Maximum G4mer window score in the 3'UTR, as a residual on log2 3'UTR length + 3'UTR GC (stride 20).", "residual", SA,
 "proline_fraction",    "Fraction (0-1) of internal-codon residues that are proline. Collaborator note: nascent-peptide features are sequence-derived and do not directly measure ribosome pausing or elongation rates.", "0-1", paste(SA, "+", SC),
 "ppp_motif_density",   "PPP occurrences per 100 residues; overlapping matches counted.", "per 100 aa", SA,
 "max_net_charge_30aa", "Maximum (K+R) - (D+E) over any 30-amino-acid window; whole chain if shorter than 30.", "charge", SA,
 "max_consecutive_aaa", "Longest in-frame run of AAA lysine codons. Collaborator note: AAA and AAG encode the same lysine peptide but create different nucleotide tracks.", "codons", paste(SA, "+", SC),
 "max_consecutive_aag", "Longest in-frame run of AAG lysine codons.", "codons", SA,
 "term_stop_TAA", "1 if the annotated stop codon is TAA, else 0.", "0/1", SA,
 "term_stop_TAG", "1 if the annotated stop codon is TAG, else 0.", "0/1", SA,
 "term_stop_TGA", "1 if the annotated stop codon is TGA, else 0.", "0/1", SA,
 "term_plus4_A",  "1 if the base immediately 3' of the stop (+4) is A, else 0.", "0/1", SA,
 "term_plus4_C",  "1 if the +4 base is C, else 0.", "0/1", SA,
 "term_plus4_G",  "1 if the +4 base is G, else 0.", "0/1", SA,
 "term_plus4_T",  "1 if the +4 base is T, else 0.", "0/1", SA,
 "term_window_gc","Percent GC of the 9-nt window: stop codon plus +4 to +9.", "% (0-100)", SA,
 "term_readthrough_cridge", "Cridge 2018 Table S1 percent readthrough for this stop/+4/+5/+6 context (192-context lookup). Collaborator note: higher = higher reporter readthrough and lower termination fidelity in that reporter system; not a direct measurement of endogenous termination efficiency.", "% readthrough", paste(SA, "+", ST),
 "term_log1p_readthrough",  "log1p(term_readthrough_cridge).", "log(1 + %)", SA,
 "term_interval_codons",    "Codons from the annotated stop to the next in-frame stop in the 3'UTR.", "codons", SA,
 "clip_eif3_count",            "Number of eIF3 CLIP datasets calling the gene bound: Lee 2015 PAR-CLIP (subunits a/b/d/g), Cate-lab NPC undifferentiated, Cate-lab NPC differentiated.", "0-3", "code: predictive_modeling/01e_feature_extraction_clip_split.Rmd (eif3_conditions)",
 "clip_cap_binding_count_v2",  "Cap-binding CLIP (subunit-seq) datasets with evidence, counting the 2 distinct conditions. Caveat found 2026-09-30: the glucose-deprivation flag holds the complete-media gene set (CAP_BINDING_DATA.md).", "0-2", paste(SA, "+ CAP_BINDING_DATA.md"),
 "dap5_polysome_lfc",          "DAP5 polysome DESeq2 interaction log2FC (GSE115142), non-silencing control vs DAP5 shRNA, heavy fraction; positive = higher in control = promoted by DAP5.", "log2 FC", "code: predictive_modeling/01b_feature_extraction_targeted.Rmd",
 "is_dap5_promoted",           "1 if the gene is labelled 'down' (reduced TE in DAP5 KO) in Weber 2022 Table S1.", "0/1", "code: predictive_modeling/01_feature_extraction.Rmd (dap5_promoted_symbols)",
 "is_dap5_repressed",          "1 if the gene is labelled 'up' (increased TE in DAP5 KO) in Weber 2022 Table S1.", "0/1", "code: predictive_modeling/01_feature_extraction.Rmd (dap5_repressed_symbols)",
 "cnot3_ribo_enrichment",      "Zhu 2024 CNOT3 selective ribosome profiling enrichment, as -log10(minimum FDR); higher = more enriched.", "-log10 FDR", "code: predictive_modeling/22_cnot3_riboseq_slamseq.Rmd",
 "cnot3_weighted_codon_score", "Sum over sense codons of (codon count x that codon's Spearman r from the Zhu 2024 CNOT3 codon stability scores), divided by total sense codons.", "score", "code: predictive_modeling/43_feature_extraction_mcf7six1_codon.Rmd (cnot3_ws); carried from the 231 baseline",
 "zhu2024_ctrl_halflife",      "Zhu 2024 SLAM-seq mRNA half-life in control (sgNT) cells.", "as published", "code: predictive_modeling/01i_feature_extraction_external_stability.Rmd (sgNT_halflife)",
 "zhu2024_cnot3ko_halflife_logfc", "Zhu 2024 log2 ratio of half-life, CNOT3 knockout vs control.", "log2 ratio", "code: predictive_modeling/01i_feature_extraction_external_stability.Rmd (log2_ratio)",
 "karner2026_mda231_log2_ct",  "Karner 2026 MDA-MB-231 SLAM-seq C-to-T conversion ratio, log2 of the mean of the two replicates.", "log2 ratio", "code: predictive_modeling/01i_ + 25_goodarzi_231_decay_rates.Rmd",
 "hia2026_dhx29_occupancy_tx", "Hia 2026 DHX29 ribosome occupancy, joined on transcript_id_clean (the source key).", "score", SA,
 "hia2026_dhx29_kd_rna_lfc",   "Hia 2026 Table S2 mRNA logFC after DHX29 knockdown.", "log2 FC", "code: predictive_modeling/87_dhx29_rnaseq_gc3_sensitivity.Rmd")

stopifnot("a modelled feature has no description" = setequal(d$feature, full))

fam_label <- c(length = "Length", gc = "GC content", codon_optimality = "Codon optimality",
  initiation = "Start-codon context", structure_utr5 = "Structure: 5'UTR", structure_cds = "Structure: CDS",
  structure_utr3 = "Structure: 3'UTR", g4 = "RNA G-quadruplexes (G4)", nascent_peptide = "Nascent peptide",
  polya_track = "Poly(A) tracks", termination = "Stop codon / termination",
  clip = "eIF3 binding (CLIP / cap binding)", dap5 = "DAP5 dependence", stability_external = "mRNA stability & other factors")
fam_what <- c(length = "How long each region is.", gc = "Nucleotide composition by region and at the codon wobble position.",
  codon_optimality = "Codon-usage indices: CAI, FOP, tAI and iCodon CSC.", initiation = "Scores of the sequence around the start codon (Kozak PWM, Kozak tier, Noderer efficiency).",
  structure_utr5 = "RNAplfold accessibility of the 5'UTR.", structure_cds = "RNAplfold accessibility of the CDS.",
  structure_utr3 = "RNAplfold accessibility of the 3'UTR.", g4 = "G4mer RNA G-quadruplex scores, by region.",
  nascent_peptide = "Sequence-derived properties of the encoded peptide (proline content, PPP motifs, local charge).", polya_track = "Longest runs of AAA and AAG lysine codons in the CDS.",
  termination = "Stop codon identity and its surrounding sequence.", clip = "Counts of published eIF3 CLIP / cap-binding datasets calling the gene bound.",
  dap5 = "DAP5 dependence from published polysome and ribosome-profiling data.", stability_external = "Published mRNA stability (SLAM-seq), CNOT3 and DHX29 data.")
fam_order <- names(fam_label)

tab <- d %>% left_join(ann %>% dplyr::select(feature, family), by = "feature") %>%
  mutate(family_label = fam_label[family],
         full_model = "✓", intrinsic = ifelse(feature %in% intr, "✓", ""),
         type = ifelse(feature %in% intr, "Intrinsic (sequence/structure)", "External (measured data)")) %>%
  arrange(match(family, fam_order), feature) %>%
  transmute(`Feature family` = family_label, Feature = feature, `What it measures` = description,
            Units = units, `Full model (69)` = full_model, `Intrinsic model (57)` = intrinsic, Type = type,
            `Description source` = source)
stopifnot("family labels missing" = !any(is.na(tab$`Feature family`)))

summ <- d %>% left_join(ann %>% dplyr::select(feature, family), by = "feature") %>%
  group_by(family) %>% summarise(full = n(), intrinsic = sum(feature %in% intr), .groups = "drop") %>%
  arrange(match(family, fam_order)) %>%
  transmute(`Feature family` = fam_label[family], `What the family captures` = fam_what[family],
            `Full model` = full, `Intrinsic model` = intrinsic)
summ <- bind_rows(summ, tibble(`Feature family` = "TOTAL", `What the family captures` = "",
                               `Full model` = sum(summ$`Full model`), `Intrinsic model` = sum(summ$`Intrinsic model`)))
stopifnot(sum(summ$`Full model`[summ$`Feature family` != "TOTAL"]) == length(full),
          sum(summ$`Intrinsic model`[summ$`Feature family` != "TOTAL"]) == length(intr))

notes <- tibble(Note = c(
  sprintf("Full model: %d features. Intrinsic model: %d features (computed from each mRNA's own sequence/structure). The %d external features are measured data from other experiments and are excluded from the intrinsic model.", length(full), length(intr), length(full) - length(intr)),
  "Membership is read from the saved models (si3d hypoxia 1hr, full and intrinsic); every heatmap column uses the same feature set, asserted in the modelling code.",
  "All features are computed on each dataset's own dominant transcript where available (MDA-MB-231, MCF7-SIX1; HeLa and HEK293T in progress).",
  "Lengths are log2-transformed; raw lengths are not modelled.",
  "Structure: 5'UTR accessibility uses a smaller RNAplfold window (80/40) than CDS and 3'UTR (150/100), so 5'UTR values are not directly comparable to the other regions.",
  "G4 'max' features are adjusted for region length and GC so they measure G4 potential beyond composition.",
  "Codon features (CAI, FOP, tAI, CSC, GC3) exclude the start and stop codons; CAI/FOP weights are shared across cell lines so values are comparable."))

wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#E8E8E8", border = "Bottom", wrapText = TRUE, valign = "top")
ws <- createStyle(wrapText = TRUE, valign = "top")
for (nm in c("Summary by family", "All features", "Notes")) addWorksheet(wb, nm)
writeData(wb, "Summary by family", summ, headerStyle = hs)
setColWidths(wb, "Summary by family", 1:4, c(30, 60, 12, 15))
addStyle(wb, "Summary by family", createStyle(textDecoration = "bold"), rows = nrow(summ) + 1, cols = 1:4, gridExpand = TRUE)
writeData(wb, "All features", tab, headerStyle = hs)
setColWidths(wb, "All features", 1:8, c(26, 42, 80, 14, 12, 14, 26, 50))
addStyle(wb, "All features", ws, rows = 2:(nrow(tab) + 1), cols = 1:8, gridExpand = TRUE)
freezePane(wb, "All features", firstRow = TRUE)
writeData(wb, "Notes", notes, headerStyle = hs); setColWidths(wb, "Notes", 1, 120)
addStyle(wb, "Notes", ws, rows = 2:(nrow(notes) + 1), cols = 1)
out <- file.path(o, "feature_table_for_slides.xlsx")
saveWorkbook(wb, out, overwrite = TRUE)
write_csv(tab, file.path(o, "feature_table_for_slides.csv"))
cat("wrote", out, "\n"); print(as.data.frame(summ), row.names = FALSE)
