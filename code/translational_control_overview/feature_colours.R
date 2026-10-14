# feature_colours.R - ONE feature -> colour-category mapping for every SHAP figure in this fork,
# so a feature has the same colour in the bar plots (03i_), the beeswarm (03j_) and any later plot.
#   source(here("code", "translational_control_overview", "feature_colours.R"))
#   feature_category(features, ann)  -> factor with levels names(FEATURE_COLOURS)
#
# Categories are by mRNA REGION for sequence features, and by biology for external datasets.
# Decided with Kate 2026-10-09:
#   - start context (Kozak, Noderer TIS) goes with the 5'UTR; termination has its own colour
#   - codon-weighted scores (CNOT3 codon score, DHX29 occupancy = codon-weighted CES over the CDS,
#     23_dhx29_riboseq_slamseq.Rmd) and CNOT3 ribosome enrichment go with the CDS
#   - half-lives, Karner SLAM-seq, CNOT3-KO half-life change and DHX29-KD RNA change are mRNA stability
#   - no category is grey: grey is reserved for the low end of the beeswarm's value gradient
# Okabe-Ito where possible. Orange and sky blue stay reserved for eIF3d / eIF3e, so termination,
# eIF3 binding and mRNA stability use Paul Tol colours; eIF3 binding and stability are darkened
# from Tol's olive / brown so they stay apart from CDS vermillion under deuteranopia (worst pair
# CIEDE2000 11.4 vs 5.1 with standard olive, colorspace::deutan simulation).
#
# feature_category() STOPS on any feature it cannot place, so a new family cannot silently
# inherit a default colour. Extend the rules below when the matrix gains one.

FEATURE_COLOURS <- c(
  "5' UTR"                    = "#009E73",  # Okabe-Ito bluish green
  "CDS"                       = "#D55E00",  # Okabe-Ito vermillion
  "Termination"               = "#882255",  # Tol wine
  "3' UTR"                    = "#F0E442",  # Okabe-Ito yellow
  "Full transcript"           = "#332288",  # Tol indigo
  "mRNA stability"            = "#5C3A0A",  # dark brown
  "eIF3 binding (CLIP)"       = "#666633",  # dark olive
  "eIF4E"                     = "#0072B2",  # Okabe-Ito blue
  "Other translation factors" = "#CC79A7"   # Okabe-Ito reddish purple
)

# low end of the beeswarm gradient per category: light grey, except yellow, whose brightness
# matches light grey - it starts from a darker grey so the gradient is visible
FEATURE_GRADIENT_LOW <- setNames(rep("grey85", length(FEATURE_COLOURS)), names(FEATURE_COLOURS))
FEATURE_GRADIENT_LOW[["3' UTR"]] <- "grey45"

# the one family split by feature rather than by family: external stability datasets
STABILITY_TO_CDS <- c("cnot3_ribo_enrichment", "cnot3_weighted_codon_score", "cnot3_weighted_codon_score_tx",
                      "hia2026_dhx29_occupancy", "hia2026_dhx29_occupancy_tx")
STABILITY_TO_STABILITY <- c("zhu2024_ctrl_halflife", "zhu2024_cnot3ko_halflife_logfc", "hia2026_ctrl_halflife",
                            "hia2026_dhx29_halflife_logfc", "karner2026_mda231_log2_ct", "hia2026_dhx29_kd_rna_lfc")

feature_category <- function(features, ann) {
  fam <- ann$family[match(features, ann$feature)]
  blk <- ann$block[match(features, ann$feature)]
  if (anyNA(fam)) stop("feature(s) missing from feature_annotation.csv: ",
                       paste(features[is.na(fam)], collapse = ", "))
  region <- dplyr::case_when(grepl("utr5", features) ~ "5' UTR",
                             grepl("utr3", features) ~ "3' UTR",
                             grepl("cds", features)  ~ "CDS",
                             TRUE ~ NA_character_)
  cat <- dplyr::case_when(
    blk == "eif4e"                                            ~ "eIF4E",
    fam == "clip"                                             ~ "eIF3 binding (CLIP)",
    fam == "dap5"                                             ~ "Other translation factors",
    fam == "stability_external" & features %in% STABILITY_TO_CDS       ~ "CDS",
    fam == "stability_external" & features %in% STABILITY_TO_STABILITY ~ "mRNA stability",
    fam == "termination"                                      ~ "Termination",
    fam %in% c("initiation", "structure_utr5")                ~ "5' UTR",
    fam == "structure_utr3"                                   ~ "3' UTR",
    fam %in% c("structure_cds", "codon_optimality", "nascent_peptide", "polya_track") ~ "CDS",
    fam == "g4"                                               ~ region,
    fam %in% c("gc", "length") & grepl("^(log2_)?total_length$|^transcript_gc$", features) ~ "Full transcript",
    fam == "gc" & grepl("^gc3", features)                     ~ "CDS",
    fam %in% c("gc", "length")                                ~ region,
    TRUE                                                      ~ NA_character_)
  if (anyNA(cat)) stop("no colour category for: ", paste(features[is.na(cat)], collapse = ", "),
                       " - add a rule to feature_colours.R")
  factor(cat, levels = names(FEATURE_COLOURS))
}
