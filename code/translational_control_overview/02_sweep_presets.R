# Sweep presets shared by 02b_neg_subsample_sweep.R (runs them) and
# 02c_sweep_stability.Rmd (reports them).
#
# Defined in one place because a sweep is only interpretable if the arms that produced it are
# still on disk. Previously 02b's CONFIGS was overwritten by each new question, so the
# termination arms survived only in one commit, the peptide arms in another - every result in
# CLAUDE.md needed a different checkout to reproduce. Presets sit side by side instead.
#
# Every arm in a preset is rendered at the same 20 neg_seeds. neg_seed seeds ONLY the negative
# draw - the split, the forest and the CV stay pinned at 9 - so within a seed every arm sees
# identical genes in an identical split, and paired differences isolate the arm alone.

# A family sweep is always the same three arms. The shuffled arm is the point: impurity
# importance is non-negative, so an uninformative family never scores zero, and whatever the
# permuted copy still absorbs is the floor the real family has to clear. Permuting preserves
# every marginal, the NA pattern and the within-family correlation structure.
family_arms <- function(families, stem) {
  list(
    none = list(feature_set_label = paste0(stem, "none"), exclude_families = families),
    real = list(feature_set_label = paste0(stem, "real")),
    shuf = list(feature_set_label = paste0(stem, "shuf"), shuffle_families = families)
  )
}

SWEEP_PRESETS <- list(
  # Definition corrections: one member of each variant pair (00_'s FIX rows).
  variants = list(
    current   = list(feature_set_label = "tco", variant_set = "current"),
    corrected = list(feature_set_label = "tco", variant_set = "corrected")
  ),
  # Two readings of start-context strength. Kept as a preset because the inventory's
  # "never model both together" turned out to be wrong - they correlate at only 0.15.
  initiation = list(
    kozak_pwm   = list(feature_set_label = "tco_init", initiation_score = "pwm"),
    noderer_tis = list(feature_set_label = "tco_init", initiation_score = "noderer")
  ),
  termination = family_arms("termination", "tco_term"),
  peptide     = family_arms("nascent_peptide,polya_track", "tco_pep"),
  g4          = family_arms("g4", "tco_g4")
)

# Which family each preset is about, for the share/rank reporting in 02c. NA = not a family
# sweep, so only the AUC comparison applies.
PRESET_FAMILIES <- list(
  variants    = NA_character_,
  initiation  = NA_character_,
  termination = "termination",
  peptide     = c("nascent_peptide", "polya_track"),
  g4          = "g4"
)

SWEEP_SEEDS <- 1:20

# Mirrors the suffix built in 02_rf_model.Rmd. If that construction changes, change it here
# too, or 02c silently finds nothing and reports an empty sweep.
sweep_suffix <- function(arm, seed) {
  p <- utils::modifyList(
    list(gene = "3d", direction = "promotes", readout = "TE", condition = "hypoxia",
         timepoint = "1hr", lfc = "0.5", feature_set_label = "tco",
         variant_set = "corrected", initiation_score = "both",
         clip_mode = "aggregate", neg_tag = ""),
    arm)
  paste0(
    "eif3", p$gene, "_", p$direction, "_",
    if (p$readout != "TE") paste0(tolower(p$readout), "_") else "",
    p$condition, "_", p$timepoint, "_lfc", p$lfc,
    "_", p$feature_set_label,
    if (p$variant_set != "current") paste0("_", p$variant_set) else "",
    if (p$initiation_score != "pwm") paste0("_", p$initiation_score) else "",
    "_clip", p$clip_mode,
    if (nzchar(p$neg_tag)) paste0("_", p$neg_tag) else "",
    if (seed != 9L) paste0("_negseed", seed) else ""
  )
}
