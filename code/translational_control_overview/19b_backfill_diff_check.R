# 19b_backfill_diff_check.R
# After the 2026-10-01 backfill rebuild, assert that feature_matrix_tco.rds differs from the
# pre-backfill copy ONLY where intended:
#   - G4 columns, hia2026_dhx29_occupancy_tx and tco_struct_accessibility_* may gain values
#     (NA -> value) but never lose or change one;
#   - g4mer_max_resid_* may also shift on existing values (coefficients refit on more rows);
#   - cnot3_weighted_codon_score_tx is new;
#   - every other column is identical.
suppressMessages({ library(dplyr); library(here) })
o <- here("output", "translational_control_overview")
old <- readRDS(file.path(o, "_pre_backfill_20261001", "feature_matrix_tco.rds"))
new <- readRDS(file.path(o, "feature_matrix_tco.rds"))
stopifnot("transcript set or order changed" = identical(old$transcript_id_clean, new$transcript_id_clean),
          "unexpected new columns" = setequal(setdiff(names(new), names(old)), "cnot3_weighted_codon_score_tx"),
          "columns were removed" = all(names(old) %in% names(new)))
fill_ok  <- function(f) grepl("^g4mer_|^hia2026_dhx29_occupancy_tx$|^tco_struct_accessibility_", f)
shift_ok <- function(f) grepl("^g4mer_max_resid_", f)
r <- bind_rows(lapply(names(old), function(f) {
  u <- old[[f]]; v <- new[[f]]
  if (!is.numeric(u)) return(tibble(feature = f, filled = 0L, emptied = 0L, changed = as.integer(!identical(u, v)), maxdev = NA_real_))
  k <- !is.na(u) & !is.na(v)
  tibble(feature = f, filled = sum(is.na(u) & !is.na(v)), emptied = sum(!is.na(u) & is.na(v)),
         changed = sum(abs(u[k] - v[k]) > 1e-12), maxdev = if (any(k)) max(abs(u[k] - v[k])) else 0)
}))
moved <- r %>% filter(filled > 0 | emptied > 0 | changed > 0)
print(as.data.frame(moved), row.names = FALSE)
bad <- moved %>% filter(emptied > 0 | (filled > 0 & !fill_ok(feature)) | (changed > 0 & !shift_ok(feature)))
if (nrow(bad)) { cat("UNEXPECTED CHANGES:\n"); print(as.data.frame(bad), row.names = FALSE) }
stopifnot("a column changed that the backfill should not touch" = nrow(bad) == 0)
cat("backfill diff check PASSED:", nrow(moved), "columns changed, all as intended\n")
