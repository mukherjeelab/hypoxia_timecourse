#!/usr/bin/env Rscript
# Render the C-terminal HiBiT ladder over the 1-4x eGFP range.
#
#   Rscript code/cds_reporter_design/16_render_1to4x.R
#
# This is a driver, not a second copy of `15_`. The two ranges differ only in
# which rungs are built, so a duplicated notebook would be ~600 lines that have
# to be kept in step by hand - the same trap `02b_`'s overwritten config list
# fell into in the overview fork. Same pattern as `39_render_contrasts.R` and
# `26b_render_sets.R`.
#
# Why 1-4x rather than 0-3x. The no-GFP rung is 96 nt: below 99.9% of expressed
# MDA-MB-231 mRNAs (only 13 of 10,172 genes are shorter) and, at typical
# ribosome spacing, holding under one ribosome at a time. It is not a shorter
# version of the same process, it is a different regime, and its 32-aa
# unstructured product is a plausible degradation substrate that would confound
# a HiBiT readout. Dropping it and adding a 4th eGFP copy gives
# 855 / 1,614 / 2,373 / 3,132 nt - the 31st to 89th percentile - with the 1x
# rung landing essentially on the eIF3d-promoted median of 891 nt.

suppressMessages(library(rmarkdown))
suppressMessages(library(here))

nb  <- here("code", "cds_reporter_design", "15_cterm_hibit_ladder.Rmd")
stopifnot("15_cterm_hibit_ladder.Rmd not found" = file.exists(nb))

out <- "15_cterm_hibit_ladder_1to4x.html"

render(nb,
       output_file = out,
       params = list(min_gfp_copies = 1, max_gfp_copies = 4),
       envir  = new.env())

cat("\nRendered", out, "\n")
cat("Outputs and plots carry the run tag `1to4x_expressed`,",
    "so the 0-3x set is untouched.\n")
