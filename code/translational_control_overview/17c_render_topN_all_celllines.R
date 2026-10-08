# 17c_render_topN_all_celllines.R
# The top-N mean-feature-value heatmap (17_) for exactly the columns of 15b_'s intrinsic figures, at
# any N. The experiment list is taken from 15b_ itself (its T / mk / FOOT definitions are evaluated
# from that file, nothing is copied), so the two figures cannot drift apart.
#   Rscript code/translational_control_overview/17c_render_topN_all_celllines.R <top_n> [capbind]
a <- commandArgs(TRUE); top_n <- as.integer(a[1]); with_capbind <- isTRUE(a[2] == "capbind")
stopifnot("usage: 17c_... <top_n> [capbind]" = !is.na(top_n))
src <- parse(here::here("code/translational_control_overview/15b_render_all_celllines.R"))
keep <- vapply(src, function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                 as.character(e[[2]]) %in% c("T", "mk", "FOOT"), logical(1))
for (e in src[keep]) eval(e)
ex <- mk("_intrinsic")
tag <- if (with_capbind) "all_celllines_intrinsic_capbind" else "all_celllines_intrinsic"
rmarkdown::render(here::here("code/translational_control_overview/17_top100_feature_means_heatmap.Rmd"),
  params = list(experiments = ex, out_tag = tag, top_n = top_n,
                footnote = if (grepl("*", ex, fixed = TRUE)) FOOT else ""),
  output_file = paste0("17_top", top_n, "_feature_means_heatmap_", tag, ".html"),
  intermediates_dir = tempfile(), quiet = TRUE, envir = new.env())
