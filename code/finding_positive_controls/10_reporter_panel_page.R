# Renders output/finding_positive_controls/08_reporter_panel_table.csv into a shareable HTML page
# using 10_reporter_panel_page_template.tpl (named .tpl because the repo gitignores *.html). Every value comes from the CSV, none is typed here.
# Usage: Rscript code/finding_positive_controls/10_reporter_panel_page.R
# Writes output/finding_positive_controls/reporter_panel.html (publish that file as the artifact).
suppressPackageStartupMessages({library(tidyverse); library(here)})
d <- read_csv(here("output", "finding_positive_controls", "08_reporter_panel_table.csv"), show_col_types = FALSE)
esc <- function(x) { x <- coalesce(as.character(x), ""); x <- gsub("&", "&amp;", x); x <- gsub("<", "&lt;", x); gsub(">", "&gt;", x) }
# "1.52 (6.5e-09)" -> value + small padj; bold when padj < 0.05
cell <- function(x) {
  m <- str_match(x, "^(-?[0-9.]+) \\((.+)\\)$")
  sig <- suppressWarnings(as.numeric(m[, 3])) < 0.05
  ifelse(is.na(m[, 1]), esc(x),
         sprintf('<span class="v%s">%s</span><span class="p">%s</span>', ifelse(sig, " sig", ""), m[, 2], m[, 3]))
}
cat_pill <- function(x) sprintf('<span class="cat cat-%s">%s</span>', tolower(x), esc(gsub("_", " ", x)))
role_cls <- c("selected eIF3d target" = "target", "positive control (Roiuk et al)" = "pos", "negative control" = "neg")
role_lab <- c("selected eIF3d target" = "Target", "positive control (Roiuk et al)" = "Positive control", "negative control" = "Negative control")
stopifnot(all(d$role %in% names(role_cls)))
gene <- function(i) sprintf('<th scope="row"><span class="dot %s"></span>%s<span class="role">%s</span></th>',
                            role_cls[d$role[i]], esc(d$symbol[i]), role_lab[d$role[i]])
t1 <- map_chr(seq_len(nrow(d)), function(i) paste0("<tr>", gene(i),
  "<td>", cell(d$te_231[i]), "</td><td>", cell(d$ribo_231[i]), "</td><td>", cell(d$rna_231[i]), "</td><td>", cat_pill(d$category_231[i]), "</td>",
  '<td class="split">', cell(d$te_mcf[i]), "</td><td>", cell(d$ribo_mcf[i]), "</td><td>", cell(d$rna_mcf[i]), "</td><td>", cat_pill(d$category_mcf[i]), "</td></tr>"))
t4 <- map_chr(seq_len(nrow(d)), function(i) paste0("<tr>", gene(i),
  "<td>", cell(d$si3e_te_231[i]), "</td><td>", cell(d$si3e_ribo_231[i]), "</td><td>", cell(d$si3e_rna_231[i]), "</td><td>", cat_pill(d$si3e_category_231[i]), "</td>",
  '<td class="split">', cell(d$si3e_te_mcf[i]), "</td><td>", cell(d$si3e_ribo_mcf[i]), "</td><td>", cell(d$si3e_rna_mcf[i]), "</td><td>", cat_pill(d$si3e_category_mcf[i]), "</td>",
  '<td class="split tgt"><span class="call">', esc(sub("^target", "Target", sub("^not", "Not", d$eif3e_target[i]))),
  '</span><span class="rna">', esc(sub("^mRNA after si3e: ", "mRNA: ", d$eif3e_rna_change[i])), "</span></td></tr>"))
t2 <- map_chr(seq_len(nrow(d)), function(i) paste0("<tr>", gene(i),
  "<td>", cell(d$te_nor_231[i]), "</td><td>", cell(d$te_nor_mcf[i]), "</td></tr>"))
t5 <- map_chr(seq_len(nrow(d)), function(i) paste0("<tr>", gene(i),
  "<td>", cell(d$te_hela_si3d[i]), "</td><td>", cell(d$te_hela_si3e[i]), "</td>",
  '<td class="split tgt"><span class="call">', esc(d$hela_target[i]), "</span></td></tr>"))
nc <- function(x) ifelse(grepl("^not called$|^0 of", x), sprintf('<span class="muted">%s</span>', esc(x)), esc(x))
t3 <- map_chr(seq_len(nrow(d)), function(i) paste0("<tr>", gene(i),
  "<td>", nc(d$cap_binding_subunitseq[i]), "</td><td>", nc(d$lee2015_parclip[i]), "</td><td>", nc(d$desilva_eif3d_band[i]),
  "</td><td>", nc(d$desilva_other_bands[i]), "</td><td>", nc(d$cate_irclip[i]), "</td></tr>"))
notes <- map_chr(seq_len(nrow(d)), function(i) sprintf('<div class="note"><h3><span class="dot %s"></span>%s</h3><p>%s</p></div>',
                                                        role_cls[d$role[i]], esc(d$symbol[i]), esc(d$notes[i])))
tpl <- readLines(here("code", "finding_positive_controls", "10_reporter_panel_page_template.tpl"))
out <- paste(tpl, collapse = "\n")
out <- sub("%%T1%%", paste(t1, collapse = "\n"), out, fixed = TRUE)
out <- sub("%%T2%%", paste(t2, collapse = "\n"), out, fixed = TRUE)
out <- sub("%%T4%%", paste(t4, collapse = "\n"), out, fixed = TRUE)
out <- sub("%%T5%%", paste(t5, collapse = "\n"), out, fixed = TRUE)
out <- sub("%%T3%%", paste(t3, collapse = "\n"), out, fixed = TRUE)
out <- sub("%%NOTES%%", paste(notes, collapse = "\n"), out, fixed = TRUE)
stopifnot(!grepl("%%", out))
out_file <- here("output", "finding_positive_controls", "reporter_panel.html")
writeLines(out, out_file)
cat("wrote", out_file, "\n")
