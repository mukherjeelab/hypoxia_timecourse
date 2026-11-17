# One dot per mRNA: si3d TE LFC vs a sequence feature, normoxia 1+4hr and hypoxia 1hr side by side.
# Same layout as 54_te_lfc_gc3_scatter.pdf, extended to GC2, GC3, CDS length and UTR lengths.
# Features come from feature_matrix_positional_gc.rds (01v_), 10,172 transcripts.
# Run: Rscript code/predictive_modeling/54b_te_lfc_feature_scatters.R
suppressPackageStartupMessages({library(tidyverse); library(here)})

fm <- readRDS(here("output", "predictive_modeling", "feature_matrix_positional_gc.rds"))
stopifnot("matrix has duplicate genes" = !anyDuplicated(fm$gene_id_clean))

read_te <- function(f, label) {
  d <- read_csv(here("output", f), show_col_types = FALSE) %>%
    mutate(gene_id_clean = sub("\\..*", "", gene_id))
  # Columns are log2(sictrl / si3d) despite the filename; EIF3D's own rna_lfc must be strongly positive.
  eif3d <- d$rna_lfc[d$gene_id_clean == "ENSG00000100353"]
  stopifnot("EIF3D missing" = length(eif3d) == 1, "contrast direction flipped" = eif3d > 2)
  d %>% transmute(gene_id_clean, te_lfc, condition = label)
}
te <- bind_rows(
  read_te("translation_categories_si3d_vs_sictrl_normoxia_1and4hr.csv", "Normoxia"),
  read_te("translation_categories_si3d_vs_sictrl_hypoxia_1hr.csv", "Hypoxia - 1 hr")
) %>% mutate(condition = factor(condition, levels = c("Normoxia", "Hypoxia - 1 hr")))

features <- list(
  gc1  = list(col = "gc1",              x = "GC1, first codon position GC content (%)",  pct = TRUE),
  gc2 = list(col = "gc2",              x = "GC2, second codon position GC content (%)", pct = TRUE),
  gc3  = list(col = "gc3",              x = "GC3, third codon position GC content (%)",  pct = TRUE),
  cds  = list(col = "log2_cds_length",  x = "log2 CDS length (nt)",                     pct = FALSE),
  utr5 = list(col = "log2_utr5_length", x = "log2 5' UTR length (nt)",                  pct = FALSE),
  utr3 = list(col = "log2_utr3_length", x = "log2 3' UTR length (nt)",                  pct = FALSE)
)

for (nm in names(features)) {
  f <- features[[nm]]
  d <- fm %>% transmute(gene_id_clean, x = .data[[f$col]]) %>%
    inner_join(te, by = "gene_id_clean") %>% filter(!is.na(x), !is.na(te_lfc))
  if (f$pct) stopifnot("GC out of [0,100]" = all(d$x >= 0 & d$x <= 100))
  stopifnot("non-finite values" = all(is.finite(d$x) & is.finite(d$te_lfc)))
  st <- d %>% group_by(condition) %>%
    summarise(rho = round(cor(x, te_lfc, method = "spearman"), 2), n = n(), .groups = "drop")
  cat(nm, ":", paste(sprintf("%s rho=%.2f n=%d", st$condition, st$rho, st$n), collapse = " | "), "\n")
  p <- ggplot(d, aes(x, te_lfc)) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.5) +
    geom_point(size = 0.5, alpha = 0.2, colour = "grey40") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#E69F00", linewidth = 0.9) +
    geom_text(data = st, aes(x = Inf, y = Inf, label = paste0("rho=", rho, "\nn=", n)),
              colour = "black", hjust = 1.1, vjust = 1.3, size = 4.5, inherit.aes = FALSE) +
    facet_wrap(~ condition, nrow = 1) +
    labs(x = f$x, y = "TE LFC (positive = eIF3d promotes)") +
    theme_classic(base_size = 16) +
    coord_cartesian(ylim = c(-1, 2))  # not ylim(): that would drop rows before the lm fit
  ggsave(here("plots", sprintf("54b_te_lfc_%s_scatter_mda231.pdf", nm)), p, width = 8, height = 4)
}
