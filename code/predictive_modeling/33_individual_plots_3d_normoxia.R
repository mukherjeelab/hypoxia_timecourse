library(tidyverse)
library(here)

# Load feature matrix + log2 lengths
features <- readRDS(here("output", "predictive_modeling", "feature_matrix_external_stability.rds")) %>%
  mutate(
    log2_utr5_length = log2(utr5_length + 1),
    log2_cds_length  = log2(cds_length  + 1),
    log2_utr3_length = log2(utr3_length + 1)
  )
cat("Feature matrix:", nrow(features), "genes x", ncol(features), "columns\n")

# --- Gene sets ---
d_hyp1_ids <- read_csv(
  here("output", "genesets", "hypoxia_3d_promotes_TE_1hr_lfc0.5.csv"),
  show_col_types = FALSE
) %>% pull(ensembl_gene)

d_norm14_ids <- read_csv(
  here("output", "translation_categories_si3d_vs_sictrl_normoxia_1and4hr.csv"),
  show_col_types = FALSE
) %>%
  filter(te_lfc > 0.5, te_padj < 0.05) %>%
  mutate(ensembl_gene = sub("\\..*", "", gene_id)) %>%
  pull(ensembl_gene)

e_hyp1_ids <- read_csv(
  here("output", "genesets", "hypoxia_3e_promotes_TE_1hr_lfc0.5.csv"),
  show_col_types = FALSE
) %>% pull(ensembl_gene)

e_norm14_ids <- read_csv(
  here("output", "te_si3e_vs_sictrl_normoxia_1and4hr.csv"),
  show_col_types = FALSE
) %>%
  filter(te_lfc > 0.5, te_padj < 0.05) %>%
  mutate(ensembl_gene = sub("\\..*", "", gene_id)) %>%
  pull(ensembl_gene)

# --- Negative controls ---
neg_pool_3d <- read_csv(
  here("output", "predictive_modeling", "negative_control_genes.csv"),
  show_col_types = FALSE
) %>%
  mutate(gene_id_clean = sub("\\..*", "", gene_id)) %>%
  filter(gene_id_clean %in% features$gene_id_clean)

neg_pool_3e <- read_csv(
  here("output", "predictive_modeling", "negative_control_genes_si3e.csv"),
  show_col_types = FALSE
) %>%
  mutate(gene_id_clean = sub("\\..*", "", gene_id)) %>%
  filter(gene_id_clean %in% features$gene_id_clean)

set.seed(9); neg_d_hyp1   <- sample(neg_pool_3d$gene_id_clean, size = min(length(d_hyp1_ids),   nrow(neg_pool_3d)))
set.seed(9); neg_d_norm14 <- sample(neg_pool_3d$gene_id_clean, size = min(length(d_norm14_ids), nrow(neg_pool_3d)))
set.seed(9); neg_e_hyp1   <- sample(neg_pool_3e$gene_id_clean, size = min(length(e_hyp1_ids),   nrow(neg_pool_3e)))
set.seed(9); neg_e_norm14 <- sample(neg_pool_3e$gene_id_clean, size = min(length(e_norm14_ids), nrow(neg_pool_3e)))

# --- Density plot helper ---
plot_density <- function(all_df, promoted_ids, neg_ids, col, xlabel, title_str) {
  bg   <- all_df %>% filter(!is.na(.data[[col]]))
  fg_p <- bg %>% filter(gene_id_clean %in% promoted_ids)
  fg_n <- bg %>% filter(gene_id_clean %in% neg_ids)
  ks   <- ks.test(
    fg_p[[col]],
    bg %>% filter(!gene_id_clean %in% promoted_ids) %>% pull(col)
  )
  ggplot() +
    geom_density(data = bg,   aes(x = .data[[col]]), color = "grey60",  fill = NA, linewidth = 0.9) +
    geom_density(data = fg_n, aes(x = .data[[col]]), color = "#000000", fill = NA, linewidth = 0.9) +
    geom_density(data = fg_p, aes(x = .data[[col]]), color = "#E69F00", fill = NA, linewidth = 0.9) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = paste0("KS p = ", formatC(ks$p.value, format = "e", digits = 2)),
             size = 3.5) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 3.2,
             label = paste0("All genes (n = ", nrow(bg), ")"),
             size = 3.5, color = "grey50") +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 4.9,
             label = paste0("Promoted (n = ", nrow(fg_p), ")"),
             size = 3.5, color = "#E69F00") +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 6.6,
             label = paste0("Neg ctrl (n = ", nrow(fg_n), ")"),
             size = 3.5, color = "#000000") +
    labs(title = title_str, x = xlabel, y = "Density") +
    theme_classic(base_size = 13)
}

# --- Conditions ---
conditions <- list(
  list(ids = d_hyp1_ids,   neg = neg_d_hyp1,   label = "eIF3d hypoxia 1 hr",      prefix = "33_3d_hyp1"),
  list(ids = d_norm14_ids, neg = neg_d_norm14,  label = "eIF3d normoxia 1+4 hr",   prefix = "33_3d_norm14"),
  list(ids = e_hyp1_ids,   neg = neg_e_hyp1,   label = "eIF3e hypoxia 1 hr",      prefix = "33_3e_hyp1"),
  list(ids = e_norm14_ids, neg = neg_e_norm14,  label = "eIF3e normoxia 1+4 hr",   prefix = "33_3e_norm14")
)

# --- Features ---
features_specs <- list(
  # list(col = "log2_utr5_length", xlabel = "5' UTR length (log2 nt)", suffix = "utr5_length", title = "5' UTR length"),
  list(col = "log2_cds_length",  xlabel = "CDS length (log2 nt)",    suffix = "cds_length",  title = "CDS length"),
  # list(col = "log2_utr3_length", xlabel = "3' UTR length (log2 nt)", suffix = "utr3_length", title = "3' UTR length"),
  # list(col = "utr5_gc",          xlabel = "5' UTR GC content (%)",   suffix = "utr5_gc",     title = "5' UTR GC"),
  list(col = "cds_gc",           xlabel = "CDS GC content (%)",      suffix = "cds_gc",      title = "CDS GC")
  # list(col = "utr3_gc",          xlabel = "3' UTR GC content (%)",   suffix = "utr3_gc",     title = "3' UTR GC")
)

# --- Generate and save all 24 plots ---
for (cond in conditions) {
  for (feat in features_specs) {
    fname <- paste0(cond$prefix, "_", feat$suffix)
    title_str <- paste0(feat$title, "\n", cond$label)
    cat("Saving:", fname, "\n")
    g <- plot_density(features, cond$ids, cond$neg, feat$col, feat$xlabel, title_str)
    ggsave(here("plots", paste0(fname, ".pdf")), g, width = 5, height = 4.5)
  }
}

cat("Done. 24 PDFs saved to plots/\n")
