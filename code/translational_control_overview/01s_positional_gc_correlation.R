# 01s_positional_gc_correlation.R - Spearman correlation of GC1 / GC2 / GC3 (internal codons,
# from 01r_) in three gene groups: the si3d eIF3d-promoted set, a size-matched negative draw
# (seed 9), and every gene in the matrix. MDA-MB-231 hypoxia/normoxia and MCF7-SIX1
# hypoxia/normoxia (each on its own-isoform matrix).
# Normoxia is MDA-MB-231's 1+4hr time-regressed contrast (never plain 1hr). Its gene sets come
# from the tco4e run: the positive set and the seed-9 negative draw depend only on the matrix
# rows, which feature_matrix_tco_posgc.rds shares with feature_matrix_tco.rds.
#   Rscript code/translational_control_overview/01s_positional_gc_correlation.R [hypoxia|normoxia|mcf7six1_hypoxia|mcf7six1_normoxia]
suppressPackageStartupMessages({ library(tidyverse); library(here) })

O   <- here("output", "translational_control_overview")
args <- commandArgs(trailingOnly = TRUE)
condition <- if (length(args) >= 1) args[1] else "hypoxia"
# MDA-MB-231 sets come from the classifier's saved gene sets (positives + its seed-9 negative
# draw). MCF7-SIX1 has no classifier run, so its negative draw is rebuilt the way 02_ does it:
# negative set restricted to the matrix, positives removed, set.seed(9), size-matched sample.
# MCF7-SIX1 has only 1hr; its gene-set te_lfc is already negated to promotes-positive.
cfg <- switch(condition,
  hypoxia  = list(cell = "mda231", ctp = "hypoxia_1hr", label = "tco4eposgc", name = "hypoxia",
                  desc = "MDA-MB-231 si3d hypoxia 1hr", mat = "feature_matrix_tco_posgc.rds"),
  normoxia = list(cell = "mda231", ctp = "normoxia_1and4hr", label = "tco4e", name = "normoxia",
                  desc = "MDA-MB-231 si3d normoxia 1+4hr (time-regressed)",
                  mat = "feature_matrix_tco_posgc.rds"),
  mcf7six1_hypoxia  = list(cell = "mcf7six1", ctp = "mcf7six1_hypoxia_1hr", name = "hypoxia",
                           cond = "hypoxia", desc = "MCF7-SIX1 si3d hypoxia 1hr",
                           mat = "feature_matrix_tco_mcf7six1_posgc.rds"),
  mcf7six1_normoxia = list(cell = "mcf7six1", ctp = "mcf7six1_normoxia_1hr", name = "normoxia",
                           cond = "normoxia", desc = "MCF7-SIX1 si3d normoxia 1hr",
                           mat = "feature_matrix_tco_mcf7six1_posgc.rds"),
  stop("condition must be hypoxia, normoxia, mcf7six1_hypoxia or mcf7six1_normoxia"))
m    <- readRDS(file.path(O, cfg$mat))
cols <- c(GC1 = "gc1_internal", GC2 = "gc2_internal", GC3 = "gc3_internal")
stopifnot(all(cols %in% names(m)))

if (cfg$cell == "mda231") {
  gs <- read_csv(file.path(O, paste0("rf_model_gene_sets_eif33d_promotes_", cfg$ctp, "_lfc0.5_",
                                     cfg$label, "_corrected_both_clipaggregate.csv")),
                 show_col_types = FALSE)
  stopifnot(!anyDuplicated(gs$gene_id_clean),
            setequal(unique(gs$eif3d_target), c("promotes", "does_not_promote")))
  pos_ids <- gs$gene_id_clean[gs$eif3d_target == "promotes"]
  neg_ids <- gs$gene_id_clean[gs$eif3d_target == "does_not_promote"]
} else {
  G   <- here("output", "genesets")
  pos <- read_csv(file.path(G, paste0("mcf7six1_", cfg$cond, "_3d_promotes_TE_1hr_lfc0.5.csv")),
                  show_col_types = FALSE)
  neg <- read_csv(file.path(G, paste0("mcf7six1_", cfg$cond, "_3d_negative_controls_lfc0.5.csv")),
                  show_col_types = FALSE) %>% mutate(gene_id_clean = sub("\\..*", "", gene_id))
  stopifnot("positive set is not promotes-positive" = all(pos$te_lfc > 0.5),
            "negative set is not flat" = all(abs(neg$te_lfc) < 0.5))
  pos_ids  <- intersect(pos$ensembl_gene, m$gene_id_clean)
  neg_pool <- sort(setdiff(intersect(neg$gene_id_clean, m$gene_id_clean), pos_ids))
  stopifnot("negative pool smaller than the positive set" = length(neg_pool) > length(pos_ids))
  set.seed(9)
  neg_ids <- sample(neg_pool, length(pos_ids))
  cat(sprintf("negative draw: %d of a %d-gene pool (seed 9)\n", length(neg_ids), length(neg_pool)))
}

# matrix is transcript-keyed; one row per gene
genes <- m %>% distinct(gene_id_clean, .keep_all = TRUE) %>%
  dplyr::select(gene_id_clean, all_of(cols)) %>% drop_na()
neg_name <- if (cfg$cell == "mda231") "Negative draw (model, seed 9)" else "Negative draw (size-matched, seed 9)"
groups <- rlang::list2(
  !!paste("eIF3d", cfg$name, "promotes") := genes %>% filter(gene_id_clean %in% pos_ids),
  !!neg_name                              := genes %>% filter(gene_id_clean %in% neg_ids),
  "All genes in matrix"                   = genes)

cm <- imap_dfr(groups, function(d, g) {
  r <- cor(d[names(cols)], method = "spearman")
  as_tibble(r, rownames = "x") %>% pivot_longer(-x, names_to = "y", values_to = "rho") %>%
    mutate(group = sprintf("%s\nn = %d", g, nrow(d)))
}) %>% mutate(group = fct_inorder(group), y = factor(y, levels = rev(names(cols))))
print(cm %>% filter(x < as.character(y)) %>% mutate(rho = round(rho, 3)) %>% arrange(group, x, y), n = Inf)

p <- ggplot(cm, aes(x, y, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 4.2) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
                       limits = c(-1, 1), name = "Spearman rho") +
  facet_wrap(~ group) + coord_fixed() +
  labs(x = NULL, y = NULL,
       title = paste0("GC content by codon position: Spearman correlation, ", cfg$desc),
       subtitle = "Internal codons (initiator and stop removed); one transcript per gene",
       caption = paste0("eIF3d-promoted: TE promotes-positive LFC > 0.5, padj < 0.05. Negatives: ",
                        if (cfg$cell == "mda231") "the size-matched draw the classifier used"
                        else "size-matched seed-9 draw from the |TE LFC| < 0.5 negative set")) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"),
        plot.caption.position = "plot", plot.title.position = "plot")
out <- here("plots", paste0("01s_positional_gc_correlation_si3d_", cfg$ctp, ".pdf"))
ggsave(out, p, width = 10, height = 4.4)
cat("wrote", out, "\n")
