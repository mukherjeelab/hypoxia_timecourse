# Kate, 2026-11-19: plot the shortest- vs longest-5% CDS pathway ORA from 26h_.
# Style follows make_ora_plot_3grp() in code/predictive_modeling/24_thesis_committee_2026.Rmd:
# dodged fold-enrichment bars per group, GeneRatio labels, bars faded where padj >= 0.05,
# dashed line at 1. Terms shown = the N_TOP most significant for each list (union), and both
# tails are drawn for every term shown so you can see the other tail is not enriched.
#
# Run after both tails of 26h_:
#   Rscript code/translational_control_overview/26h_shortest_cds_pathway_ora.R shortest
#   Rscript code/translational_control_overview/26h_shortest_cds_pathway_ora.R longest
#   Rscript code/translational_control_overview/26i_cds_length_tails_ora_plots.R

suppressPackageStartupMessages({ library(tidyverse); library(here) })
O <- here("output", "translational_control_overview")
PADJ <- 0.05; N_TOP <- 15

# Okabe-Ito hues that are not reserved for si3d (#E69F00) or si3e (#56B4E9).
tail_colors <- c("Shortest 5% CDS" = "#CC79A7", "Longest 5% CDS" = "#0072B2")

files <- c(shortest = "shortest_cds_5pct", longest = "longest_cds_5pct")
stopifnot("run 26h_ for both tails first" =
            all(file.exists(file.path(O, paste0(files, "_pathway_ora.csv")))))
ora <- imap_dfr(files, ~ read_csv(file.path(O, paste0(.x, "_pathway_ora.csv")), show_col_types = FALSE) %>%
                  mutate(tail = .y))
genes <- imap(files, ~ read_csv(file.path(O, paste0(.x, "_genes.csv")), show_col_types = FALSE))
cds_nt <- map_chr(genes, ~ sprintf("n = %d, CDS %s-%s nt", nrow(.x),
                                   scales::comma(round(2^min(.x$log2_cds_length))),
                                   scales::comma(round(2^max(.x$log2_cds_length)))))
lab <- c(shortest = "Shortest 5% CDS", longest = "Longest 5% CDS")
ora <- ora %>% mutate(group = factor(lab[tail], levels = names(tail_colors)))

prefixes <- c(Hallmark = "^HALLMARK_", KEGG = "^KEGG_(MEDICUS_)?", C2 = "^", GO_BP = "^GOBP_")
# Universe passed to enricher = every expressed gene (feature_matrix_tco.rds, one row per symbol).
# enricher then keeps only universe genes annotated in the collection, so the effective background
# is per collection: it is the denominator of BgRatio, and both numbers go on the plot.
fm <- readRDS(file.path(O, "feature_matrix_tco.rds"))
n_universe <- tibble(symbol = fm$symbol, l = as.numeric(unlist(fm$log2_cds_length))) %>%
  filter(!is.na(symbol), !is.na(l)) %>% distinct(symbol) %>% nrow()

plot_coll <- function(cn) {
  d <- ora %>% filter(collection == cn)
  keep <- d %>% filter(p.adjust < PADJ) %>% group_by(tail) %>% slice_min(p.adjust, n = N_TOP, with_ties = FALSE) %>%
    pull(pathway) %>% unique()
  if (!length(keep)) return(NULL)
  # A term tested in only one tail (<10 members in the other foreground's universe slice is not
  # the reason - enricher drops terms with no foreground hit) is filled in as 0 genes, FE 0.
  pd <- d %>% filter(pathway %in% keep) %>%
    complete(pathway, group, fill = list(fold_enrichment = 0, p.adjust = 1, GeneRatio = "0")) %>%
    mutate(term = str_replace(pathway, prefixes[[cn]], "") %>% str_replace_all("_", " ") %>% str_to_title(),
           term = str_trunc(term, 70), sig = p.adjust < PADJ) %>%
    group_by(pathway) %>%
    # order: shortest-enriched terms at the top, longest-enriched below, each by fold enrichment
    mutate(lead = group[which.min(p.adjust)], best = max(fold_enrichment)) %>% ungroup() %>%
    mutate(term = fct_reorder2(term, lead == "Longest 5% CDS", -best, .fun = function(x, y) -first(x) * 1e3 + first(y)))
  n_sig <- d %>% filter(p.adjust < PADJ) %>% dplyr::count(tail) %>% { setNames(.$n, .$tail) }
  p <- ggplot(pd, aes(x = fold_enrichment, y = term, fill = group, alpha = sig)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = GeneRatio), position = position_dodge(width = 0.8),
              hjust = -0.15, size = 2.6, color = "grey25", alpha = 1, show.legend = FALSE) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = tail_colors, drop = FALSE) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(x = "Fold enrichment", y = NULL, fill = NULL,
         title = sprintf("%s ORA: shortest vs longest 5%% CDS, MDA-MB-231 expressed mRNAs", sub("_", " ", cn)),
         subtitle = paste0(sprintf("Terms shown: the %d most significant for the shortest-5%% list plus the %d most significant for the longest-5%% list",
                                   N_TOP, N_TOP),
                           "\nDashed line = 1 (no enrichment); faded bars are padj >= 0.05; labels = genes in term / foreground in collection",
                           "\nShortest: ", cds_nt["shortest"], "; longest: ", cds_nt["longest"],
                           "\nUniverse = all expressed genes in MDA-MB-231 (dominant isoform, siCTRL normoxia), n = ",
                           scales::comma(n_universe), "; of these, annotated in ", sub("_", " ", cn), " (test background): n = ",
                           scales::comma(as.numeric(sub(".*/", "", d$BgRatio[1]))))) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top", plot.title.position = "plot", plot.subtitle = element_text(size = 9))
  h <- 2.2 + 0.32 * n_distinct(pd$term)
  f <- here("plots", paste0("26i_ora_cds_length_tails_", tolower(cn), ".pdf"))
  ggsave(f, p, width = 11, height = h)
  cat(sprintf("%s: %d terms -> %s\n", cn, n_distinct(pd$term), f))
  p
}

walk(names(prefixes), plot_coll)
