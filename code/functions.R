# Helper Functions for RiboPipe Shiny App
# Based on Mukherjee Lab reference analysis workflow
# Reference: https://github.com/mukherjeelab/2024_translationInhib_HNSCC

library(DESeq2)
library(dplyr)
library(tibble)

#' Create Normalized Count Matrix with Spike-In Support
#'
#' Generates DESeq2-normalized counts for ribosome profiling data with optional
#' Drosophila spike-in normalization. This function follows the reference workflow
#' from the Mukherjee Lab 2024 HNSCC study.
#'
#' @param allCounts Count matrix (genes x samples) with both human/mouse and spike-in genes
#' @param metadata Data frame with columns: sample_name, library (rna/ribo), comparison_variable
#' @param control_condition String specifying the reference level for comparison_variable
#' @param spikein Optional spike-in count matrix (spike-in genes x samples)
#' @param keepSpikes Optional character vector of spike-in gene IDs to use for normalization
#'
#' @return Matrix of normalized counts (same dimensions as allCounts)
#'
#' @details
#' The function performs the following steps:
#' 1. Combines allCounts with spike-ins if provided
#' 2. Creates DESeq2 dataset with design: ~ comparison_variable + library + comparison_variable:library
#' 3. Estimates size factors using spike-ins (via controlGenes parameter) if provided
#' 4. Returns normalized counts for all genes (including spike-ins)
#'
#' Reference: 01-QC.Rmd lines 46-111
#'
createNormMatrix <- function(allCounts, metadata, control_condition,
                              spikein = NULL, keepSpikes = NULL) {

  # Round counts to integers
  allCounts <- round(as.matrix(allCounts), 0)

  # Check if we have valid spike-in data
  has_spike <- !is.null(spikein) &&
    ((is.matrix(spikein) || is.data.frame(spikein)) && nrow(spikein) > 0)

  if (has_spike) {
    message("Using spike-in normalization with ", length(keepSpikes), " control genes")

    spikein <- round(as.matrix(spikein), 0)

    # Ensure spike-in samples match allCounts samples
    if (!identical(colnames(allCounts), colnames(spikein))) {
      if (!all(colnames(allCounts) %in% colnames(spikein))) {
        stop("spikein is missing some sample columns found in allCounts.")
      }
      spikein <- spikein[, colnames(allCounts), drop = FALSE]
    }

    # Combine counts with spike-ins
    counts_all <- rbind(allCounts, spikein)

    # Identify spike-in indices for controlGenes parameter
    if (!is.null(keepSpikes) && length(keepSpikes) > 0) {
      Spikeindices <- which(rownames(counts_all) %in% keepSpikes)
      message("Found ", length(Spikeindices), " spike-in genes in combined matrix")
    } else {
      Spikeindices <- NULL
    }
  } else {
    message("Using standard DESeq2 normalization (no spike-ins)")
    counts_all <- allCounts
    Spikeindices <- NULL
  }

  # Prepare metadata for DESeq2
  # Ensure metadata rows match count columns
  metadata <- metadata %>%
    filter(sample_name %in% colnames(counts_all)) %>%
    arrange(match(sample_name, colnames(counts_all)))

  if (!identical(metadata$sample_name, colnames(counts_all))) {
    stop("Metadata sample names do not match count matrix column names")
  }

  # Create sample info with proper factor levels
  sample_info <- metadata %>%
    dplyr::select(library, comparison_variable) %>%
    mutate(
      SeqType = factor(library, levels = c("rna", "ribo")),
      Condition = factor(comparison_variable)
    )

  # Relevel Condition to use control as reference
  if (!is.null(control_condition) && length(control_condition) > 0 &&
      control_condition %in% levels(sample_info$Condition)) {
    sample_info$Condition <- relevel(sample_info$Condition, ref = control_condition)
  } else {
    warning("control_condition not provided or not found in metadata. Using first level as reference.")
  }

  # Create DESeq2 dataset
  # Design includes interaction term for TE analysis
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_all,
    colData = sample_info,
    design = ~ Condition + SeqType + Condition:SeqType
  )

  # Estimate size factors - matches Mukherjee Lab reference workflow
  # No tryCatch/poscounts fallback to avoid NA size factors
  if (is.null(Spikeindices) || length(Spikeindices) == 0) {
    dds <- DESeq2::estimateSizeFactors(dds)
    message("Size factors estimated using standard DESeq2 method")
  } else {
    dds <- DESeq2::estimateSizeFactors(dds, controlGenes = Spikeindices)
    message("Size factors estimated using ", length(Spikeindices), " spike-in control genes")
  }

  # Extract normalized counts
  norm <- DESeq2::counts(dds, normalized = TRUE)

  message("DESeq2 normalization complete. Matrix dimensions: ",
          nrow(norm), " genes x ", ncol(norm), " samples")

  return(norm)
}


#' Differential Expression Analysis with Spike-In Support
#'
#' Performs DESeq2 differential expression analysis for ribosome profiling data
#' with support for TE (translation efficiency), RNA, and RIBO modes.
#' Follows the reference deFunction workflow from Mukherjee Lab 2024 HNSCC study.
#'
#' @param countsMatrix Raw count matrix (genes x samples)
#' @param allCounts_Norm Pre-normalized count matrix from createNormMatrix()
#' @param metadata Data frame with sample metadata
#' @param control_condition String specifying control condition name
#' @param test_condition String specifying test condition name
#' @param spikeInMatrix Optional spike-in count matrix
#' @param keepSpikes Optional character vector of spike-in gene IDs
#' @param mode Analysis mode: "TE" (default), "RNA", or "RIBO"
#'
#' @return Data frame with DESeq2 results merged with normalized counts
#'
#' @details
#' The function performs differential expression analysis with the following workflow:
#' 1. Filters samples to selected conditions
#' 2. Combines counts with spike-ins (if provided)
#' 3. Creates DESeq2 dataset with appropriate design:
#'    - TE mode: ~ Condition + SeqType + Condition:SeqType (interaction term)
#'    - RNA/RIBO mode: ~ Condition (simple comparison)
#' 4. Estimates size factors using spike-ins if provided
#' 5. **REMOVES spike-ins from dataset** before running DESeq2 (critical step!)
#' 6. Runs DESeq2 differential analysis
#' 7. Merges results with normalized counts
#'
#' Reference: 02-APL_analysis.Rmd lines 44-164
#'
deFunction <- function(countsMatrix, allCounts_Norm, metadata,
                       control_condition, test_condition,
                       spikeInMatrix = NULL, keepSpikes = NULL,
                       mode = c("TE", "RNA", "RIBO")) {

  mode <- match.arg(mode)
  message("\n=== Running deFunction in ", mode, " mode ===")
  message("Control: ", control_condition, " vs Test: ", test_condition)

  # Round counts to integers
  allCounts <- round(as.matrix(countsMatrix), 0)

  # Remove any rows or columns that are all NA
  allCounts <- allCounts[!apply(is.na(allCounts), 1, all), , drop = FALSE]
  allCounts <- allCounts[, !apply(is.na(allCounts), 2, all), drop = FALSE]

  # Check if we have valid spike-in data
  is_empty_spike <- is.null(spikeInMatrix) ||
    (is.matrix(spikeInMatrix) && nrow(spikeInMatrix) == 0) ||
    (is.data.frame(spikeInMatrix) && nrow(spikeInMatrix) == 0)

  if (!is_empty_spike) {
    message("Processing spike-in data for normalization")
    spikeInMatrix <- round(as.matrix(spikeInMatrix), 0)

    # Ensure spike-in samples match allCounts samples
    if (!identical(colnames(allCounts), colnames(spikeInMatrix))) {
      spikeInMatrix <- spikeInMatrix[, colnames(allCounts), drop = FALSE]
    }

    # Combine counts with spike-ins
    countsMatrix <- rbind(allCounts, spikeInMatrix)

    # Get spike-in indices
    if (!is.null(keepSpikes) && length(keepSpikes) > 0) {
      Spikeindices <- which(rownames(countsMatrix) %in% keepSpikes)
      message("Found ", length(Spikeindices), " spike-in genes for size factor estimation")
    } else {
      Spikeindices <- NULL
    }
  } else {
    countsMatrix <- allCounts
    Spikeindices <- NULL
  }

  # Filter metadata to selected conditions - REMOVE NAs FIRST
  metadata_filtered <- metadata %>%
    filter(!is.na(comparison_variable)) %>%
    filter(!is.na(library)) %>%
    filter(comparison_variable %in% c(control_condition, test_condition))

  # Validate we have samples after filtering
  if (nrow(metadata_filtered) == 0) {
    stop("No samples found for conditions: ", control_condition, ", ", test_condition)
  }

  # Filter counts to selected conditions AND library type based on mode
  if (mode == "RNA") {
    samples_use <- metadata_filtered %>%
      filter(library == "rna") %>%
      pull(sample_name)
    message("Filtering to RNA samples only (n = ", length(samples_use), ")")
  } else if (mode == "RIBO") {
    samples_use <- metadata_filtered %>%
      filter(library == "ribo") %>%
      pull(sample_name)
    message("Filtering to RIBO samples only (n = ", length(samples_use), ")")
  } else {
    # TE mode: use both RNA and RIBO
    samples_use <- metadata_filtered %>% pull(sample_name)
    message("Using both RNA and RIBO samples (n = ", length(samples_use), ")")
  }

  # Validate we have samples after filtering by library type
  if (length(samples_use) == 0) {
    stop("No ", mode, " samples found for the selected conditions")
  }

  # Filter counts and metadata
  counts_use <- countsMatrix[, samples_use, drop = FALSE]
  metadata_use <- metadata_filtered %>%
    filter(sample_name %in% samples_use) %>%
    arrange(match(sample_name, colnames(counts_use)))

  # Prepare sample info for DESeq2
  sample_info <- metadata_use %>%
    mutate(
      Condition = factor(comparison_variable, levels = c(control_condition, test_condition)),
      SeqType = factor(library, levels = c("rna", "ribo"))
    ) %>%
    dplyr::select(Condition, SeqType)

  # Create DESeq2 dataset with appropriate design
  if (mode == "TE") {
    message("Using interaction design: ~ Condition + SeqType + Condition:SeqType")
    ddsMat <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts_use,
      colData = sample_info,
      design = ~ Condition + SeqType + Condition:SeqType
    )
  } else {
    message("Using simple design: ~ Condition")
    ddsMat <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts_use,
      colData = sample_info %>% dplyr::select(Condition),
      design = ~ Condition
    )
  }

  # Estimate size factors
  if (is.null(Spikeindices) || length(Spikeindices) == 0) {
    ddsMat <- DESeq2::estimateSizeFactors(ddsMat)
    message("Size factors estimated using standard DESeq2 method")
  } else {
    ddsMat <- DESeq2::estimateSizeFactors(ddsMat, controlGenes = Spikeindices)
    message("Size factors estimated using ", length(Spikeindices), " spike-in control genes")
  }

  # CRITICAL: Remove spike-ins after size factor estimation
  if (!is.null(Spikeindices) && length(Spikeindices) > 0) {
    message("Removing spike-ins from dataset before differential analysis")
    dds_noSpikes <- ddsMat[-Spikeindices, ]
    # Preserve size factors from spike-in normalization
    sizeFactors(dds_noSpikes) <- sizeFactors(ddsMat)
    ddsMat <- dds_noSpikes
    message("Dataset now contains ", nrow(ddsMat), " genes (spike-ins removed)")
  }

  # Run DESeq2 differential analysis
  message("Running DESeq2 differential analysis...")
  ddsMat <- DESeq2::DESeq(ddsMat)

  # Get available result names
  all_results <- resultsNames(ddsMat)
  message("Available result names from DESeq2:")
  print(all_results)

  # Select appropriate coefficient based on mode
  if (mode == "TE") {
    # TE mode: use interaction term
    # Format: Condition{test}_vs_{control}.SeqTyperibo
    # Or: Condition{test}.SeqTyperibo (depending on DESeq2 version)

    # Try both possible naming conventions
    te_coef_name1 <- paste0("Condition", test_condition, ".SeqTyperibo")
    te_coef_name2 <- paste0("Condition", test_condition, "_vs_", control_condition, ".SeqTyperibo")

    if (te_coef_name1 %in% all_results) {
      coef_name <- te_coef_name1
    } else if (te_coef_name2 %in% all_results) {
      coef_name <- te_coef_name2
    } else {
      # Fallback: look for any interaction term with SeqTyperibo
      interaction_terms <- grep("SeqTyperibo", all_results, value = TRUE)
      if (length(interaction_terms) > 0) {
        coef_name <- interaction_terms[1]
        warning("Using fallback interaction term: ", coef_name)
      } else {
        stop("TE interaction term not found. Available: ", paste(all_results, collapse = ", "))
      }
    }
  } else {
    # RNA or RIBO mode: use main effect
    # Format: Condition_test_vs_control
    coef_name <- paste0("Condition_", test_condition, "_vs_", control_condition)

    if (!coef_name %in% all_results) {
      # Fallback: use index 2 (typically the main comparison)
      if (length(all_results) >= 2) {
        coef_name <- all_results[2]
        warning("Using fallback coefficient: ", coef_name)
      } else {
        stop("Expected coefficient not found. Available: ", paste(all_results, collapse = ", "))
      }
    }
  }

  message("Using coefficient: ", coef_name)

  # Extract results
  res <- DESeq2::results(ddsMat, name = coef_name, independentFiltering = FALSE)

  # Apply LFC shrinkage
  res <- DESeq2::lfcShrink(dds = ddsMat, coef = coef_name, res = res, type="apeglm")

  # Convert to data frame
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id")

  # Merge with normalized counts
  # Filter normalized counts to relevant samples
  norm_counts_filtered <- allCounts_Norm[, samples_use, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("gene_id")

  # Merge DESeq2 results with normalized counts
  res_final <- norm_counts_filtered %>%
    right_join(res_df, by = "gene_id")

  message("Analysis complete. Results contain ", nrow(res_final), " genes")
  message("Significant genes (padj < 0.05): ", sum(res_final$padj < 0.05, na.rm = TRUE))

  return(res_final)
}


#' Categorize Translational Regulation Changes
#'
#' Categorizes genes based on RNA-seq, Ribo-seq, and TE changes following
#' Harnett et al. 2022 framework for translational regulation.
#'
#' @param res_rna DESeq2 results from RNA-seq analysis (output from deFunction with mode="RNA")
#' @param res_ribo DESeq2 results from Ribo-seq analysis (output from deFunction with mode="RIBO")
#' @param res_te DESeq2 results from TE analysis (output from deFunction with mode="TE")
#' @param padj_cutoff Adjusted p-value threshold (default: 0.05)
#' @param lfc_cutoff Log2 fold-change threshold (default: 0.58496250072, ~1.5-fold)
#'
#' @return List with three elements:
#'   - plot: ggplot2 scatter plot of RNA vs Ribo log2FC colored by category
#'   - categorized_table: Data frame with gene_id, LFCs, and category assignment
#'   - summary_table: Count table of genes per category
#'
#' @details
#' Categories:
#' - Forwarded: RNA and Ribo both significant, TE not significant (transcriptional regulation)
#' - Exclusive: RNA not significant, Ribo and TE significant (pure translational regulation)
#' - Buffered: RNA, Ribo, and TE all significant OR RNA and TE significant but not Ribo (compensatory regulation)
#' - TE_only: TE significant but RNA and Ribo not individually significant
#' - Not_sig: No significant changes
#'
#' Reference: Harnett et al. 2022
#'
categorize_translation_changes <- function(res_rna, res_ribo, res_te,
                                          padj_cutoff = 0.05,
                                          lfc_cutoff = 0.58496250072) {

  # Ensure gene_id alignment across datasets
  common_genes <- Reduce(intersect, list(res_rna$gene_id, res_ribo$gene_id, res_te$gene_id))
  res_rna  <- res_rna  %>% filter(gene_id %in% common_genes) %>% arrange(gene_id)
  res_ribo <- res_ribo %>% filter(gene_id %in% common_genes) %>% arrange(gene_id)
  res_te   <- res_te   %>% filter(gene_id %in% common_genes) %>% arrange(gene_id)

  stopifnot(identical(res_rna$gene_id, res_ribo$gene_id),
            identical(res_rna$gene_id, res_te$gene_id))

  # Identify significant genes per dataset
  sig_rna  <- res_rna  %>% filter(padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)  %>% pull(gene_id)
  sig_ribo <- res_ribo %>% filter(padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)  %>% pull(gene_id)
  sig_te   <- res_te   %>% filter(padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)  %>% pull(gene_id)

  # Combine into a single data frame
  merged <- data.frame(
    gene_id = res_rna$gene_id,
    rna_lfc = res_rna$log2FoldChange,
    ribo_lfc = res_ribo$log2FoldChange,
    te_lfc = res_te$log2FoldChange
  )

  # Categorize genes
  merged <- merged %>%
    mutate(
      category = case_when(
        gene_id %in% sig_rna  & gene_id %in% sig_ribo & !(gene_id %in% sig_te) ~ "Forwarded",
        !(gene_id %in% sig_rna) & gene_id %in% sig_ribo & gene_id %in% sig_te  ~ "Exclusive",
        gene_id %in% sig_rna & gene_id %in% sig_ribo & gene_id %in% sig_te     ~ "Buffered",
        gene_id %in% sig_rna & !(gene_id %in% sig_ribo) & gene_id %in% sig_te  ~ "Buffered",
        gene_id %in% sig_te & !(gene_id %in% sig_ribo | gene_id %in% sig_rna)  ~ "TE_only",
        TRUE                                                                   ~ "Not_sig"
      )
    ) %>%
    mutate(category = factor(category,
                             levels = c("Not_sig","Forwarded","Exclusive","Buffered","TE_only")))

  # Summarize counts per category
  cat_table <- table(merged$category) %>% as.data.frame()
  colnames(cat_table) <- c("Category","Count")

  # Create scatter plot (RNA vs Ribo)
  myCols <- c("Not_sig" = "grey80",
              "Forwarded" = "blue",
              "Exclusive" = "red",
              "Buffered" = "green",
              "TE_only" = "purple")

  p <- ggplot(merged %>% filter(category != "Not_sig"),
              aes(x = rna_lfc, y = ribo_lfc, color = category)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = myCols) +
    xlim(-5, 5) + ylim(-5, 5) +
    xlab("RNA-seq log2FC") +
    ylab("Ribo-seq log2FC") +
    theme_minimal(base_size = 13) +
    theme(legend.title = element_blank())

  # Return list of outputs
  list(
    plot = p,
    categorized_table = merged,
    summary_table = cat_table
  )
}
