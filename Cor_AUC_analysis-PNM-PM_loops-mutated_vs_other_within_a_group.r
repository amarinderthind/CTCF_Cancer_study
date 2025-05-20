## Title: CTCF Loop Differential Expression and Mutation Analysis
## Author: Amarinder Thind

# Load libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DESeq2)
library(pROC)
library(infotheo)

# Load mutation matrix
table_path <- "Dropbox/Amarinder/Amarinder_main_projects/CTCF_motif_AT/2024-ctcf-dragen-additional-samples/loops-vs-mutational-matrix-EXCLUDE_GC_FAILED_3SAMPLES.csv"
Mut_matrix <- read_csv(table_path) |> as.data.frame()
rownames(Mut_matrix) <- Mut_matrix$...1
Mut_matrix <- Mut_matrix[, -1]

# Load and filter metadata
meta <- read.csv("Dropbox/Amarinder/Sample_level_data_v11_short.csv")
rownames(meta) <- meta$sample
meta <- subset(meta, Group %in% c("PRI_NonMet", "PRI_Met"))
gc_sample <- c("CSCC_0015-M1", "CSCC_0023-M1", "CSCC_0021-P1")
meta <- meta[!(meta$SampleName %in% gc_sample), ]
meta <- meta[meta$SampleName %in% colnames(Mut_matrix), ]
meta <- meta[order(meta$Group), ]
rownames(meta) <- meta$SampleName

# Filter mutation matrix to matched samples
Mut_matrix <- Mut_matrix[, meta$SampleName, drop = FALSE]

# Subset loops with >3 mutations
mut <- Mut_matrix[rowSums(Mut_matrix) > 3, ]

# Prepare data for statistical testing
mut1_pval <- t(mut) |> as.data.frame()
named_group <- setNames(meta$Group, rownames(meta))
mut1_pval$group <- named_group[rownames(mut1_pval)]

# Binarize mutation data
var_names <- setdiff(names(mut1_pval), "group")
mut1_pval[var_names] <- lapply(mut1_pval[var_names], function(x) ifelse(x >= 1, 1, 0))

# Perform chi-square or Fisher's test
compare_groups <- function(data, var) {
  tab <- table(data$group, data[, var])
  test <- if (min(tab) >= 5) chisq.test(tab) else fisher.test(tab)
  data.frame(variable = var, p_value = test$p.value)
}
results <- lapply(var_names, function(var) compare_groups(mut1_pval, var)) |> bind_rows()
results$p_adjusted <- p.adjust(results$p_value, method = "BH")

# Combine mutation and results
temp_binary <- t(mut1_pval)
combined <- cbind(temp_binary[-nrow(temp_binary), ], results)

# Compute group-wise mutation sums
primary_samples <- meta$SampleName[meta$Group == "PRI_NonMet" & !is.na(meta$SampleName)]
met_samples <- meta$SampleName[meta$Group == "PRI_Met" & !is.na(meta$SampleName)]


# Identify sample columns in the combined data
primary_cols <- which(colnames(combined) %in% primary_samples)
met_cols <- which(colnames(combined) %in% met_samples)

# Check column validity
stopifnot(length(primary_cols) > 0)
stopifnot(length(met_cols) > 0)

# Extract primary and metastasis columns and ensure binary format
temp_primary <- combined[, primary_cols, drop = FALSE]
temp_metastasis <- combined[, met_cols, drop = FALSE]

temp_primary <- data.frame(lapply(temp_primary, as.numeric))
temp_metastasis <- data.frame(lapply(temp_metastasis, as.numeric))

# Compute row sums for mutation counts (0 or 1 as binary)
combined$sumPNM <- rowSums(temp_primary, na.rm = TRUE)
combined$sumPM <- rowSums(temp_metastasis, na.rm = TRUE)


# Identify top loops
plot_loops <- combined %>%
  arrange(p_value) %>%
  filter(p_value < 0.3 & (sumPNM <= 1 | sumPM <= 1)) %>%
  mutate(loop_id = rownames(.))

# Plotting
plot_data <- plot_loops %>%
  pivot_longer(cols = c(sumPNM, sumPM), names_to = "Category", values_to = "Count") %>%
  mutate(Count_mod = ifelse(Category == "sumPNM", -Count, Count))

max_count <- ceiling(max(abs(plot_data$Count_mod)))

ggplot(plot_data, aes(x = reorder(loop_id, -abs(Count_mod)), y = Count_mod, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(x = "Loop ID", y = "Count", title = "PNM vs PM Signal in DE Loops") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  coord_flip() +
  scale_y_continuous(breaks = seq(-max_count, max_count, by = 1), labels = abs)

# Prepare data for expression analysis
loops <- combined
loops2 <- loops[, 1:(ncol(loops) - 5)]

# Align metadata to RNA-seq samples
common <- intersect(meta$wgs_samples, meta$totalRNAseq)
meta <- subset(meta, SampleName %in% common)

# Load loop-specific gene lists and raw counts
# Load static data
all_genes <- read.csv("Dropbox/Amarinder/Amarinder_main_projects/CTCF_motif_AT/2024-ctcf/all_gene_biomart.csv")
gene_list_loop <- read.csv("Dropbox/Amarinder/Amarinder_main_projects/CTCF_motif_AT/2024-ctcf-dragen-additional-samples/Step1.b-extraction-of-loops-gene-list/Loops_with_gene_list_combine-inside-and-1000bp-outside-.csv")
rawcount_ori <- read_csv("Dropbox/Amarinder/Amarinder_main_projects/CTCF_motif_AT/2024-ctcf-dragen-additional-samples/feature_count-v3-80samples/ftype-exon-fattr-gene_id/countmatrix_genes_Modified_samp_names-unique-ensemble-id-v2.csv") |> as.data.frame()


rownames(rawcount_ori) <- rawcount_ori[[1]]
rawcount_ori <- rawcount_ori[, -1]
rawcount_ori[is.na(rawcount_ori)] <- 0

loop_ids <- c("loop_70","loop_95","loop_657","loop_53","loop_45","loop_1083","loop_284","loop_390","loop_937","loop_818","loop_115")
#loop_ids <- c("loop_70","loop_45")

final_results <- list()

meta2 <- meta
loops3 <- loops2

for (loop_id in loop_ids) {
  cat("\n\nProcessing", loop_id, "...\n")
  
  # Filter mutation and metadata
  
  
  
  # if(loop_id %in% c("loop_70","loop_95","loop_657","loop_53")) {
  #   meta <- subset(meta2, Group %in% c("PRI_Met"))  # special case filtering
  # } else {
  #   meta <- subset(meta2, Group %in% c("PRI_NonMet"))  # default case filtering (same here)
  # }
  
  loops2 <- loops3[, meta$SampleName, drop = FALSE]
  
  # Gene list for this loop
  gene_list_loop$hgnc_symbol <- gsub("\\s+", "", gene_list_loop$hgnc_symbol)
  gene_list_loop$ensembl_gene_id <- gsub("\\s+", "", gene_list_loop$ensembl_gene_id)
  intersect_genes <- gene_list_loop[gene_list_loop$identifier %in% loop_id,]
  cleaned_genes <- intersect_genes$ensembl_gene_id[intersect_genes$ensembl_gene_id != ""]
  all_Loopgenes <- unique(trimws(unlist(strsplit(cleaned_genes, ",\\s*"))))
  
  # Filter rawcount
  keep <- rowSums(rawcount_ori > 0) >= round(ncol(rawcount_ori) * 0.15)
  rawcount <- rawcount_ori[keep,]
  
  # Mutation status vector
  mut_status <- as.numeric(loops2[loop_id, ])
  names(mut_status) <- colnames(loops2)
  mut_binary <- ifelse(mut_status >= 1, 1, 0)
  
  meta$loop_status <- NA
  meta[names(mut_binary), "loop_status"] <- mut_binary
  meta$loop_status <- ifelse(meta$loop_status == 1, "mut", ifelse(meta$loop_status == 0, "nonmut", NA))
  
  # Subset and align
  anno <- subset(meta, loop_status %in% c("nonmut", "mut"))
  anno <- anno[anno$SampleName != "CSCC_0023-M1", ]
  rawcount <- rawcount[, names(rawcount) %in% anno$SampleName]
  anno <- anno[match(colnames(rawcount), anno$SampleName),]
  rownames(anno) <- anno$SampleName
  
  # DESeq2
  dds <- DESeqDataSetFromMatrix(countData = round(rawcount), colData = anno, design = ~totalRNA_Batch + loop_status)
  dds <- DESeq(dds)
  vst <- varianceStabilizingTransformation(dds, blind = TRUE)
  normalized_counts <- as.data.frame(assay(vst))
  
  # Extract gene expression for loop genes
  rc <- normalized_counts[rownames(normalized_counts) %in% all_Loopgenes, ]
  rc$ensembl_gene_id <- rownames(rc)
  rc <- rc %>%
    left_join(all_genes %>% select(hgnc_symbol, ensembl_gene_id), by = "ensembl_gene_id") %>%
    mutate(hgnc_symbol = ifelse(hgnc_symbol == "" | is.na(hgnc_symbol), ensembl_gene_id, hgnc_symbol))
  rownames(rc) <- rc$hgnc_symbol
  rc <- rc[, !(names(rc) %in% c("ensembl_gene_id", "hgnc_symbol"))]
  rc <- t(rc)
  
  # Prepare data for analysis
  glm_data <- data.frame(label = factor(anno$loop_status, levels = c("nonmut", "mut")), rc)
  glm_data$label_numeric <- ifelse(glm_data$label == "mut", 1, 0)
  
  # Correlation
  gene_correlations <- apply(rc, 2, function(gene_expr) {
    cor(gene_expr, glm_data$label_numeric, method = "spearman")
  })
  
  # AUC
  auc_scores <- apply(rc, 2, function(expr) {
    roc_obj <- roc(glm_data$label, expr, quiet = TRUE)
    auc(roc_obj)
  })
  
  # 
  
  # Combine and store
  loop_result <- data.frame(
    Gene = names(gene_correlations),
    SpearmanR = gene_correlations,
    AUC = auc_scores[names(gene_correlations)]
  )
  final_results[[loop_id]] <- loop_result
}

# Combine all loop results into one table with loop_id
summary_table <- do.call(rbind, lapply(names(final_results), function(lid) {
  cbind(LoopID = lid, final_results[[lid]])
}))

# View combined summary
head(summary_table)

#write.csv(summary_table,"step-6-PriNM-PriMet-progression/correlatino_AUC_for_mut_vs_other_slected_loops_PNM_PM.csv")




