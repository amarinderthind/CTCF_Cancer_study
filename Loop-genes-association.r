## Title: "Differential Expression Analysis of RNA-seq Data Associated with Genomic Loops"

## Load Required Packages
library(readr)
library(dplyr)
library(DESeq2)
library(edgeR)
library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(foreach)
library(doParallel)

# Load loop annotations
loops <- read_csv("/path/to/loops-vs-mutational-matrix.csv")

dim(loops)

# Load sample annotation
anno2 <- read_csv("Sample_level_data_v11_short.csv")
 
anno2 <- anno2 %>%
  filter(!is.na(Group) & Group %in% c("Met", "PRI_Met", "PRI_NonMet"))

# Find common samples between RNA-seq and WGS
common <- intersect(anno2$wgs_samples, anno2$totalRNAseq)
common1 <- common[!grepl("MCL", common, ignore.case = TRUE) & !is.na(common)]
print(paste("Number of shared samples (excluding cell lines):", length(common1)))

#### Prepare Raw Count Matrix

# Load raw counts
rawcount_ori <- read_csv("/path/to/countmatrix_genes.csv")
rawcount_ori <- as.data.frame(rawcount_ori)

# Set rownames
rownames(rawcount_ori) <- rawcount_ori[[1]]
rawcount_ori <- rawcount_ori[, -1]

# Replace NA with 0
rawcount_ori[is.na(rawcount_ori)] <- 0

# Filter genes expressed in ≥15% of samples
keep <- rowSums(rawcount_ori > 0) >= round(ncol(rawcount_ori) * 0.15)
rawcount_ori <- rawcount_ori[keep,]

### Load Gene Lists Associated with Loops
gene_list_loop <- read.csv("Step1.b-extraction-of-loops-gene-list/Loops_with_gene_list_combine-inside-and-1000bp-outside-.csv")


gene_list_loop$hgnc_symbol <- gsub("\\s+", "", gene_list_loop$hgnc_symbol)
gene_list_loop$ensembl_gene_id <- gsub("\\s+", "", gene_list_loop$ensembl_gene_id)

all_genes <- read_csv("all_gene_biomart.csv")
 
### Preprocess Loops for DE Analysis

# Keep only columns corresponding to shared RNA-seq samples
loops <- loops[, colnames(loops) %in% c(anno2$totalRNAseq, "...1")]

# Filter loops with mutations present in at least 4 samples
loops2 <- loops[rowSums(loops > 0) > 3, ]

#### Differential Expression Analysis Using DESeq2

# Parallelization Setup
numCores <- detectCores() - 2
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Chunk Loops for Parallel Processing
split_into_chunks <- function(data, chunk_size) {
  split(data, ceiling(seq_along(data) / chunk_size))
}
loop_chunks <- split_into_chunks(loops2$...1, chunk_size = numCores)

# Initialize Result Storage
Final_summary2 <- data.frame()
meta_list <- list()

# Loop over chunks
for (loop_chk in loop_chunks) {
  print(paste("Starting analysis for loops:", paste(loop_chk, collapse=", ")))

  Final_summary <- foreach(value = loop_chk, .combine = 'rbind', .packages = c('dplyr')) %dopar% {
    
    # Subset loop-specific mutation matrix
    loopSelec <- loops2[loops2$...1 == value, ]
    loopSelec <- loopSelec[, colnames(loopSelec) %in% anno2$totalRNAseq]
    
    # Prepare metadata
    meta <- t(loopSelec) %>% as.data.frame()
    meta$V1 <- ifelse(meta$V1 > 0, "mut", "no-mut")
    colnames(meta) <- "group"
    meta$SampleName <- rownames(meta)
    meta <- merge(meta, anno2[, c("SampleName", "totalRNA_Batch")], by = "SampleName", all.x = TRUE)
    
    # Remove problematic sample ## Low overall count and GC bias
    meta <- meta[meta$SampleName != "CSCC_0023-M1", ]
    
    # Match and subset counts
    rawcount <- rawcount_ori[, colnames(rawcount_ori) %in% meta$SampleName]
    meta <- meta[match(colnames(rawcount), meta$SampleName), ]
    rownames(meta) <- meta$SampleName

    # DESeq2 Analysis
    dds <- DESeqDataSetFromMatrix(countData = round(rawcount), colData = meta, design = ~ totalRNA_Batch + group)
    dds <- DESeq(dds)
    contrast <- c("group", "no-mut", "mut")
    res <- results(dds, contrast = contrast)
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)

    # Select genes related to the loop
    intersect_genes <- gene_list_loop[gene_list_loop$identifier == value, ]
    loop_genes <- unique(trimws(unlist(strsplit(paste(intersect_genes$ensembl_gene_id, collapse=","), ","))))
    
    # Subset DE results
    filtered_df <- res_df[res_df$gene %in% loop_genes, ]


# Identify significant genes
    significant_genes <- filtered_df %>%
      filter(abs(log2FoldChange) > 1 & padj  < 0.05)

# Summary Table

   data.frame(
      Loop_Name = value,
      Number_of_Significant_Genes = nrow(significant_genes),
      Gene_Names = if (nrow(significant_genes) > 0) paste(significant_genes$gene, collapse = ", ") else "None",
      Number_of_tested_Genes = nrow(filtered_df),
      Gene_RNAseq_tested = if (nrow(filtered_df) > 0) paste(filtered_df$gene, collapse = ", ") else "None",
      Number_of_Genes_inLoop = length(loop_genes),
      Gene_loop = paste(loop_genes, collapse = ", "),
      Mut_num = sum(meta$group == "mut"),
      No_mut = sum(meta$group == "no-mut"),
      stringsAsFactors = FALSE
    )
  }

  Final_summary2 <- rbind(Final_summary2, Final_summary)
}

# Stop parallel processing
stopCluster(cl)
### Post-processing and Summary Table Creation

# Format output
colnames(Final_summary2) <- c("Loop", "No_Sig_Genes", "Sig_Genes", "No_Genes_RNAseq",
                              "RNAseqGenes", "No_Loops_genes", "Loop_gene_names",
                              "SampleSize_Mut", "SampleSize_noMut")

# Map Ensembl IDs to HGNC Symbols
mapping <- setNames(all_genes$hgnc_symbol, all_genes$ensembl_gene_id)

replace_with_hgnc <- function(ids) {
  id_list <- unlist(strsplit(ids, ","))
  hgnc_list <- sapply(id_list, function(x) mapping[x])
  paste(hgnc_list, collapse = ",")
}

Final_summary2$Sig_Genes_hgnc <- sapply(Final_summary2$Sig_Genes, replace_with_hgnc)
Final_summary2$RNAseqGenes_hgnc <- sapply(Final_summary2$RNAseqGenes, replace_with_hgnc)
Final_summary2$loopGenes_hgnc <- sapply(Final_summary2$Loop_gene_names, replace_with_hgnc)

# Identify loops with significant DEGs
loops_sigDEG <- Final_summary2$Loop[Final_summary2$No_Sig_Genes > 0]
loops_with_sigDEG <- loops2[loops2$...1 %in% loops_sigDEG, ]

# Save final outputs (if needed)
# write.csv(Final_summary2, "Final_summary3_association_genesDE_vs_mut.csv")
# write.csv(loops_with_sigDEG, "Loops_with_Significant_DEGs.csv")
