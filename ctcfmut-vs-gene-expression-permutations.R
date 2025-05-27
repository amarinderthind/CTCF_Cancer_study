## Permutation test to test the significance of number of loops which are associated with gene-expression
## 1000 permuatition were performed (1000*284 DGE analysis)

# Load necessary libraries
library(readr)
library(readxl)

library(foreach)
library(doParallel)

setwd("/ctcf-simulations/")
# Set file paths
loops_path <- "loops-vs-mutational-matrix-EXCLUDE_GC_FAILED_3SAMPLES.csv"
#anno_path <- "Sample_level_data_v11_short.csv"
gene_list_loop <- "Loops_with_gene_list_combine-inside-and-1000bp-outside-.csv"
rawcount_path <- "countmatrix_genes_Modified_samp_names-unique-ensemble-id-v2.csv"

# Read data
loops <- read_csv(loops_path)

anno2 <- read_csv("Sample_level_data_v11_short.csv")
common <- intersect(anno2$wgs_samples, anno2$totalRNAseq)
#cat("Number of shared samples b/w RNAseq and WGS (overall including cell-lines):", length(common) - 1, "\n")
# Exclude entries containing "MCL"
common1 <- common[!grepl("MCL", common, ignore.case = TRUE) & !is.na(common)]

cat("Number of shared samples b/w RNAseq and WGS (overall including cell-lines):", length(common1) - 1, "\n")

loops <- loops[,colnames(loops) %in% c(anno2$totalRNAseq,"...1")]  #intersecrt with samples shared with RNaseq
# Remove rows without any "yes" values
loops2 <- loops[rowSums(loops > 0) > 3, ]  ## minimum 2 Samples, per Group required for Deseq2


# Prepare gene lists
gene_list_loop <- read.csv(gene_list_loop)
# Remove trailing commas from the 'genes' column
gene_list_loop$hgnc_symbol <- gsub("\\s+", "", gene_list_loop$hgnc_symbol)
gene_list_loop$ensembl_gene_id <- gsub("\\s+", "", gene_list_loop$ensembl_gene_id)



all_genes <- read.csv("all_gene_biomart.csv")
rawcount_ori <- read_csv(rawcount_path)

rawcount_ori <- as.data.frame(rawcount_ori)
names3 <- rawcount_ori$...1
rawcount_ori <- rawcount_ori[, -1]
rownames(rawcount_ori) <- names3
## Replace NAs by zero and changing the input to required format
rawcount_ori[is.na(rawcount_ori)] <- 0

# Discard genes expressed in less than 15% of all samples
keep <- rowSums(rawcount_ori > 0) >= round(ncol(rawcount_ori) * 0.15)
rawcount_ori <- rawcount_ori[keep,]

# Initialize summary table
summary_table <- data.frame(
  Loop_Name = character(),
  Number_of_Significant_Genes = integer(),
  stringsAsFactors = FALSE
)


# Parallel processing setup
numCores <- detectCores() 

cat("number of core", numCores,"\n")

cl <- makeCluster(numCores)
registerDoParallel(cl)

# Chunking and processing loops
chunk_size <- 48

#loops2 <- loops2[1:5,] ## for quick test

loop_chunks <- split(loops2$...1, ceiling(seq_along(loops2$...1) / chunk_size))



set.seed(511)
Final_summary3 <- data.frame()

for (ite in 1:1000) {

      	Final_summary2 <- data.frame()

cat("Starting Iteration:", ite, "\n")
  
  


for (loop_chk in loop_chunks) {

    #cat("Starting Analysis for Loop:", loop_chk, "\n")
    
    Final_summary <- foreach(value = loop_chk, .combine = 'rbind', .packages = c('dplyr')) %dopar% {
#value="loop_355"
            
      library(DESeq2)
      library(edgeR)
      library(tidyverse)
      library(ggrepel)
      library(RColorBrewer)
      
      loopSelec <- loops2[loops2$...1 == value,]
      loopSelec <- loopSelec[, colnames(loopSelec) %in% anno2$totalRNAseq]
      
      meta <- as.data.frame(t(loopSelec))
      meta$V1[meta$V1 > 0] <- "mut"
      meta$V1[meta$V1 == 0] <- "no-mut"
      colnames(meta) <- "group"
      
      group_sizes <- table(meta$group)
      shuffled_groups <- unlist(mapply(function(group, size) rep(group, size), names(group_sizes), group_sizes))
      shuffled_groups
      
      
      shuffled_groups <- sample(shuffled_groups)
      meta$Shuffled <- shuffled_groups
      
      meta$SampleName <- row.names(meta)
      meta <- merge(meta, anno2[, c("SampleName", "totalRNA_Batch")], by = "SampleName", all.x = TRUE)
      
      rawcount <- rawcount_ori[, names(rawcount_ori) %in% meta$SampleName]
      

      meta <- meta[match(colnames(rawcount), meta$SampleName),]
      rownames(meta) <- meta$SampleName
      
      dds <- DESeqDataSetFromMatrix(countData = round(rawcount), colData = meta, design = ~ Shuffled)
      dds <- DESeq(dds)
      
      contrast <- c("Shuffled", "no-mut", "mut")
      res <- results(dds)
      res_df <- data.frame(res)
      res_df$gene <- row.names(res_df)
      
      
      intersect_genes <- gene_list_loop[gene_list_loop$identifier %in% value,]
      
      
      ## extract clean gene list
      DE_loops_related_genes <- intersect_genes$ensembl_gene_id 
      # Remove empty entries
      cleaned_genes <- DE_loops_related_genes[DE_loops_related_genes != ""]
      
      # Split entries with multiple genes into individual genes and trim whitespaces
      split_genes <- unlist(strsplit(cleaned_genes, ",\\s*"))
      #split_genes
      # Remove any extra spaces around gene names
      split_genes <- trimws(split_genes)
      # Remove duplicate entries
      all_Loopgenes <- unique(split_genes)
      
      
      
      filtered_df <- res_df[res_df$gene %in% all_Loopgenes & !is.na(res_df$padj), ]

  # Filter based on log2FoldChange and padj
significant_genes <- filtered_df %>%
filter((log2FoldChange < -1 | log2FoldChange > 1) & padj < 0.05)

      
      summary_table <- data.frame(
        Loop_Name = value,
        Number_of_Significant_Genes = nrow(significant_genes),
        stringsAsFactors = FALSE
      )
     
# Clean up memory
      rm(loopSelec, meta, dds, res, res_df, intersect_genes_inside, intersect_genes_outside, DE_loops_related_genes, DE_loops_related_genes_outside, all_Loopgenes, filtered_df, significant_genes)
      gc()

      return(summary_table)
    }
    
    Final_summary2 <- rbind(Final_summary2, Final_summary)
  }
  
  if (ncol(Final_summary2) >= 2) {
    colnames(Final_summary2)[2] <- paste("iter", ite, sep = "_")
  }
  
  if (ite == 1) {
    Final_summary3 <- Final_summary2
  } else {
    Final_summary3 <- cbind(Final_summary3, Final_summary2[, 2, drop = FALSE])
  }
}

stopCluster(cl)

# Write results to CSV
write.csv(Final_summary3, "Permutations_1000_seed511_ctcf-vs-genes-results.csv", row.names = FALSE)

