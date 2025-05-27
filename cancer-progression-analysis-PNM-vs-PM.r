---
title: "Extraction of Loop vs Mutational Matrix"
author: "Amarinder Thind"
date: "`r Sys.Date()`"
output:
  html_document: 
    default: true
  pdf_document: 
    default: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "/Users/athind/")

```

```{r}
#### DE for loop 
library(readr) 
Mut_matrix <- read_csv("Dropbox/Amarinder/Amarinder_main_projects/CTCF_motif_AT/2024-ctcf-dragen-additional-samples/loops-vs-mutational-matrix-EXCLUDE_GC_FAILED_3SAMPLES.csv")
Mut_matrix <- as.data.frame(Mut_matrix)
rownames(Mut_matrix) <- Mut_matrix$...1
Mut_matrix <- Mut_matrix[,-1]
 
#Mut_matrix[is.na(Mut_matrix)] <- 0

final_data <- Mut_matrix
```


```{r}
# Ensure metadata is read and row names are set correctly
meta <- read.csv("Dropbox/Amarinder/Sample_level_data_v11_short.csv", sep=',')
row.names(meta) <- meta$sample
meta <- subset(meta, Group %in% c("PRI_NonMet", "PRI_Met"))


## EXCLUDE - Extreme GC Bias Sample
gc_sample <- c("CSCC_0015-M1","CSCC_0023-M1","CSCC_0021-P1")
meta <- meta[!(meta$SampleName %in% gc_sample), ]
# intersect meta with Mut_matrix
meta <- meta[meta$SampleName %in% colnames(Mut_matrix), ]

 
# Sort the filtered dataframe by the 'Value' column (ascending order)
meta <- meta[order(meta$Group), ]
rownames(meta) <- meta$SampleName
```
###  

```{r}

# Filter Mut_matrix to include only the matched samples from meta
Mut_matrix <- Mut_matrix[,  meta$SampleName, drop = FALSE]

 


dim(Mut_matrix)
mut <- Mut_matrix[rowSums(Mut_matrix[,1:ncol(Mut_matrix)])>3,] 
dim(mut)

#mut <-mut[1:118,1:36] ## sort based on group

## check data is sorted
rownames(meta) == colnames(mut)
meta$Group ## check the order of the Groups

#mut1_pval <- ifelse(mut == 0, "no", "yes")
mut1_pval <- t(mut)
mut1_pval <- as.data.frame(mut1_pval)
dim(mut1_pval)

setdiff(rownames(mut1_pval), rownames(meta))

# Create named vector
named_group <- setNames(meta$Group, rownames(meta))
# Transfer group values
mut1_pval$group <- named_group[rownames(mut1_pval)]


#mut1_pval$group <- c(rep("Met", 39), rep("Pm", 13), rep("PNM", 13))

head(mut1_pval[,(ncol(mut1_pval)-10):ncol(mut1_pval)])
```


```{r}
# Function to perform chi-square or Fisher's test

compare_groups <- function(data, var) {
  tab <- table(data$group, data[, var])
  if (min(tab) >= 5) {  ## minimum 5 value in sub of confusion matrix is aspected
    test <- chisq.test(tab)
  } else {
    test <- fisher.test(tab)
  }
  return(data.frame(variable = var, p_value = test$p.value))
}


# Ensure binary data conversion (0 for no mutation, 1 for mutation)
Primary_only <- mut1_pval
var_names <- setdiff(names(Primary_only), "group")

# Convert all mutation values to binary (0 for no mutation, 1 for mutation)
Primary_only[var_names] <- lapply(Primary_only[var_names], function(x) ifelse(x >= 1, 1, 0))

# Perform statistical tests between primary and metastasis groups
results <- lapply(var_names, function(var) compare_groups(Primary_only, var))
results <- do.call(rbind, results)

# Adjust p-values using BH method
results$p_adjusted <- p.adjust(results$p_value, method = "BH")

# Combine data and results
Primary_only <- t(Primary_only)
combined_primary <- cbind(Primary_only[-nrow(Primary_only), ], results)

# Ensure column names match valid samples in meta
primary_samples <- meta$SampleName[meta$Group == "PRI_NonMet" & !is.na(meta$SampleName)]
metastasis_samples <- meta$SampleName[meta$Group == "PRI_Met" & !is.na(meta$SampleName)]

# Identify sample columns in the combined data
primary_cols <- which(colnames(combined_primary) %in% primary_samples)
metastasis_cols <- which(colnames(combined_primary) %in% metastasis_samples)

# Check column validity
stopifnot(length(primary_cols) > 0)
stopifnot(length(metastasis_cols) > 0)

# Extract primary and metastasis columns and ensure binary format
temp_primary <- combined_primary[, primary_cols, drop = FALSE]
temp_metastasis <- combined_primary[, metastasis_cols, drop = FALSE]

temp_primary <- data.frame(lapply(temp_primary, as.numeric))
temp_metastasis <- data.frame(lapply(temp_metastasis, as.numeric))

# Compute row sums for mutation counts (0 or 1 as binary)
combined_primary$sumPNM <- rowSums(temp_primary, na.rm = TRUE)
combined_primary$sumPM <- rowSums(temp_metastasis, na.rm = TRUE)

# Verify binary sums and inspect final results
head(combined_primary[, c("sumPNM", "sumPM")])


```
```{r}
print("Range of Primary no Met samples mutations in different loop")
 range(combined_primary$sumPNM)

 print("Range of primary Metastases samples mutations in different loop")
 range(combined_primary$sumPM)

```


```{r}
#write.csv(combined_primary,"Dropbox/Amarinder/Amarinder_main_projects/CTCF_motif_AT/2024-ctcf-dragen-additional-samples/step2-extraction-of-loop-mutations-matrix/DE_loops_PNM_vs_MetALL_combined-Dec2024.csv")
```


```{r}
#  subset top 10 loops based on p-value 

# library(dplyr) 
# 
# top_loops <- combined_primary %>% 
#   arrange(p_value) %>% 	#Arrange rows by padj values
#   filter(p_value < 0.2) %>%   #filter based on logFC
#   pull(variable) 



library(dplyr)

# Assuming combined_primary is your data frame
top_loops <- combined_primary %>%
  arrange(p_value) %>%  # Arrange rows by p_value
  filter( p_value <0.3 & (sumPNM <= 1 | sumPM <= 1)) %>%  # Filter based on p_value and absence in PNM or MetALL
  pull(variable)  # Select the 'variable' column

top_loops
top_loops

## JUSt to show here

temp1 <- combined_primary %>% 
  arrange(p_value) %>% 	#Arrange rows by padj values
  filter( p_value <0.3 & (sumPNM <= 1 | sumPM <= 1))  

temp1[,c("sumPNM","sumPM")]
```
```{r}
# If loop names are in rownames
temp1$loop_id <- rownames(temp1)

library(tidyr)
library(ggplot2)

# Reshape for plotting
plot_data <- temp1 %>%
  pivot_longer(cols = c(sumPNM, sumPM), names_to = "Category", values_to = "Count")

plot_data$Count_mod <- ifelse(plot_data$Category == "sumPNM", -plot_data$Count, plot_data$Count)

# Calculate max value for symmetrical axis
max_count <- ceiling(max(abs(plot_data$Count_mod)))

# Plot
p <- ggplot(plot_data, aes(x = reorder(loop_id, -abs(Count_mod)), y = Count_mod, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(x = "Loop ID", y = "Count", title = "PNM vs PM Signal in DE Loops") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  coord_flip() +
  scale_y_continuous(
    breaks = seq(-max_count, max_count, by = 1),
    labels = abs
  )

print(p)



#ggsave("Dropbox/Amarinder/Amarinder_main_projects/CTCF_motif_AT/Manuscript_abstract/DEloop_PNM-vsPM.png", plot = p, dpi = 600, width = 4, height = 3, units = "in", bg="white")

```


```{r}
```
 


```{r}
loops <- combined_primary

colnames(loops)
end <- ncol(loops)-5 ## remove last five entries contain non-samples cols
loops2 <- loops[,1:end]
colnames(loops2)

print(paste("Total number of WGS samples :",nrow(meta),"out of ",ncol(loops2)))
table(meta$Group)

common <- intersect(meta$wgs_samples,meta$totalRNAseq) 
common

meta <- subset(meta, SampleName %in% common)
```


```{r}
print(paste("Total number of RNAseq samples :",nrow(meta),"out of ",ncol(loops2)))
table(meta$Group)
```



```{r}
library(DESeq2)

## Prepare the gene list for the selected loop (calling it top loop here)
gene_list_loop <- read.csv("Dropbox/Amarinder/Amarinder_main_projects/CTCF_motif_AT/2024-ctcf-dragen-additional-samples/Step1.b-extraction-of-loops-gene-list/Loops_with_gene_list_combine-inside-and-1000bp-outside-.csv")

# Remove trailing commas from the 'genes' column
gene_list_loop$hgnc_symbol <- gsub("\\s+", "", gene_list_loop$hgnc_symbol)
gene_list_loop$ensembl_gene_id <- gsub("\\s+", "", gene_list_loop$ensembl_gene_id)



 
all_genes <- read.csv("Dropbox/Amarinder/Amarinder_main_projects/CTCF_motif_AT/2024-ctcf/all_gene_biomart.csv") 

library(dplyr)
library(stringr)
library(purrr)

# Clean lookup table (remove Ensembl version suffix like ".1")
all_genes <- all_genes %>%
  mutate(
    ensembl_gene_id = str_trim(str_remove(as.character(ensembl_gene_id), "\\.\\d+$")),
    hgnc_symbol = as.character(hgnc_symbol)
  ) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# Named vector: ensembl_id -> hgnc_symbol
ensembl_to_hugo <- setNames(all_genes$hgnc_symbol, all_genes$ensembl_gene_id)

# Track unmapped IDs
unmapped_ids <- character()

# Safely map gene IDs
gene_list_loop <- gene_list_loop %>%
  mutate(gene2merged = map_chr(ensembl_gene_id, function(x) {
    if (is.na(x) || x == "") return(NA_character_)
    
    # Split and clean each ID
    ids <- str_split(x, ",")[[1]] %>%
      str_trim() %>%
      str_remove("\\.\\d+$")

    # Map each ID
    mapped <- map_chr(ids, function(id) {
      if (id %in% names(ensembl_to_hugo) && !is.na(ensembl_to_hugo[id])) {
        return(ensembl_to_hugo[id])
      } else {
        unmapped_ids <<- unique(c(unmapped_ids, id))  # collect missing ones
        return(id)  # fallback to Ensembl ID
      }
    })

    str_c(mapped, collapse = ",")
  }))

# Print summary of unmapped Ensembl IDs
if (length(unmapped_ids) > 0) {
  cat("⚠️ Unmapped Ensembl IDs (showing up to 20):\n")
  print(head(unmapped_ids, 20))
  cat(sprintf("Total unmapped IDs: %d\n", length(unmapped_ids)))
}




 ###### Prepare LOOP GENE LIST
intersect_genes <- gene_list_loop[gene_list_loop$identifier %in% top_loops,]

print("Top DE loops with gene names")
print(intersect_genes[colnames(intersect_genes) %in% c( "identifier","gene2merged")])



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
```


```{r}
library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare long-format mutation data for top loops
temp1_long <- temp1 %>%
  pivot_longer(cols = c(sumPM, sumPNM), names_to = "Group", values_to = "MutCount") %>%
  mutate(Group = recode(Group, sumPM = "PRI_Met", sumPNM = "PRI_NonMet"))

# Compute number of genes per loop
loop_gene_counts <- gene_list_loop %>%
  filter(identifier %in% temp1$loop_id) %>%
  mutate(num_genes = sapply(strsplit(gsub("\\s+", "", ensembl_gene_id), ","), length)) %>%
  group_by(identifier) %>%
  summarise(n_genes = sum(num_genes))

# Merge gene counts into long mutation data
temp1_long <- temp1_long %>%
  left_join(loop_gene_counts, by = c("loop_id" = "identifier"))

# Plot grouped bars
p <- ggplot(temp1_long, aes(x = reorder(loop_id, -MutCount), y = MutCount, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(
    aes(
      x = loop_id,
      y = 5.5,  # Position text slightly above the bar
      label = ifelse(Group == "PRI_Met", paste0("Genes: ", n_genes), "")
    ),
    position = position_dodge(width = 0),
    vjust = 0,  # Adjust vertical position relative to the bar
    size = 3.5,
    color = "red",
    angle = 90  # Rotate text by 90 degrees
  ) +
  ylim(0, 6.5) +
  labs(
    title = "Top DE Loops with Mutation Counts and Gene Count per Loop",
    x = "Loop ID",
    y = "Mutation Count",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)


#ggsave("Dropbox/Amarinder/Amarinder_main_projects/CTCF_motif_AT/Manuscript_abstract/DEloop_PNM-vsPM_with_genes_info.png", plot = p, dpi = 600, width = 4, height = 3, units = "in", bg="white")

```


```{r}
 

# Specify the order of loop identifiers
ordered_loops <- c("loop_818", "loop_937", "loop_70", "loop_390", "loop_284", 
                   "loop_1083", "loop_95", "loop_657", "loop_53", "loop_45", "loop_115")

# Reorder intersect_genes based on the given order of loop_id
ordered_intersect_genes <- intersect_genes %>%
  filter(identifier %in% ordered_loops) %>%
  mutate(identifier = factor(identifier, levels = ordered_loops)) %>%
  arrange(identifier)

library(gridExtra)
library(grid)  # <-- Load this to use gpar()



# View the reordered dataset with loop_id and hgnc_symbol
shorttb <-ordered_intersect_genes[, c("identifier", "hgnc_symbol")]
# Create the tableGrob
table_plot <- tableGrob(shorttb)

# Loop over data rows (skip header)
for (i in 1:nrow(shorttb)) {
  row_index <- i + 1  # tableGrob's header is row 1
  col_index <- 2      # second column

  # Find the index in layout that matches the cell
  layout_match <- which(table_plot$layout$t == row_index & table_plot$layout$l == col_index)

  if (length(layout_match) > 0) {
    grob_name <- table_plot$layout$name[layout_match]
    grob_index <- which(names(table_plot$grobs) == grob_name)
    
    if (length(grob_index) > 0) {
#      table_plot$grobs[[grob_index]]$gp <- gpar(fontface = "italic")
      table_plot$grobs[[grob_index]]$gp <- gpar(fontface = "italic", fontfamily = "Times")

    }
  }
}

# 
# Draw the table
grid.newpage()
grid.draw(table_plot)
# Save the table as an image
ggsave("Dropbox/Amarinder/Amarinder_main_projects/CTCF_motif_AT/Manuscript_abstract/ordered_intersect_genes_table.png", table_plot, width = 10, height = 8, dpi = 300)
dev.off()
 
# Save to PDF to verify if italics render properly
grid.newpage()
pdf("Dropbox/Amarinder/Amarinder_main_projects/CTCF_motif_AT/Manuscript_abstract/ordered_intersect_genes_table.png", width = 10, height = 8)
grid.draw(table_plot)
dev.off()

```


```{r}
library(readr)
rawcount_ori <- read_csv("Dropbox/Amarinder/Amarinder_main_projects/CTCF_motif_AT/2024-ctcf-dragen-additional-samples/feature_count-v3-80samples/ftype-exon-fattr-gene_id/countmatrix_genes_Modified_samp_names-unique-ensemble-id-v2.csv")
# Set the first column as rownames
rawcount_ori <- as.data.frame(rawcount_ori)
rownames(rawcount_ori) <- rawcount_ori[[1]]
rawcount_ori <- rawcount_ori[, -1]


## Replace NAs by zero and changing the input to required format
#rawcount_ori <- round(rawcount_ori) 
rawcount_ori[is.na(rawcount_ori)] <- 0

# Discard genes expressed in less than 15% of all samples
keep <- rowSums(rawcount_ori > 0) >= round(ncol(rawcount_ori) * 0.15)
rawcount <- rawcount_ori[keep,]



colnames(rawcount)

```




```{r}

#rawcount <- rawcount[rownames(rawcount) %in% all_Loopgenes,]
head(rawcount)
#rawcountCoding <- rawcount[row.names(rawcount) %in%  all_coding_genes$hgnc_symbol,]

```


```{r}
######################################3
#anno <- read_csv("Dropbox/Amarinder/Sample_level_data_v11_short.csv")

anno <- meta

intersect(colnames(rawcount), anno$SampleName)
setdiff(colnames(rawcount), anno$SampleName)
setdiff( anno$SampleName,colnames(rawcount))

firstC<-"PRI_NonMet"       #case1 #case2 #case3 etc          
SecondC <-"PRI_Met"   

anno <- subset(anno, Group == firstC | Group == SecondC ) #| group == "Primary (To be updated)"


## discard sample CSCC_0023-M1
anno <- anno[anno$SampleName != "CSCC_0023-M1", ]

rawcount <- rawcount[,names(rawcount) %in% anno$SampleName]



anno <- anno[match(colnames(rawcount), anno$SampleName),] ## reordering anno row with colnmaes of rawcount
rownames(anno) <- anno$SampleName
```


```{r}

library(PCAtools)
library(edgeR)
logRaw <- as.matrix(log2(rawcount+1))
top1000.order <- head(order(matrixStats::rowVars(logRaw), decreasing = TRUE), 1000)
p <- PCAtools::pca(mat = logRaw[top1000.order,], metadata = anno, removeVar = 0.1)


biplot(p,
       lab = paste0(p$metadata$Group),
       # colby = 'group',
       hline = 0, vline = 0,
       legendPosition = 'right',
       encircle = T )
```


```{r}
table(anno$Group)

dds <- DESeqDataSetFromMatrix(countData = round(rawcount), colData = anno, design = ~totalRNA_Batch+ Group )
#dds <- DESeqDataSetFromMatrix(countData = round(rawcount), colData = anno, design = ~ Group )

# then run DESeq()
## Run DESEQ2
dds <- DESeq(dds)

##ensure your data is a good fit for the DESeq2 model
plotDispEsts(dds)
```


```{r}

#Then, we can extract the normalized count values for these top 20 genes:

## normalized counts for top 20 significant genes
#normalized_counts$gene <- row.names(normalized_counts)  
#dds1 <- estimateSizeFactors(dds)

vst <- varianceStabilizingTransformation(dds, blind = TRUE)
  ### Transform counts for data visualization #options (1) vst (2) rld

normalized_counts <- as.data.frame(assay(vst)) ##VST

```

```{r}
#In case of multiple comparisons ## we need to change the contrast for every comparision

contrast<- c("Group",firstC,SecondC)

res <- results(dds, contrast=contrast)

#View(data.frame(res))
```


```{r}
###############
library(tidyverse)
library(ggplot2)
library(ggrepel)
#library(DEGreport)
library(RColorBrewer)
library(pheatmap)

res_df <-data.frame(res)
### change names to hugo
res_df <- res_df[row.names(res_df) %in%  all_Loopgenes,]
res_df$gene <- row.names(res_df)


```








```{r}


##################### Subset the data frame
#filtered_df <- res_df[res_df$gene %in% all_Loopgenes, ]
 
filtered_df <- res_df
# View the filtered DataFrame
print(filtered_df[,c("gene", "log2FoldChange", "pvalue", "padj")])
```



```{r}

colnames(filtered_df)[colnames(filtered_df) == "gene"] <- "ensembl_gene_id"
filtered_df <- filtered_df %>%
  left_join(all_genes %>% select(hgnc_symbol, ensembl_gene_id), by = "ensembl_gene_id")

filtered_df <- filtered_df %>%
  mutate(hgnc_symbol = ifelse(hgnc_symbol == "", ensembl_gene_id, hgnc_symbol))

# Check the result
head(filtered_df)
 
filtered_df2 <- filtered_df
```

```{r}

# Recalculate adjusted p-values (padj) for the filtered genes
#filtered_df$new_padj <- p.adjust(filtered_df$pvalue, method = "BH")

# Assuming your DataFrame is named 'df'
 

filtered_df <- filtered_df %>%
  filter((log2FoldChange > 1 | log2FoldChange < -1) & pvalue < 0.01)



# View the filtered DataFrame
print(filtered_df[,c("hgnc_symbol", "log2FoldChange", "pvalue", "padj")])
```

```{r}
library(ggplot2)
library(reshape2)
library(dplyr)
library(pheatmap)

# Ensure rownames are gene names for better visualization
selected_genes <- filtered_df$ensembl_gene_id
subset_counts <- normalized_counts[rownames(normalized_counts) %in% selected_genes, ]
rownames(subset_counts) <- filtered_df$hgnc_symbol[match(rownames(subset_counts), filtered_df$ensembl_gene_id)]
subset_counts <- subset_counts[, meta$SampleName]  # Ensure sample order matches meta

# Reshape for ggplot
melted_counts <- as.data.frame(subset_counts) %>%
  rownames_to_column(var = "Gene") %>%
  melt(id.vars = "Gene", variable.name = "Sample", value.name = "Expression")

# Join with metadata
melted_counts <- left_join(melted_counts, meta, by = c("Sample" = "SampleName"))

# Create separate facet boxplots for each gene
p1 <- ggplot(melted_counts, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(
    title = "Gene Expression Levels per Sample Group",
    x = "Sample Group",
    y = "Expression Level"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("PRI_Met" = "red", "PRI_NonMet" = "blue")) +
  facet_wrap(~Gene, scales = "free_y", ncol = 3)

# Scale counts for the heatmap visualization
scaled_counts <- t(scale(t(subset_counts)))

# Create heatmap with sample annotations
annotation_df <- data.frame(Group = meta$Group)
rownames(annotation_df) <- meta$SampleName
annotation_colors <- list(Group = c("PRI_Met" = "red", "PRI_NonMet" = "blue"))

pheatmap(
  scaled_counts,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  annotation_col = annotation_df,
  annotation_colors = annotation_colors,
  main = "Heatmap of Selected Gene Expression",
  fontsize_row = 10,
  fontsize_col = 8,
  border_color = NA
)

pheatmap(
  scaled_counts,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  annotation_col = annotation_df,
  annotation_colors = annotation_colors,
  main = "Heatmap of Selected Gene Expression",
  fontsize_row = 10,
  fontsize_col = 8,
  border_color = NA
)

# Print the boxplot separately
print(p1)
```


```{r}
library(ggplot2)

# library(ggplot2)
library(dplyr)

# Ensure `pvalue` does not contain NA and is log-transformed correctly
filtered_df2 <- filtered_df2 %>%
  mutate(
    neg_log10_pval = -log10(pvalue),
    significant = ifelse(pvalue < 0.01 & abs(log2FoldChange) > 1, "Significant", "Not Significant")
  )

# Create volcano plot using p-values instead of adjusted p-values
library(ggplot2)
library(ggrepel)

p <- ggplot(filtered_df2, aes(x = log2FoldChange, y = neg_log10_pval, color = significant)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  labs(
    title = "Volcano Plot",
  x = "Log2 Fold Change",
  y = "-Log10 p-value",
  color = "Significance"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  
  # Add gene labels for significant points
  geom_text_repel(
    data = subset(filtered_df2, significant == "Significant"),
    aes(label = hgnc_symbol),
    size = 3,
    box.padding = 0.5,
    max.overlaps = 10
  )

print(p)
# Save the table as an image
#ggsave("Dropbox/Amarinder/Amarinder_main_projects/CTCF_motif_AT/Manuscript_abstract/vacano.png", plot = p, width = 4, height = 3, dpi = 600, bg="white")

```

