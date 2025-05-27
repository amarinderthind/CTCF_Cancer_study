# Load necessary libraries
library(GenomicRanges)
library(dplyr)
library(ggplot2)

# Set working directory
setwd("/caculations-tmb-extract-Overlap-co-ordinates-then-mut-them-TMB/")

# Function to load BED files
read_bed <- function(file_path) {
  df <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(df) == 6) {
    colnames(df) <- c("chrom", "start", "end", "strand", "loopsite", "loopName")
  } else {
    colnames(df) <- c("chrom", "start", "end")
  }
  return(df)
}

# Load genomic region coordinates
ctcf_regions <- read_bed("motifs_with_loops_sorted_hg38_1806.bed")
utr3_regions <- read_bed("3utr_regions_hg38_Dragon.bed")
utr5_regions <- read_bed("5utr_regions_hg38_Dragon.bed")
promoter_regions <- read_bed("promoters_regions_hg38_Dragon.bed")
lnc_regions <- read_bed("LongNC_regions_hg38_Dragon.bed")
cds_regions <- read_bed("CDS_regions_hg38_Dragon.bed")

# Load mutation positions
mutation_positions <- read.csv("mutation_positions_summary_72_samples.csv")

# Function to create GenomicRanges objects
create_gr <- function(region_df, include_loopname = FALSE) {
  gr <- GRanges(
    seqnames = region_df$chrom,
    ranges = IRanges(start = region_df$start, end = region_df$end)
  )
  if (include_loopname && "loopName" %in% colnames(region_df)) {
    mcols(gr)$loopName <- region_df$loopName
  }
  return(gr)
}

# Create GRanges
ctcf_gr <- create_gr(ctcf_regions, include_loopname = TRUE)
utr3_gr <- create_gr(utr3_regions)
utr5_gr <- create_gr(utr5_regions)
promoter_gr <- create_gr(promoter_regions)
lnc_gr <- create_gr(lnc_regions)
cds_gr <- create_gr(cds_regions)
utr_gr <- reduce(c(utr3_gr, utr5_gr))  # all non-promoter features

# Combine all annotated features into one non-overlapping region
all_annotated_features <- reduce(c(utr_gr, lnc_gr, cds_gr, promoter_gr))

# Check region merging
message("Total CTCF regions: ", length(ctcf_gr))
message("Total annotated merged regions: ", length(all_annotated_features))

# Find overlaps between CTCF and annotated features
overlaps <- findOverlaps(ctcf_gr, all_annotated_features)
message("CTCF regions overlapping with annotated features: ", length(unique(queryHits(overlaps))))

# Subset CTCF sites that do not overlap with annotated features
ctcf_nonoverlapping_gr <- ctcf_gr[-queryHits(overlaps)]
message("CTCF non-overlapping regions: ", length(ctcf_nonoverlapping_gr))

# Mutation GRanges
mut_gr <- GRanges(
  seqnames = mutation_positions$CHROM,
  ranges = IRanges(start = mutation_positions$POS, end = mutation_positions$POS),
  SampleID = mutation_positions$Sample
)

# Function to compute overlaps
compute_overlap_regions <- function(ctcf_gr, feature_gr) {
  overlaps <- findOverlaps(ctcf_gr, feature_gr)
  if (length(overlaps) == 0) return(GRanges())
  overlapping_regions <- pintersect(ctcf_gr[queryHits(overlaps)], feature_gr[subjectHits(overlaps)])
  mcols(overlapping_regions) <- mcols(ctcf_gr[queryHits(overlaps)])
  return(overlapping_regions)
}

# Overlapping regions
overlapping_regions_list <- list(
  "Promoter_CTCFbs" = compute_overlap_regions(ctcf_gr, promoter_gr),
  "LncRNA_CTCFbs"   = compute_overlap_regions(ctcf_gr, lnc_gr),
  "CDS_CTCFbs"      = compute_overlap_regions(ctcf_gr, cds_gr),
  "UTR_CTCFbs"      = compute_overlap_regions(ctcf_gr, utr_gr)
)

# Check region sizes
region_sizes <- sapply(overlapping_regions_list, function(region) sum(width(region)))
ctcf_size <- sum(width(ctcf_gr))
message("Total size of all overlapping regions (bp):")
print(region_sizes)
message("Total CTCF size (bp): ", ctcf_size)

# Function to count mutations per sample
count_mutations_per_sample <- function(mut_gr, region_gr, region_name) {
  overlaps <- findOverlaps(mut_gr, region_gr)
  overlap_indices <- queryHits(overlaps)
  if (length(overlap_indices) == 0) {
    return(data.frame(SampleID = character(), Mutation_Count = integer(), Region = region_name))
  }
  mutation_counts <- as.data.frame(table(mcols(mut_gr[overlap_indices])$SampleID))
  colnames(mutation_counts) <- c("SampleID", "Mutation_Count")
  mutation_counts$Region <- region_name
  return(mutation_counts)
}

# Calculate mutation counts and TMB
mutation_df <- do.call(rbind, lapply(names(overlapping_regions_list), function(region_name) {
  region <- overlapping_regions_list[[region_name]]
  mut_count <- count_mutations_per_sample(mut_gr, region, region_name)
  mut_count$Region_Size_bp <- region_sizes[[region_name]]
  mut_count$TMB <- ifelse(mut_count$Region_Size_bp > 0,
                          (mut_count$Mutation_Count / mut_count$Region_Size_bp) * 1e6,
                          NA)
  return(mut_count)
}))



# CTCF-only (non-overlapping)
region_sizes[["CTCFbs_Only"]] <- sum(width(ctcf_nonoverlapping_gr))
ctcf_only_mutations <- count_mutations_per_sample(mut_gr, ctcf_nonoverlapping_gr, "CTCFbs_Only")
ctcf_only_mutations <- ctcf_only_mutations %>%
  mutate(Region_Size_bp = region_sizes[["CTCFbs_Only"]],
         TMB = ifelse(Region_Size_bp > 0, (Mutation_Count / Region_Size_bp) * 1e6, NA))

# Merge all mutation data
mutation_count_sample <- bind_rows(mutation_df,   ctcf_only_mutations)

# Summarise total mutations per region
total_mutations <- mutation_count_sample %>%
  group_by(Region) %>%
  summarise(Total_Mutations = sum(Mutation_Count, na.rm = TRUE), .groups = "drop")

# Check region name consistency
stopifnot(all(mutation_count_sample$Region %in% total_mutations$Region))

# Join totals
mutation_count_sample <- mutation_count_sample %>%
  left_join(total_mutations, by = "Region")

# Plot TMB
ggplot(mutation_count_sample, aes(x = Region, y = TMB, fill = Region)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_text(data = total_mutations %>% filter(!is.na(Total_Mutations)),
            aes(x = Region, y = max(mutation_count_sample$TMB, na.rm = TRUE) + 500,
                label = paste0("Total: ", Total_Mutations)),
            size = 3, vjust = -0.5, color = "red", angle = 45) +
  theme_minimal() +
  labs(title = "TMB per Region", y = "TMB (Mutations per MB)", x = "Region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(quantile(mutation_count_sample$TMB, 0.05, na.rm = TRUE), 
                           max(mutation_count_sample$TMB, na.rm = TRUE) + 1500)) +
  scale_y_continuous(breaks = seq(0, max(mutation_count_sample$TMB, na.rm = TRUE) + 1000, by = 1000))

# Define the filename
output_file <- "TMB_per_Region_plot_SHORT.png"

# Half of A4 size in inches: A4 = 8.27 x 11.69 inches → half vertically = ~5.85 inches height
# We'll use landscape half (e.g., 8.27 x 5.85 inches)
ggsave(
  filename = output_file,
  plot = last_plot(),        # or replace with your plot object
  device = "png",
  dpi = 400,
  width = 6,
  height = 4,
  units = "in",
  bg = "white"
)


library(ggpubr)
# Plot TMB with significance tests
ggplot(mutation_count_sample, aes(x = Region, y = TMB, fill = Region)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_text(data = total_mutations %>% filter(!is.na(Total_Mutations)),
            aes(x = Region, y = max(mutation_count_sample$TMB, na.rm = TRUE) + 500,
                label = paste0("Total: ", Total_Mutations)),
            size = 3, vjust = -0.5, color = "red", angle = 45) +
  stat_compare_means(method = "t.test", 
                     comparisons = combn(unique(mutation_count_sample$Region), 2, simplify = FALSE),
                     label = "p.signif", hide.ns = TRUE)+
  
  theme_minimal() +
  labs(title = "TMB per Region", y = "TMB (Mutations per MB)", x = "Region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(quantile(mutation_count_sample$TMB, 0.05, na.rm = TRUE), 
                           max(mutation_count_sample$TMB, na.rm = TRUE) + 1500)) +
  scale_y_continuous(breaks = seq(0, max(mutation_count_sample$TMB, na.rm = TRUE) + 1000, by = 1000))
