### Genomic Regions and BackGround mutational density comparison using various scaling factors
## Either Upstream background or Downstream Background (Need to define)

library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(tidyr)

getwd()
setwd("/motif-position-analysis/") ## Path for working directory

# Read region coordinates (e.g. 3' UTR , in BED format) one at a time

regions <- read.table("../oncdriveFML_regions_for_coverage/3utr_regions_hg38_Dragon.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(regions) <- c("chrom", "start", "end", "region")

## header names differs for CTCF co-ordinates so using following option
#colnames(regions) <- c("chrom", "start", "end", "strand","utr","loop") ## if using CTCF, go to next PLOT DIRECTLY

# Read mutation positions (multiple samples) # File prepared using script named "motif-CTCFbs-specific-plots.R"
mutation_positions <- read.csv("mutation_positions_summary_72_samples.csv")

# Calculate the length (width) of each region
regions$length <- regions$end - regions$start
print(range(regions$length))  # Print the range of lengths to understand the distribution

# Remove regions with length <= 3
regions <- regions[regions$length > 3 ,]

# Plot the density of region lengths
ggplot(regions, aes(x = length)) +
  geom_density(fill = "#1f77b4", color = "black", alpha = 0.7) +
  labs(
    title = "Density of Region Lengths",
    x = "Region Length (bp)",
    y = "Density"
  ) +
  theme_minimal()

# Convert regions into a GenomicRanges object
region_gr <- GRanges(seqnames = regions$chrom,
                     ranges = IRanges(start = regions$start, end = regions$end) )

# Create a GenomicRanges object for mutations
mut_gr <- GRanges(seqnames = mutation_positions$CHROM,
                  ranges = IRanges(start = mutation_positions$POS, end = mutation_positions$POS))

####################################################################
# Step: Define the surrounding regions with multiple scaling factors
####################################################################
library(data.table)
library(BiocParallel)

scaling_factors <- c(10,50,100,200,400,600,800)

# Convert mutation_positions to data.table
mutation_positions_dt <- as.data.table(mutation_positions)

# Placeholder for density results
all_surrounding_density <- list()

# Batch processing of regions to reduce memory usage
batch_size <- 500  # Adjust based on available memory
num_batches <- ceiling(length(region_gr) / batch_size)

# Predefined chromosome lengths based on GRCh38.p14
chromosome_lengths <- c(
  chr1 = 248956422, chr2 = 242193529, chr3 = 199396827, chr4 = 191154276,
  chr5 = 180915260, chr6 = 171115067, chr7 = 159138663, chr8 = 146364022,
  chr9 = 141213431, chr10 = 135534747, chr11 = 134452384, chr12 = 132289534,
  chr13 = 114127980, chr14 = 106360585, chr15 = 100338915, chr16 = 88822254,
  chr17 = 78654742, chr18 = 76117153, chr19 = 63806651, chr20 = 62435965,
  chr21 = 46944323, chr22 = 49528953, chrX = 156040895, chrY = 57274415
)

compute_density_for_batch <- function(batch_idx, region_gr, mutation_positions, scaling_factors, chromosome_lengths) {

  # Subset regions for the current batch
  batch_start <- (batch_idx - 1) * batch_size + 1
  batch_end <- min(batch_idx * batch_size, length(region_gr))
  region_gr_batch <- region_gr[batch_start:batch_end]
  
  batch_results <- list()
  
  for (scale in scaling_factors) {
    
    # Compute valid start and end ranges for the surrounding regions           ## DOWNSTREAM
    start_pos <- pmax(end(region_gr_batch) + 1, 1)
    end_pos <- pmin(
      chromosome_lengths[as.character(seqnames(region_gr_batch))],
      end(region_gr_batch) + (scale * width(region_gr_batch))
    )
    
    ## for upstream regions                                                     ## UPSTREAM
   # start_pos <- pmax(start(region_gr_batch) - (scale * width(region_gr_batch)), 1)
    #end_pos <- pmax(start(region_gr_batch) - 1, 1)  # Ensure valid ends
    
    # Ensure the ranges are valid
    valid_indices <- which(start_pos <= end_pos)
    if (length(valid_indices) == 0) next
    
    surrounding_gr <- GRanges(
      seqnames = seqnames(region_gr_batch)[valid_indices],
      ranges = IRanges(start = start_pos[valid_indices], end = end_pos[valid_indices])
    )
    
    # Compute overlaps
    surrounding_mutations <- findOverlaps(mut_gr, surrounding_gr)
    surrounding_hits <- unique(queryHits(surrounding_mutations))
    
    if (length(surrounding_hits) > 0) {
      # Filter mutation positions for valid hits
      temp_df <- mutation_positions[mutation_positions$POS %in% mutation_positions$POS[surrounding_hits], ]
      surrounding_sample_density <- aggregate(POS ~ Sample, data = temp_df, FUN = length)
      
      total_width <- sum(width(surrounding_gr))
      surrounding_sample_density$surrounding_density <- (surrounding_sample_density$POS / total_width) * 1e6
      surrounding_sample_density$scaling_factor <- paste(scale, "x width")
    } else {
      surrounding_sample_density <- data.frame(
        Sample = character(), POS = integer(), surrounding_density = numeric(), scaling_factor = character()
      )
      surrounding_sample_density$scaling_factor <- paste(scale, "x width")
    }
    
    batch_results[[paste(batch_idx, scale, sep = "_")]] <- surrounding_sample_density
  }
  
  return(batch_results)
}


# Running the sequential loop instead of parallel
all_surrounding_density <- list()

for (batch_idx in seq_len(num_batches)) {
  batch_results <- compute_density_for_batch(batch_idx, region_gr, mutation_positions, scaling_factors, chromosome_lengths)
  all_surrounding_density <- c(all_surrounding_density, batch_results)  # Append results
}

# After loop, combine the results into a single dataframe
combined_surrounding_density <- do.call(rbind, all_surrounding_density)


# Calculate the mean mutation density for each surrounding region across all samples
mean_surrounding_density <- combined_surrounding_density %>%
  group_by(scaling_factor) %>%
  summarise(mean_density = mean(surrounding_density), .groups = "drop")

# Define "is_in_region" for mutations in the target region
region_hits <- findOverlaps(mut_gr, region_gr)
region_mutation_indices <- queryHits(region_hits)

# Annotate mutations as being in the region
mutation_positions <- mutation_positions %>%
  mutate(is_in_region = POS %in% mutation_positions$POS[region_mutation_indices])

# Calculate the mean mutation density for the actual region

region_sample_density <- mutation_positions %>%
  filter(is_in_region == TRUE) %>%
  group_by(Sample) %>%
  summarise(region_mutations = n(), .groups = "drop") %>%
  mutate(region_density = (region_mutations / sum(width(region_gr))) * 1e6)  # Per Mb unit

# Calculate the mean region density across all samples
mean_region_density <- mean(region_sample_density$region_density)

# Add the region density as a comparison entry to the result

mean_surrounding_density <- bind_rows(
  mean_surrounding_density,
  data.frame(scaling_factor = "Region", mean_density = mean_region_density)
)

# Extract scaling factors as numbers, excluding "Region"
scaling_values <- gsub(" x width", "", mean_surrounding_density$scaling_factor)
scaling_values <- as.numeric(scaling_values[!scaling_values %in% "Region"])

# Sort and create proper levels, placing "Region" at the beginning
sorted_levels <- c("Region", paste0(sort(scaling_values), " x width"))

# Assign ordered factor levels
mean_surrounding_density$scaling_factor <- factor(mean_surrounding_density$scaling_factor, 
                                                  levels = sorted_levels)

# Visualize the comparison of mean mutation density across the region and different scaling factors

ggplot(mean_surrounding_density, aes(x = scaling_factor, y = mean_density, fill = scaling_factor)) +
  geom_bar(stat = "identity") +
  labs(
    #title = "Mean Mutation Density in Region vs Surrounding Regions at Different Scaling Factors",
    x = "Region Type (Scaling Factor)",
    y = "mutations per Mb"
  ) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", 
                               "#9467bd", "#8c564b", "#e377c2", "#17becf","black","green","gray","blue","red","orange")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

