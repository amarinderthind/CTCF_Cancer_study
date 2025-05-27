### Script is for 
## Mutation density Across ±1kb CTCF Binding Sites and CTCFbs itself (wihtout substitution type) and plot 2 per sample multiplet
## Length of loops distribution
## TMB CTCFbs vs Surround 1k (Stack) per sample (one plot not multiplet)

################################################################################################3#

## For coputaional purpose files were loaded and saved before to use it again and again

# # Step: List all VCF files in your directory (assuming all VCFs are in the same folder)
# #vcf_files <- list.files(path = "../step1-intersect-vcf-with-loop-motif/", pattern = "\\.vcf$", full.names = TRUE)
# 
# vcf_files <- list.files(path = "/mnt/NAS_rrg/cSCC-omics-data/DNA/WGS-illumina/purple-mutech-pipeline-outputs/Vcf-first-pipeline-mutect2-purple/ALL_WGS_purple_somatic_till_feb_2023_except_CL", pattern = "\\.vcf$", full.names = TRUE)
# vcf_files
# 
# # Initialize an empty data frame to store mutation data
# mutation_positions <- data.frame()
# 
# # Step: Loop through each VCF file and extract mutation positions
# for (vcf_file in vcf_files) {
#   
#   print("reading file: ")
#   print(vcf_file)
#   # Read the VCF file using vcfR
#   vcf_data <- read.vcfR(vcf_file)
#   
#   # Step 6: Extract the variant positions (CHROM, POS) from the VCF file
#   vcf_variants <- data.frame(
#     CHROM = vcf_data@fix[, "CHROM"],
#     POS = as.integer(vcf_data@fix[, "POS"]),
#     REF = vcf_data@fix[, "REF"],  # REF column
#     ALT = vcf_data@fix[, "ALT"]   # ALT column
#     
#   )
#   
#   # Add the sample name to the variant data (sample name is derived from the VCF filename)
#   sample_name <- gsub("\\.vcf$", "", basename(vcf_file))
#   sample_name = gsub("_vs.*", "", sample_name)
#   vcf_variants$Sample <- sample_name
#   
#   # Append the extracted mutation data for this sample to the main data frame
#   mutation_positions <- rbind(mutation_positions, vcf_variants)
# }
# 
# 
# head(mutation_positions)
# unique(mutation_positions$Sample)
# 
# samples_to_exclude <- c("CSCC_0015-M1",
#                         "CSCC_0023-M1",
#                         "CSCC_0021-P1")  # Specify the sample identifiers to exclude
# 
# mutation_positions <- mutation_positions[!(mutation_positions$Sample %in% samples_to_exclude),]
# unique(mutation_positions$Sample)
# 
# #write.csv(mutation_positions,"mutation_positions_summary_initial-ctcfbs-surrenounding_72_samples.csv")
# #write.csv(mutation_positions,"mutation_positions_PURPLE_49samples_reading_via_vcfR_vcfData.csv")

######################################################################################################################
###############################
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(vcfR)
library(GenomicRanges)
library(tidyr)

setwd("/home/singh/Dropbox/Amarinder/Amarinder_main_projects/CTCF_motif_AT/2024-ctcf-dragen-additional-samples/")

# Step: Read CTCF region coordinates (e.g., in BED format)
ctcf_regions <- read.table("motifs_with_loops_sorted_hg38_1806.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(ctcf_regions) <- c("chrom", "start", "end","strand","5or3prime","loop")


# Step: Sort data by loop and start position to ensure correct order
df_sorted <- ctcf_regions %>%
  arrange(loop, start)

# Step: Calculate the difference between the end of one CTCF binding site and the start of the next one
df_sorted <- df_sorted %>%
  group_by(loop) %>%  # Group by loop to calculate differences within each loop
  mutate(
    next_start = lead(start),  # Get the start of the next CTCF binding site
    diff = next_start - end  # Calculate the difference between end of current and start of next
  ) %>%
  ungroup()  # Remove grouping after mutation


# Assuming df_sorted is the data with the diff column calculated (from previous steps)

library(scales)

# Step 1: Plot the distribution of loop lengths (differences) as a density plot with custom x-axis labels
ggplot(df_sorted, aes(x = diff)) +
  geom_density(fill = "skyblue", color = "black", alpha = 0.7) +  # Fill with color and outline in black
  scale_x_continuous(
    breaks = seq(50000, max(df_sorted$diff, na.rm = TRUE), by = 100000),  # Start from 50,000
    labels = label_number(scale = 1e-3, suffix = "k"),  # Convert numbers to "k" format
    limits = c(50000, max(df_sorted$diff, na.rm = TRUE))  # Set x-axis to start from 50,000
  ) +
  labs(
    title = "Density of Loop Lengths (Differences between CTCF Binding Sites)",
    x = "Loop Length (bp)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,size = 10)  # Rotate x-axis labels by 45 degrees
  )



###############################

# Convert CTCF regions into a GenomicRanges object
ctcf_gr <- GRanges(seqnames = ctcf_regions$chrom,
                   ranges = IRanges(start = ctcf_regions$start, end = ctcf_regions$end))

# Step: Define the surrounding region ±1kb (1000 bp) around each CTCF region
surrounding_gr <- GRanges(
  seqnames = seqnames(ctcf_gr),
  ranges = IRanges(
    start = start(ctcf_gr) - 990,  # 1kb upstream of CTCF
    end = end(ctcf_gr) + 990      # 1kb downstream of CTCF
  )
)

# Step: Merge overlapping regions (CTCF regions + surrounding ±1kb)
merged_gr <- reduce(c(ctcf_gr, surrounding_gr))

##################
## Load mutational data ##################

mutation_positions <- read.csv("motif-position-analysis/mutation_positions_summary_72_samples.csv")
mutation_positions2 <- mutation_positions
 



# Step 6: Create a GenomicRanges object for mutations
mut_gr <- GRanges(seqnames = mutation_positions$CHROM,
                  ranges = IRanges(start = mutation_positions$POS, end = mutation_positions$POS))

# Step 7: Find mutations in merged CTCF + surrounding regions
overlaps <- findOverlaps(mut_gr, merged_gr)
#View(as.data.frame(overlaps)) ## this has index of mut_gr && merged_gr which hit against each other
#View(as.data.frame(merged_gr))
#View(as.data.frame(mut_gr))

# Extract mutation positions that fall within the CTCF regions
mutations_in_ctcf <- mutation_positions[queryHits(overlaps), ]
#View(as.data.frame(mutations_in_ctcf))

# Normalize positions to the ±1kb scale
normalized_positions <- data.frame(
  CHROM = seqnames(merged_gr)[subjectHits(overlaps)],
  Normalized_POS = start(mut_gr)[queryHits(overlaps)] - start(merged_gr)[subjectHits(overlaps)] - 1000,
  Sample = mutations_in_ctcf$Sample
)


#################### PLOT 1
# Step 8: Summarize mutation recurrence
total_samples <- length(unique(normalized_positions$Sample))  # Count the unique sample IDs

recurrence_summary <- normalized_positions %>%
  group_by(Normalized_POS) %>%
  summarise(
    Recurrence = n(),
    Total_Loops = 1806 * total_samples,  # Adjust total loops by number of samples
    Recurrence_Rate = Recurrence / Total_Loops,
    Recurrence_Rate_per_mb = (Recurrence / Total_Loops) * 1e6  # Correct recurrence rate per Mb
  )


# Step 9: Visualize the mutation recurrence across the ±1kb regions
p <- ggplot(recurrence_summary, aes(x = Normalized_POS, y = Recurrence_Rate_per_mb)) +
  geom_line(color = "blue", linewidth = 1) +
  labs(
    title = "Mutation density Across ±1kb CTCF Binding Sites",
    x = "Position Relative to center of CTCF Binding Site (bp)",
    y = "Mean mutation per Mb"
  ) +
  theme_minimal()


# Save the plot with high resolution and custom size
ggsave("motif-position-analysis/out-figures/1000bp-ctcf-mean-all-samples-per-mb.png", plot = p, dpi = 600, width = 7, height = 5, units = "cm", bg="white")

################ PLOT 2
# Step 8: Summarize mutation recurrence for each sample
recurrence_summary2 <- normalized_positions %>%
  group_by(Sample, Normalized_POS) %>%
  summarise(
    Recurrence = n(),  # Count of mutations at this position
    Total_Loops = 1806,  # NO NEED to adjust per sample; plot is sample based,
    Recurrence_Rate = Recurrence / Total_Loops,  # Calculate the recurrence rate
   Recurrence_Rate_per_mb = (Recurrence / Total_Loops) * 1e6 
  )

# Step 9: Visualize the mutation recurrence across the ±1kb regions for each sample
# ggplot(recurrence_summary2, aes(x = Normalized_POS, y = Recurrence_Rate, color = Sample, group = Sample)) +
#   #geom_line(size = 1) +
#   geom_point(size = 2, alpha = 0.7) +
#   labs(
#     title = "Mutation Recurrence Across ±1kb CTCF Binding Sites",
#     x = "Position Relative to CTCF Binding Site (bp)",
#     y = "Recurrence Rate",
#     color = "Sample"
#   ) +
#   theme_minimal() +
#   theme(legend.position = "none")


ggplot(recurrence_summary2, aes(y = Normalized_POS, x = Sample, fill = Recurrence_Rate_per_mb)) +
  geom_tile() +
  scale_fill_gradient(
    low = "white", high = "blue", trans = "log", 
    breaks = scales::trans_breaks("log10", function(x) 10^x), 
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  labs(
    title = "Mutation Recurrence Across ±1kb CTCF Binding Sites",
    y = "Position Relative to CTCF Binding Site (bp)",
    x = "Sample",
    fill = "Recurrence Rate per mb (log-scaled)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggplot(recurrence_summary2, aes(y = Normalized_POS, x = Sample, fill = Recurrence_Rate_per_mb)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = median(recurrence_summary2$Recurrence_Rate_per_mb, na.rm = TRUE)
  ) +
  labs(
    title = "Mutation Recurrence Across ±1kb CTCF Binding Sites",
    y = "Position Relative to CTCF Binding Site (bp)",
    x = "Sample",
    fill = "Recurrence Rate"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

### multple plots
# ggplot(recurrence_summary2, aes(x = Normalized_POS, y = Recurrence_Rate_per_mb)) +
#   geom_line(size = 1) +
#   facet_wrap(~Sample, scales = "free_y") +
#   labs(
#     title = "Mutation Recurrence Across ±1kb CTCF Binding Sites by Sample",
#     x = "Position Relative to CTCF Binding Site (bp)",
#     y = "Recurrence Rate"
#   ) +
#   theme_minimal()


#############################################
## PLOT CTCFbs vs surroundinig region mutation density for all samples; not position specific
##############################################

##write.csv(data,"mutation_density_surrunding-vs-ctcfbs-all-samples-summary.csv")
#File is already save load for skipping the computation and go to end of the script
 

# Step 1: Read CTCF region coordinates (e.g., in BED format)
# Example: CTCF regions are in a BED file (chrom, start, end)
ctcf_regions <- read.table("motifs_with_loops_sorted_hg38_1806.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(ctcf_regions) <- c("chrom", "start", "end","Strand","motif_type","loop_number")



# Convert CTCF regions into a GenomicRanges object
ctcf_gr <- GRanges(seqnames = ctcf_regions$chrom,
                   ranges = IRanges(start = ctcf_regions$start, end = ctcf_regions$end))

# Step 2: Define the surrounding region ±1kb (1000 bp) around each CTCF region
surrounding_gr <- GRanges(
  seqnames = seqnames(ctcf_gr),
  ranges = IRanges(
    start = start(ctcf_gr) - 990,  # 1kb upstream of CTCF
    end = end(ctcf_gr) + 990       # 1kb downstream of CTCF
  )
)

# Step 3: List all VCF files in your directory (assuming all VCFs are in the same folder)
vcf_files <- list.files(path = "step1-intersect-vcf-with-loop-motif/", pattern = "\\.vcf$", full.names = TRUE)


# Initialize an empty data frame to store mutation density for all samples
mutation_density_summary <- data.frame()

# Step 4: Loop through each VCF file and calculate mutation density
for (vcf_file in vcf_files) {
  
  # Read the VCF file using vcfR
  vcf_data <- read.vcfR(vcf_file)
  
  # Step 5: Extract the variants from the VCF file
  # Extract the genotype information (GT matrix) from the VCF file
  vcf_df <- extract.gt(vcf_data)
  
  # Extract the variant positions (CHROM, POS) and combine with the genotype info
  vcf_variants <- data.frame(
    CHROM = vcf_data@fix[, "CHROM"],
    POS = as.integer(vcf_data@fix[, "POS"])
  )
  
  #### 
  
  # Step 6: Convert VCF mutations into a GenomicRanges object
  mut_gr <- GRanges(seqnames = vcf_variants$CHROM,
                    ranges = IRanges(start = vcf_variants$POS, end = vcf_variants$POS))
  
  # Step 7: Find mutations in CTCF region
  ctcf_mutations <- findOverlaps(mut_gr, ctcf_gr)
  
  # Step 8: Find mutations in the surrounding ±1kb region
  surrounding_mutations <- findOverlaps(mut_gr, surrounding_gr)
  
  # Step 9: Calculate mutation density (mutations per region length)
  ctcf_density <- length(ctcf_mutations) / sum(width(ctcf_gr))
  surrounding_density <- length(surrounding_mutations) / sum(width(surrounding_gr))
  
  # Step 10: Extract the sample name from the VCF filename (remove the file extension)
  sample_name <- gsub("\\.vcf$", "", basename(vcf_file))
  
  # Step 11: Store the results in the summary data frame
  mutation_density_summary <- rbind(mutation_density_summary, data.frame(
    sample = sample_name,
    ctcf_density = ctcf_density,
    surrounding_density = surrounding_density
  ))
}


library(tidyr)
# Step 12: Summarize mutation densities
mutation_density_summary <- mutation_density_summary %>%
  pivot_longer(cols = c(ctcf_density, surrounding_density), 
               names_to = "region", 
               values_to = "mutation_density")

library(dplyr)

data <- mutation_density_summary
# Transform the data
data <- data %>%
  mutate(
    Sample = gsub("CSCC_", "", gsub("_vs_.*", "", sample)),  # Extract and trim Sample
    Mutation_Density_per_Mb = mutation_density * 1e6        # Convert density to per Mb
  )

getwd()
#write.csv(data,"mutation_density_surrunding-vs-ctcfbs-all-samples-summary.csv")


# Step 13: Visualize the mutation density across all samples
ggplot(data, aes(x = Sample, y = Mutation_Density_per_Mb, fill = region)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    #title = "Mutation Density in CTCFbs Region vs ±1kb Surrounding Region",
    x = "Sample",
    y = "Mutations per Mb"
  ) +
  scale_fill_manual(values = c("ctcf_density" = "#1f77b4", "surrounding_density" = "#ff7f0e")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size =6))  # Rotate x-axis labels for better readability

ggplot(data, aes(x = Sample, y = Mutation_Density_per_Mb, fill = region)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = NULL,  # Remove x-axis label
    y = "Mutations per Mb"
  ) +
  scale_fill_manual(values = c("ctcf_density" = "#1f77b4", "surrounding_density" = "#ff7f0e")) +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.y = element_text(size = 8),  # Adjust y-axis text size
    axis.title.y = element_text(size = 10),  # Adjust y-axis title size
    plot.title = element_blank(),  # Remove the plot title
    legend.position = "top",  # Place legend on top
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 10),  # Adjust legend text size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove panel background
    plot.background = element_blank()  # Remove plot background
  )



 