# Title: "CTCF loop genes extraction both inside and outside loop upto 1000bp either side with ensemble-id"
# Author: "Amarinder Thind"

library(biomaRt)
library(parallel)
library(dplyr)


setwd('Path')

## Load CTCFbs motif bed file with 3' and 5' Anchor 

motif_1806_hg38 <- read.delim("motifs_with_loops_sorted_hg38_1806.bed", header=FALSE)

# In here, the first side is end of the motif and second side is starting of the motif ## since we captured inside part only 
# Add additional 1000 bp on left side of first motif (of Loop XX) and 1000 bp on right side of second motif ( of Loop Xx)

# separate the 5 prime-motif site from 3 prime-motif (will end up in 903/1806) entries for each

CTCF_5site <- motif_1806_hg38[motif_1806_hg38$V5 =="5_CTCF",] 
CTCF_3site <- motif_1806_hg38[motif_1806_hg38$V5 =="3_CTCF",]

CTCF_5site$start_1000 <- CTCF_5site$V2 -1000
CTCF_5site$end_1000 <- CTCF_5site$V2
 
CTCF_3site$start_1000 <- CTCF_3site$V3 
CTCF_3site$end_1000 <- CTCF_3site$V3 + 1000
 
# Merge the two data frames based on column V6
merged_df <- merge(CTCF_5site[, c("V1",   "V6","V5", "start_1000")], #take start from 5site
                   CTCF_3site[, c("V1",   "V6", "V5", "end_1000")],  ##take end from 3site
                   by = c("V6"))

# View the merged data frame
head(merged_df)

# Connect to the Ensembl database
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)

## confirm loaded reference is hg38 or not?
listDatasets(mart)
listAttributes(mart)

# Define the attributes and filters
attributes <- c("chromosome_name", "start_position", "end_position", "strand","ensembl_gene_id", "hgnc_symbol")
filters <- c("chromosome_name", "start", "end")

# List of genomic regions
regions_5 <- merged_df[,c(2,4,7,1)]

colnames(regions_5) <- c("chromosome_name", "start", "end","identifier")
# Remove "chr" prefix from the chromosome_name column
regions_5$chromosome_name <- gsub("^chr", "", regions_5$chromosome_name)
# Convert dataframe to list of lists
regions_list5 <- apply(regions_5, 1, function(row) {
  list(
    chromosome_name = as.character(row["chromosome_name"]),
    start = as.numeric(row["start"]),
    end = as.numeric(row["end"]),
    identifier = as.character(row["identifier"])
  )
})

#write.csv(regions_5,"modified-Co-ordinates-loop_with-1000k.csv")
 
# Function to query BioMart and add identifiers

get_data_with_identifier <- function(region, mart, attributes, filters) {
  values <- list(
    chromosome_name = region$chromosome_name,
    start = region$start,
    end = region$end
  )
  
  # Query BioMart
  results <- tryCatch({
    getBM(attributes = attributes, filters = filters, values = values, mart = mart)
  }, error = function(e) {
    message(paste("Error in query for identifier:", region$identifier, e$message))
    return(NULL)
  })
  
  # Add identifier to results
  if (!is.null(results) && nrow(results) > 0) {
    results$identifier <- region$identifier
  } else {
    # Create an empty data frame with the same structure as `results` and add identifier
    results <- data.frame(matrix(ncol = length(attributes) + 1, nrow = 1))
    colnames(results) <- c(attributes, "identifier")
    results[] <- NA # Fill with NA
    results$identifier <- region$identifier
  }
  
  return(results)
}

# Parallelize BioMart queries
parallel_get_data <- function(regions_list, mart, attributes, filters) {
  # Use parallel processing (you can adjust the number of cores depending on your system)
  num_cores <- detectCores() - 1 # Use one less than the number of cores for safety
  results_list <- mclapply(regions_list, get_data_with_identifier, mart = mart, attributes = attributes, filters = filters, mc.cores = num_cores)
  
  # Combine results into one data frame
  all.genes <- do.call(rbind, results_list)
  return(all.genes)
}

# Query BioMart for each region in parallel and combine results
all.genes_5 <- parallel_get_data(regions_list5, mart, attributes, filters)

# Group by 'identifier' and summarize other columns
aggregated_data_5 <- all.genes_5 %>%
  group_by(identifier) %>%
  summarize(
    chromosome_name = first(chromosome_name),
    start_position = first(start_position),
    end_position = first(end_position),
    strand = first(strand),
    ensembl_gene_id = paste(unique(ensembl_gene_id), collapse = ", "),
    hgnc_symbol = paste(unique(hgnc_symbol), collapse = ", "),
    .groups = 'drop'
  )

# Clean function for NA values

# Clean function
clean_column <- function(column) {
  # Replace standalone "NA" with empty strings and trim whitespace
  cleaned_column <- gsub("\\bNA\\b", "", column)
  
  # Remove leading/trailing commas and extra commas
  cleaned_column <- gsub("^,\\s*|,\\s*$", "", cleaned_column) # Remove leading/trailing commas
  cleaned_column <- gsub(",\\s*,", ", ", cleaned_column) # Remove extra commas
  cleaned_column <- gsub(",\\s*$", "", cleaned_column) # Remove trailing comma
  
  # Remove leading comma if it's the only entry
  cleaned_column <- gsub("^\\s*,", "", cleaned_column) # Remove leading comma
  cleaned_column <- gsub(",\\s*,$", "", cleaned_column) # Remove trailing comma if it's the only entry
  
  return(cleaned_column)
}


# Apply cleaning function to relevant columns
df_cleaned <- aggregated_data_5 %>%
  mutate(
    ensembl_gene_id = clean_column(ensembl_gene_id),
    hgnc_symbol = clean_column(hgnc_symbol)
  )

# ### We lost few chromesom/start/end entries inthe beiging of region file, We can fill this from previous one
# ## Ideally, we can skip these cols from most of the cleaning
# 
# colnames(df_cleaned)
# colnames(regions_5)
# 
# df_cleaned2 <- full_join( regions_5, df_cleaned[,c("ensembl_gene_id","hgnc_symbol","identifier") ], by = "identifier", suffix = c(".df1", ".df2"))

write.csv(df_cleaned,"Loops_with_gene_list_combine-inside-and-1000bp-outside-.csv")





