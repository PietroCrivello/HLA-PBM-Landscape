# Load necessary libraries
library(dplyr)
library(readr)

# Import list of alleles
List.alleles <- read_delim(file.choose(), delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Extract alleles and corresponding motif values
alleles <- List.alleles$`PWM`

# Create file list based on allele names
file_list <- paste0(alleles, ".txt")

# Initialize list for peptides
peptides_list <- list()

# Read and filter peptide data for each file
for (i in seq_along(file_list)) {
  
  file <- file_list[i]
  
  # Read the file
  df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Filter and select peptides/core column
  extracted_values <- df %>%
    filter(!!sym("X1") > 0.5) %>%
    filter(!grepl("[sty]", .[[1]])) %>%  # To exclude peptides containing lowercase s, t, or y
    pull(1)
  
  # Use ID prefix (e.g., first 2 characters) as group key
  id_prefix <- substr(alleles[i], 1, 19)
  
  # Append values to the appropriate group in the list
  if (!id_prefix %in% names(peptides_list)) {
    peptides_list[[id_prefix]] <- extracted_values
  } else {
    peptides_list[[id_prefix]] <- c(peptides_list[[id_prefix]], extracted_values)
  }
}

# Find max vector length for padding
max_length <- max(sapply(peptides_list, length))

# Pad vectors to equal length with NAs
peptides_list <- lapply(peptides_list, function(x) {
  length(x) <- max_length
  return(x)
})

# Convert to data frame
peptides_df <- as.data.frame(peptides_list)

# Replace NAs with empty strings
peptides_df[is.na(peptides_df)] <- ""

# Print and optionally save
print(peptides_df)
write.table(peptides_df, "peptides_lists_merged.txt", sep = "\t", row.names = FALSE, quote = FALSE)
