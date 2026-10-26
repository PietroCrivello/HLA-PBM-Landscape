library(dplyr)
library(readr)

# Import the input table from step 3a
Input <- read_delim(
  file.choose(),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# Extract allotype IDs and manually selected motifs
alleles <- Input$ID
N.motif <- Input$Selected_for_analysis

# Create file list based on IDs
file_list <- paste0(alleles, ".txt")

# Initialize list for peptide sequences
peptides_list <- list()

# Column containing the motif assignment in the MoDec output files
condition_column <- "Motif_1st"

# Read each file and retain sequences assigned to the selected motif
for (i in seq_along(file_list)) {
  
  file <- file_list[i]
  condition_value <- N.motif[i]
  
  df <- read.table(
    file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  extracted_values <- df %>%
    filter(.data[[condition_column]] == condition_value) %>%
    select(1) %>%
    pull()
  
  # Group sequences by the first 10 characters of the ID, adjust depending on the size of the allotype name
  id_prefix <- substr(alleles[i], 1, 10)
  
  if (!id_prefix %in% names(peptides_list)) {
    peptides_list[[id_prefix]] <- extracted_values
  } else {
    peptides_list[[id_prefix]] <- c(
      peptides_list[[id_prefix]],
      extracted_values
    )
  }
}

# Pad lists to equal length
max_length <- max(lengths(peptides_list))

peptides_list <- lapply(peptides_list, function(x) {
  length(x) <- max_length
  x
})

# Convert to a data frame and replace missing values with empty strings
peptides_df <- as.data.frame(peptides_list)
peptides_df[is.na(peptides_df)] <- ""

# Save the merged peptide table
print(peptides_df)

write.table(
  peptides_df,
  "peptide_lists_merged.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)