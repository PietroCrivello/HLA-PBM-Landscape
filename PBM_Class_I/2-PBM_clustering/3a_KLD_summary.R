# Load required libraries
library(fs)
library(readr)
library(dplyr)
library(tibble)
library(readr)

# Select any file within the input folder to get the folder path
Input_folder <- dirname(file.choose())

# List all directories (output folders) inside the input folder
output_folders <- dir_ls(Input_folder, type = "directory", recurse = FALSE)

# Extract and clean folder names
folder_names <- path_file(output_folders)
clean_names <- gsub("_unnested_output", "", folder_names)

# Initialize list to collect KLD values
kld_columns <- list()

for (i in seq_along(output_folders)) {
  # Path to best_ncl.txt
  kld_file <- path(output_folders[i], "KLD", "best_ncl.txt")
  
  if (file_exists(kld_file)) {
    # Read and extract 2nd column
    kld_data <- read_delim(kld_file, delim = "\t", skip = 1, col_names = FALSE, show_col_types = FALSE)
    kld_column <- kld_data[[2]]
    kld_columns[[clean_names[i]]] <- kld_column
  } else {
    cat("KLD file not found in folder:", folder_names[i], "\n")
    kld_columns[[clean_names[i]]] <- rep(NA, 5)
  }
}

# Create data frame
kld_results <- as.data.frame(kld_columns)
rownames(kld_results) <- 1:5

# Transpose the data frame
kld_transposed <- as.data.frame(t(kld_results))

# Add "Max" column: row index of the maximum value
kld_transposed$Max <- apply(kld_transposed, 1, function(x) which.max(x))

# Add "PWM_IDs" column
kld_transposed$PWM_IDs <- sapply(kld_transposed$Max, function(max_val) {
  if (is.na(max_val)) return(NA)
  paste0("PWM_", max_val, "_", 1:max_val, collapse = ";")
})
kld_transposed <- tibble::rownames_to_column(kld_transposed, var = "Allele")

# Convert PWM_IDs column into multiple rows
expanded_df <- kld_transposed %>%
  mutate(PWM_IDs = strsplit(as.character(PWM_IDs), ";")) %>%  # split strings
  tidyr::unnest(PWM_IDs)

# Save the  data
output_file_expanded <- path(Input_folder, paste0(gsub("_outputs", "", basename(Input_folder)), "_KLD_summary.txt"))
write_delim(expanded_df, output_file_expanded, delim = "\t")

cat("Expanded KLD summary saved to:", output_file_expanded, "\n")