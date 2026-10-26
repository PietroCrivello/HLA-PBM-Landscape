# Load necessary libraries
library(future)
library(furrr)
library(readr)
library(tidyverse)

# Define the function to run MODEC
run_modec <- function(input_file) {
  output_file <- file.path(output_folder, paste0(tools::file_path_sans_ext(basename(input_file)), "_output.txt"))
  
  command <- sprintf(
    '"MoDec.exe" --Kmax 5 --MHC2 --nruns 10 --makeReport --outAdd -i %s -o %s',
    shQuote(input_file),
    shQuote(output_file)
  )
  
  # Run the command and capture errors
  tryCatch(
    {
      system(command, intern = TRUE)
      cat(sprintf("Processing of %s completed successfully.\n", input_file))
    },
    error = function(e) {
      cat(sprintf("Error processing %s: %s\n", input_file, e$message))
    }
  )
}
# Update this path to the new input folder
input_folder <- "path_to_input_peptide_lists"
output_folder <- "path_to_output_folder"

# List of input files
input_files <- list.files(
  input_folder,
  pattern = "\\unnested.csv$",  # Ensure proper escaping of special characters
  full.names = TRUE
)

# Ensure input files exist
if (length(input_files) == 0) {
  stop("No input files found matching the pattern '\\unnested.csv$'")
}

# Set up parallel processing
plan(multisession, workers = 6) # or alternative parallel::detectCores() - 1

# Run MODEC in parallel
future_map(input_files, run_modec)
