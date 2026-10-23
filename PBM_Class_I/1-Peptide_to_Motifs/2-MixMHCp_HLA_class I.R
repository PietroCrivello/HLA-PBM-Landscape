# Load required libraries
library(future)
library(furrr)
library(readr)
library(tidyverse)  
library(stringr)    

# Define input and output folders
input_folder <- "path to the folder containing the input file"
output_folder <- "path to the folder where the ouput files will be saved"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Function to convert CSV to TXT
prepare_input_txt <- function(csv_file) {
  peptides <- read_csv(csv_file, col_names = FALSE, show_col_types = FALSE) %>%
    pull(1) %>%
    str_trim() %>%
    .[nchar(.) >= 8]
  
  if (length(peptides) < 50) {
    warning(sprintf("⚠️ Skipped %s: fewer than 50 valid 9-mer peptides", basename(csv_file)))
    return(NULL)
  }
  
  # Print input peptides
  cat(sprintf("\n📄 Input peptides from %s:\n", basename(csv_file)))
  print(peptides)
  
  txt_file <- file.path(tempdir(), paste0(tools::file_path_sans_ext(basename(csv_file)), ".txt"))
  write_lines(peptides, txt_file)
  return(txt_file)
}

# Function to run MixMHCp on a given CSV file
run_mixMHCp <- function(csv_file) {
  txt_input <- prepare_input_txt(csv_file)
  if (is.null(txt_input)) return(NULL)
  
  base_name <- tools::file_path_sans_ext(basename(csv_file))
  output_file <- file.path(output_folder, paste0(base_name, "_output"))
  dir.create(output_file, showWarnings = FALSE, recursive = TRUE)
  
  # Build MixMHCp command
  command <- sprintf(
    'path_to_MixMHCp_software --all-lengths -i %s -o %s -m 5 -a ACDEFGHIKLMNPQRSTVWYsty -b path_to_lib_bias_phospho.txt',
    shQuote(txt_input), shQuote(output_file)
  )
  
  # Print the command to be executed
  cat(sprintf("🔧 Running MixMHCp for %s:\n%s\n", base_name, command))
  
  # Run command and capture output
  tryCatch({
    result <- system(command, intern = TRUE)
    cat(sprintf("✅ Finished: %s (%d peptides)\n", base_name, length(read_lines(txt_input))))
    return(result)
  }, error = function(e) {
    cat(sprintf("❌ Error processing %s: %s\n", base_name, e$message))
    return(NULL)
  })
}

# List all input CSV files
input_files <- list.files(input_folder, pattern = "_unnested\\.csv$", full.names = TRUE)
if (length(input_files) == 0) stop("❌ No input files found.")

# Set up parallel processing (adjust number of workers if needed)
plan(multisession, workers = parallel::detectCores() - 1)

# Run MixMHCp in parallel
results <- future_map(input_files, run_mixMHCp)

# Save all output from MixMHCp to a log file
log_file <- file.path(output_folder, "mixmhcp_log.txt")
write_lines(unlist(results), log_file)
