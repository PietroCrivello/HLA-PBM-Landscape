# Load required library
library(fs)
library(readr)

# Load the input file (KLD summary)
Input_file <- file.path(file.choose())
Input <- read_delim(Input_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
Input_folder <- dirname(Input_file)

# Define the list of folders from the Input file
folders <- paste0(c(Input$Allele), "_unnested_output")
PWM_ID <- c(Input$PWM_IDs)
Allele <- c(Input$Allele)

# Construct the full path for each folder
folders.name <- path(Input_folder, folders, "Multiple_PWMs")

# Define the output folder
output_folder <- path(getwd(), gsub("outputs", "PWM", basename(Input_folder)))

# Define the specific file names to be copied
files_to_copy <- paste0(PWM_ID, ".txt")

# Create the output folder if it doesn't exist
if (!dir_exists(output_folder)) {
  dir_create(output_folder)
}

# Copy one specified file from each folder to the output folder using the folder names as file names
for (i in seq_along(folders.name)) {
  # Create the full path to the specific file in the folder
  file_path <- path(folders.name[i], files_to_copy[i])
  
  # Check if the file exists in the folder
  if (file_exists(file_path)) {
    # Use the folder name as the new file name
    new_file_name <- paste0(Allele[i], ".txt")
    
    # Copy the file to the output folder with the new name
    file_copy(file_path, path(output_folder, new_file_name), overwrite = TRUE)
    cat("Copied:", files_to_copy[i], "from", folders.name[i], "to", new_file_name, "in", output_folder, "\n")
  } else {
    cat("File not found:", files_to_copy[i], "in folder:", folders.name[i], "\n")
  }
}
