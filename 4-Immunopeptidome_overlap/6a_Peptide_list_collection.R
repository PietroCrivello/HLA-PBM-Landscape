# Load required library
library(fs)
library(readr)

# Load the input file
Input <- read_delim(file.choose(), delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Define the list of folders from the Input file
folders <- c(Input$`Folder_output`)

# Construct the full path for each folder
folders.name <- paste0(folders, "_unnested_output.txt\\Responsibilities")

# Define the output folder
output_folder <- "path_to_MoDec_output_fodlers"

# Get the corresponding file information from the Input file
file.K <- c(Input$`N_of_motifs`)
file.selected_motif <- c(Input$`Selected_for_analysis`)
ID <- c(Input$ID)

# Define the specific file names to be copied
files_to_copy <- paste0("bestPepResp_K", file.K, ".txt")

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
    new_file_name <- paste0(ID[i], ".txt")
    
    # Copy the file to the output folder with the new name
    file_copy(file_path, path(output_folder, new_file_name), overwrite = TRUE)
    cat("Copied:", files_to_copy[i], "from", folders.name[i], "to", new_file_name, "in", output_folder, "\n")
  } else {
    cat("File not found:", files_to_copy[i], "in folder:", folders.name[i], "\n")
  }
}

cat("Specified files copied and renamed successfully to", output_folder, "\n")
