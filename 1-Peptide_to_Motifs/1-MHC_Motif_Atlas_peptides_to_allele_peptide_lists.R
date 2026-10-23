# Load the necessary libraries
library(readr)
library(tidyverse)

# Input of data from MHC Motif Atlas (http://mhcmotifatlas.org/class1)
ligands_input <- read_delim(file.choose(), 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

ligands_input$length <- NULL

# Remove duplicate rows based on all columns
ligands_input_unique <- ligands_input[!duplicated(ligands_input), ]

peptide_list <- ligands_input_unique %>%
  pivot_wider(names_from = Allele, values_from = Peptide, )

peptide_list <- peptide_list %>% select(-Locus)

df_list <- peptide_list %>%
  map(~ tibble(value = .x) %>% filter(!is.na(value))) %>%
  set_names(names(peptide_list))

# Replace non-alphanumeric characters with underscores in file names
sanitize_filename <- function(name) {
  name <- gsub("[^a-zA-Z0-9]", "_", name)
  return(name)
}

# Unnest each dataframe and save to CSV files
walk2(df_list, names(df_list), ~ {
  unnested_df <- .x %>% unnest(cols = value)
  
  # Remove column names and save to CSV file
  file_name <- paste0(sanitize_filename(.y), "_unnested.csv")
  write_csv(unnested_df, file_name, col_names = FALSE)  # Set col_names to FALSE to omit column names
})

