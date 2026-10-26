# Load the necessary libraries
library(readr)
library(tidyverse)

# Input dataset
ligands_input <- read_csv(file.choose())
ligands_input <- ligands_input %>% select(-1)

# Convert all columns to character type
ligands_input <- ligands_input %>% mutate(across(everything(), as.character))

# Remove duplicate rows based on all columns
ligands_input_unique <- ligands_input[!duplicated(ligands_input), ]

peptide_list <- ligands_input_unique %>%
  pivot_wider(names_from = allele, values_from = sequence, )

df_list <- peptide_list %>%
  map(~ tibble(value = .x) %>% filter(!is.na(value))) %>%
  set_names(names(peptide_list))

# Function to sanitize file names
sanitize_filename <- function(name) {
  name <- gsub("[^a-zA-Z0-9]", "_", name)  # Replace non-alphanumeric characters with underscores
  return(name)
}

# Unnest each dataframe and save to CSV files
walk2(df_list, names(df_list), ~ {
  unnested_df <- .x %>% unnest(cols = value)  # Unnest the 'value' column
  
  # Remove column names and save to CSV file
  file_name <- paste0(sanitize_filename(.y), "_unnested.csv")
  write_csv(unnested_df, file_name, col_names = FALSE)  # Set col_names to FALSE to omit column names
})

# Combine all dataframes into one dataframe with metadata
combined_df <- bind_rows(df_list, .id = "source_column")

# Rename the `source_column` column to reflect the original column names
combined_df <- combined_df %>%
  rename(column_name = source_column)

# Inspect the combined dataframe
print(combined_df)

# Optionally, save the combined dataframe to a CSV file
write_csv(combined_df, "combined_data_with_source.csv")
