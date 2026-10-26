# Load necessary libraries
library(readr)    # For reading the data
library(Rtsne)    # For t-SNE computation
library(ggplot2)  # For plotting
library(RColorBrewer)
library(ggrepel)
library(scales)
library(svglite)

#initial setup (choose PBM for manually assigned PBM groups or Cluster for computed clusters from step 4)

Plot_type <- "Cluster"

custom_colors <- c(
  "#999999",  # Grey
  "#2A5B88",  # Dark Blue
  "#317132",  # Dark Green
  "#B35900",  # Dark Orange
  "#6B3368",  # Dark Purple
  "#B2B229",  # Dark Yellow
  "#6BAED6",  # instead of "#377EB8",  # Blue
  "#E41A1C",  # Red
  "#4D8372",  # Dark Teal
  "#B35745",   # Dark Salmon
  "#A31414"  # Dark Red
)

# Load the matrix data from the .txt file (assuming it's tab-delimited)
print("input data matrix")
file_path <- file.choose()  # This will open a file dialog to choose the .txt file
data <- read_delim(file_path, delim = "\t", col_names = TRUE)

#input info
print("input clustering info")
Input_info <- read_delim(file.choose(), delim = "\t", escape_double = FALSE, trim_ws = TRUE)


# Remove the first 2 columns if they are non-numeric identifiers (e.g., "AA" and "Position")
data <- data[, -c(1, 2)]

data_matrix<-t(data)
rownames(data_matrix)<- colnames(data)

# Store row names (sample names) before converting to matrix
sample_names <- rownames(data_matrix)

# define groups based on hierarchical clustering or PBM
Group_PBM_data <- as.character(c(Input_info$PBM))
Cluster_data <- as.character(c(Input_info$Cluster))

# Convert to matrix and ensure numeric values
data_matrix <- as.matrix(data_matrix)
data_matrix <- apply(data_matrix, 2, as.numeric)

# Reassign row names (sample names) to the matrix
rownames(data_matrix) <- sample_names

# Perform t-SNE dimensionality reduction
set.seed(20)  # Set seed for reproducibility
tsne_result <- Rtsne(data_matrix, dims = 2, perplexity = 5, verbose = TRUE, max_iter = 100000)

# Create a data frame with t-SNE results
tsne_data <- as.data.frame(tsne_result$Y)
colnames(tsne_data) <- c("tSNE1", "tSNE2")
tsne_data$Sample <- sample_names
tsne_data$Group_PBM <- Group_PBM_data
tsne_data$Group_cluster <- Cluster_data


if (Plot_type == "PBM") {
  selected_group <- Group_PBM_data
} else if (Plot_type == "Cluster") {
  selected_group <- Cluster_data
}

# Plot t-SNE using ggplot2
tsne_plot <- ggplot(tsne_data, aes(x = tSNE1, y = tSNE2, label = Sample, color = selected_group)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2.5, family = "Arial", color = "black") +
  labs(x = "t-SNE Dimension 1", y = "t-SNE Dimension 2", title = NULL) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(labels = scientific) +
  scale_y_continuous(labels = scientific) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Display the plot
print(tsne_plot)

# Save as JPEG
jpeg_file <- paste0("tsne_plot_by_sample_", Plot_type, ".jpeg")
ggsave(
  filename = jpeg_file,
  plot = tsne_plot,
  device = "jpeg",
  width = 5,
  height = 4,
  units = "in",
  dpi = 300
)

# Save as SVG (vector graphics)
svg_file <- paste0("tsne_plot_by_sample_", Plot_type, ".svg")
ggsave(
  filename = svg_file,
  plot = tsne_plot,
  device = svglite::svglite,
  width = 5,
  height = 4,
  units = "in"
)
