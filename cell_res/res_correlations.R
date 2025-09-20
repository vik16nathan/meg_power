# Load required libraries
library(R.matlab)
library(ggplot2)
library(ggpubr)


#find the average across n=4 subjects
for(subject in c("sub-0002", "sub-0004", "sub-0006", "sub-0007")) {
  # Extract variables
  vertex_res_map_old <- as.numeric(mat_data[[1]])
  vertex_resistance_map <- as.numeric(mat_data[[2]])
  mat_data <- readMat(paste0(subject, "_region_res_cell_column.mat"))
  # Create a data frame
  data <- data.frame(vertex_res_map_old = vertex_res_map_old, vertex_resistance_map = vertex_resistance_map)
  
  #filter zeroes
  data <- data[which(data$vertex_resistance_map != 0),]
  
  subavg_vert_res_old <- c(subavg_vert_res_old, mean(data$vertex_res_map_old))
  subavg_vert_res_new <- c(subavg_vert_res_new, mean(data$vertex_resistance_map))
  
}

###summary stats
mean(subavg_vert_res_new)
min(subavg_vert_res_new)
max(subavg_vert_res_new)

mean(subavg_vert_res_old)
min(subavg_vert_res_old)
max(subavg_vert_res_old)

###repeat all for sub-0002 for visualization purposes
# Load the .mat file
setwd("C:/Users/vik16/Documents/Baillet Lab Manuscript/")
mat_data <- readMat("sub-0002_region_res_cell_column.mat")

subavg_vert_res_old <- c()
subavg_vert_res_new <- c()
# Extract variables
vertex_res_map_old <- as.numeric(mat_data[[1]])
vertex_resistance_map <- as.numeric(mat_data[[2]])

# Create a data frame
data <- data.frame(vertex_res_map_old = vertex_res_map_old, vertex_resistance_map = vertex_resistance_map)

# Filter zeros
data <- data[which(data$vertex_resistance_map != 0),]

# Correlation with full dataset
correlation <- cor(data$vertex_res_map_old, data$vertex_resistance_map, method = "spearman")

# Find top-right outlier (for display only)
nrml <- function(x) (x - min(x)) / (max(x) - min(x))
score <- nrml(data$vertex_res_map_old) + nrml(data$vertex_resistance_map)
idx_outlier <- which.max(score)

display_data <- data[-idx_outlier, ]
all_data     <- data

# Line-of-unity segment
xr <- range(display_data$vertex_res_map_old, na.rm = TRUE)
yr <- range(display_data$vertex_resistance_map, na.rm = TRUE)
x_min <- max(min(xr), min(yr))
x_max <- min(max(xr), max(yr))
line_df <- data.frame(x = c(x_min, x_max), y = c(x_min, x_max))

p <- ggplot(display_data, aes(x = vertex_res_map_old, y = vertex_resistance_map)) +
  geom_point(show.legend = FALSE) +
  # Line of unity with single legend entry
  geom_line(data = line_df,
            aes(x = x, y = y, linetype = "Line of unity"),
            color = "blue", linewidth = 1) +
  scale_linetype_manual(values = c("Line of unity" = "dashed")) +
  labs(
    x = "R, GM resistivity in cortical column (Ohms)",
    y = "R, dendrite and neuron density (Ohms)",
    title = "Dendritic vs. Columnar Resistance (sub-0002)"
  ) +
  annotate(
    "text",
    x = max(display_data$vertex_res_map_old, na.rm = TRUE),
    y = max(display_data$vertex_resistance_map, na.rm = TRUE),
    label = paste("Spearman r =", round(correlation, 2)),
    hjust = 1.1, vjust = 1.1, size = 6, color = "darkred"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.position = "top"
  ) +
  guides(linetype = guide_legend(override.aes = list(color = "blue", linewidth = 1.2)))

print(p)
