# ====================================================
# Metabolite Association Heatmap Visualization
# Author: [liqin hu]
# Date: [Date]
# Description: Heatmap visualization of metabolite associations
#              with clinical outcomes, colored by metabolite superclass
# ====================================================

# ------------------------------
# 1. LOAD REQUIRED LIBRARIES
# ------------------------------
library(ggplot2)    # Data visualization
library(reshape2)   # Data reshaping
library(dplyr)      # Data manipulation
library(readr)      # Fast data reading
library(patchwork)  # Plot composition

# ------------------------------
# 2. DATA LOADING AND PREPARATION
# ------------------------------

# Load dataset - choose one based on your analysis
# Option 1: Monozygotic twins analysis results
# df <- read_csv("path.csv")


# ------------------------------
# 3. DATA PREPROCESSING
# ------------------------------

# Sort data: first by superclass, then alphabetically by metabolite name
df <- df %>%
  arrange(superclass, Metabolites) %>%
  mutate(
    # Create composite sorting variable
    sort_var = paste(superclass, Metabolites),
    # Convert Metabolites to factor with desired order
    Metabolites = factor(Metabolites, levels = rev(unique(Metabolites)))
  )

# ------------------------------
# 4. DEFINE COLOR SCHEMES
# ------------------------------

# Assign distinct colors for each superclass
superclass_colors <- c(
  "Lipid" = "#E31A1C",              # Red
  "Amino Acid" = "#33A02C",         # Green
  "Xenobiotics" = "#1F78B4",        # Blue
  "Carbohydrate" = "#FF7F00",       # Orange
  "Cofactors and vitamins" = "#FEE090",  # Light yellow
  "Nucleotide" = "#E0F3F8",         # Light blue
  "Energy" = "#FFB6C1"              # Light pink
)

# Keep only superclasses that exist in the data
existing_superclasses <- unique(df$superclass)
superclass_colors <- superclass_colors[names(superclass_colors) %in% existing_superclasses]

# ------------------------------
# 5. CREATE MAIN HEATMAP
# ------------------------------

# Main heatmap using -log10(p-value) for color scale (more scientific)
p_main <- ggplot(df, aes(x = reorder(outcomes, order), 
                         y = Metabolites, 
                         fill = -log10(p_value))) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradientn(
    colors = c("#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#2171B5"),
    name = expression(-log[10](FDR)),  # Scientific notation
    guide = guide_colorbar(
      barwidth = 10,
      barheight = 0.8,
      title.position = "top",
      title.hjust = 0.5
    ),
    na.value = "grey90"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 7),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "top",
    legend.title = element_text(size = 10),
    panel.grid = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "cm")
  ) +
  # Add significance stars/text
  geom_text(
    aes(label = Significance), 
    color = "black", 
    size = 3.5,
    fontface = "bold"
  ) +
  xlab("Clinical Outcomes") +
  ylab("Metabolites")

# ------------------------------
# 6. CREATE SUPERCLASS COLOR BAR
# ------------------------------

# Side color bar showing metabolite superclass
p_side <- ggplot(df, aes(x = 1, y = Metabolites, fill = superclass)) +
  geom_tile() +
  scale_fill_manual(values = superclass_colors, name = "Superclass") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0.2, "cm")
  ) +
  scale_x_discrete(expand = c(0, 0))

# ------------------------------
# 7. COMBINE PLOTS USING PATCHWORK
# ------------------------------

# Combine plots: left side shows superclass color bar, right side shows main heatmap
combined_plot <- p_side + p_main + 
  plot_layout(widths = c(0.03, 1))  # Left bar: 3% width, Main heatmap: 97% width

# Display the combined plot
print(combined_plot)

# ------------------------------
# 8. SAVE PLOT TO FILE
# ------------------------------

# Save as high-resolution image
ggsave(
  filename = "metabolite_association_heatmap.png",
  plot = combined_plot,
  width = 16,          # Adjust width based on number of outcomes
  height = 20,         # Adjust height based on number of metabolites
  dpi = 300,
  bg = "white"
)

# Save as PDF for publication quality
ggsave(
  filename = "metabolite_association_heatmap.pdf",
  plot = combined_plot,
  width = 16,
  height = 20,
  device = "pdf"
)

# ------------------------------
# 9. ADDITIONAL UTILITY FUNCTIONS
# ------------------------------

# Function to create heatmap with custom parameters
create_metabolite_heatmap <- function(data_path, output_prefix = "heatmap") {
  # Load data
  df <- read_csv(data_path)
  
  # Process data
  df <- df %>%
    arrange(superclass, Metabolites) %>%
    mutate(
      Metabolites = factor(Metabolites, levels = rev(unique(Metabolites)))
    )
  
  # Create plot (similar to above)
  p_main <- ggplot(df, aes(x = reorder(outcomes, order), 
                           y = Metabolites, 
                           fill = -log10(p_value))) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradientn(
      colors = c("#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#2171B5"),
      name = expression(-log[10](FDR))
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6),
      panel.grid = element_blank()
    ) +
    geom_text(aes(label = Significance), size = 2.5)
  
  return(p_main)
}

# ====================================================
# END OF SCRIPT
# ====================================================