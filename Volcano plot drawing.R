# ====================================================
# Volcano Plot Visualization for Metabolite Associations
# Author: [liqin hu]
# Date: [Date]
# Description: Multi-panel volcano plot visualization for 
#              metabolite association analysis results
# ====================================================

# ------------------------------
# 1. LOAD REQUIRED LIBRARIES
# ------------------------------
library(ggplot2)    # Data visualization
library(dplyr)      # Data manipulation
library(readr)      # Fast data reading
library(ggh4x)      # Enhanced ggplot2 facets

# ------------------------------
# 2. DATA LOADING
# ------------------------------
# Load the association analysis results
# The dataset should contain:
# - Correlation/Coefficient values
# - P-values and FDR-adjusted p-values
# - Group labels for faceting
# - Metabolite labels
data <- read_csv("path.csv")

# ------------------------------
# 3. DATA PREPROCESSING
# ------------------------------

# Determine number of groups for coloring
n_groups <- length(unique(data$Group))

# Create grouping labels with colors
group_labels <- data.frame(
  Group = 1:n_groups,
  Label = as.character(1:n_groups),
  Color = scales::hue_pal()(n_groups)  # Generate distinct colors
)

# Define significance based on FDR-adjusted p-values
data <- data %>%
  mutate(
    Significance = case_when(
      FDR_P_value < 0.05 & Correlation > 0 ~ "Positive",
      FDR_P_value < 0.05 & Correlation < 0 ~ "Negative",
      TRUE ~ "Not Significant"
    )
  )

# Create colored facet labels
facet_labels_colored <- setNames(
  paste0(group_labels$Label),
  1:n_groups
)

# ------------------------------
# 4. VOLCANO PLOT VISUALIZATION
# ------------------------------

volcano_plot <- ggplot(data, aes(x = Correlation, y = -log10(P_value))) +
  # Scatter points colored by significance
  geom_point(
    aes(color = Significance),
    size = 1.5,
    alpha = 0.7
  ) +
  
  # Text labels for significant points
  geom_text(
    data = subset(data, FDR_P_value < 0.05),
    aes(label = Label),
    color = "black",
    size = 2.5,
    vjust = -0.5,
    hjust = 0.5,
    check_overlap = TRUE
  ) +
  
  # Vertical reference line at x = 0
  geom_vline(
    xintercept = 0, 
    linetype = "dashed",
    color = "black", 
    alpha = 0.6,
    linewidth = 0.4
  ) +
  
  # Color scale for significance
  scale_color_manual(
    values = c(
      "Positive" = "red",
      "Negative" = "blue", 
      "Not Significant" = "gray70"
    ),
    name = "Significance",
    guide = guide_legend(position = "bottom")
  ) +
  
  # Faceting with colored strips
  facet_wrap2(
    ~ Group, 
    scales = "free", 
    nrow = 4,
    labeller = as_labeller(facet_labels_colored),
    strip = strip_themed(
      background_x = elem_list_rect(
        fill = group_labels$Color,
        color = group_labels$Color,
        linewidth = 1
      )
    )
  ) +
  
  # Axis scaling
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.1)),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6)
  ) +
  
  # Theme and labels
  theme_minimal() +
  labs(
    x = "Coefficients",
    y = expression(-log[10](italic(P))),  # Scientific notation
    title = "Volcano Plot of Metabolite Associations",
    subtitle = "Multi-panel visualization of association results"
  ) +
  
  # Custom theme settings
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 14, color = "white"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.line = element_line(color = "gray30"),
    axis.ticks = element_line(color = "gray30"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Display the plot
print(volcano_plot)

# ------------------------------
# 5. SAVE PLOT TO FILE
# ------------------------------

# Save as PDF (publication quality)
ggsave(
  "volcano_plot_multi_panel.pdf",
  plot = volcano_plot,
  width = 14,
  height = 10,
  dpi = 300
)


# ====================================================
# END OF SCRIPT
# ====================================================