# Load necessary packages
library(pROC)
library(ggplot2)
library(caret)
library(dplyr)
library(grid)  # Add grid package
library(gridExtra)

# 6-month gross motor as example
selected_metabs <- c(
  "Com_35341_pos",
  "Com_8451_pos",               
  "Com_550_pos",
  "Com_22881_pos",
  "Com_27317_pos",
  "Com_57_neg",
  "Com_45_neg"
)

# 1. Data preparation ----------------------------------------------------------
data_final1 <- data_final[, c(selected_metabs, "Group")]
data_final1$Group <- factor(data_final1$Group)

# 2. Set cross-validation parameters -------------------------------------------------
set.seed(123)
ctrl <- trainControl(method = "repeatedcv", 
                     number = 10,
                     repeats = 5,
                     summaryFunction = twoClassSummary,
                     classProbs = TRUE,
                     savePredictions = "final")

# 3. Calculate single metabolite AUC and 95% CI --------------------------------------
calculate_auc_with_ci <- function(metab, data) {
  temp_data <- data.frame(Y = data$Group, X = data[[metab]])
  
  cv_model <- train(Y ~ X,
                    data = temp_data,
                    method = "glm",
                    family = "binomial",
                    trControl = ctrl,
                    metric = "ROC")
  
  # Calculate 95% CI
  roc_obj <- roc(response = cv_model$pred$obs,
                 predictor = cv_model$pred$Group2,
                 levels = c("Group1", "Group2"))
  
  ci_obj <- ci.auc(roc_obj)
  
  return(data.frame(
    Metabolite = metab,
    AUC = ci_obj[2],
    AUC_lower = ci_obj[1],
    AUC_upper = ci_obj[3],
    AUC_SD = sd(cv_model$results$ROC)
  ))
}

# Calculate AUC and CI for all metabolites
single_auc_list <- lapply(selected_metabs, calculate_auc_with_ci, data = data_final1)
single_auc_df <- do.call(rbind, single_auc_list)
single_auc_df <- single_auc_df[order(-single_auc_df$AUC), ]

cat("\n=== Single Metabolite AUC (95% CI) ===\n")
print(single_auc_df)

# 4. Calculate combined model AUC and 95% CI -------------------------------------------
joint_data <- data_final[, c(selected_metabs, "Group")]
joint_model <- train(Group ~ .,
                     data = joint_data,
                     method = "glm",
                     family = "binomial",
                     trControl = ctrl,
                     metric = "ROC")

joint_roc <- roc(response = joint_model$pred$obs,
                 predictor = joint_model$pred$Group2,
                 levels = c("Group1", "Group2"))
joint_ci <- ci.auc(joint_roc)

# 5. Calculate score model AUC and 95% CI -------------------------------------------
score_data <- data_final[, c("Score", "Group")]
score_model <- train(Group ~ .,
                     data = score_data,
                     method = "glm",
                     family = "binomial",
                     trControl = ctrl,
                     metric = "ROC")

score_roc <- roc(response = score_model$pred$obs,
                 predictor = score_model$pred$Group2,
                 levels = c("Group1", "Group2"))
score_ci <- ci.auc(score_roc)

# 6. Create ROC curve list and CI data frame ----------------------------------------
roc_list <- list()
ci_df <- data.frame()

# Add combined model
roc_list[["Combined Model"]] <- joint_roc
ci_df <- rbind(ci_df, data.frame(
  Model = "Combined Model",
  AUC = joint_ci[2],
  AUC_lower = joint_ci[1],
  AUC_upper = joint_ci[3]
))

# Add score model
roc_list[["Metabolite Score"]] <- score_roc
ci_df <- rbind(ci_df, data.frame(
  Model = "Metabolite Score",
  AUC = score_ci[2],
  AUC_lower = score_ci[1],
  AUC_upper = score_ci[3]
))

# Add single metabolites
for(i in 1:nrow(single_auc_df)) {
  metab <- single_auc_df$Metabolite[i]
  temp_data <- data.frame(Y = data_final$Group, X = data_final[[metab]])
  temp_model <- train(Y ~ X, data = temp_data, method = "glm", 
                      trControl = trainControl(method = "none"))
  
  prob <- predict(temp_model, type = "prob")[[2]]
  roc_obj <- roc(response = temp_data$Y, predictor = prob, 
                 levels = c("Group1", "Group2"))
  
  roc_list[[metab]] <- roc_obj
  ci_df <- rbind(ci_df, data.frame(
    Model = metab,
    AUC = single_auc_df$AUC[i],
    AUC_lower = single_auc_df$AUC_lower[i],
    AUC_upper = single_auc_df$AUC_upper[i]
  ))
}

# 7. Generate plot data -------------------------------------------------------
plot_data <- data.frame()
for(model_name in names(roc_list)) {
  roc_obj <- roc_list[[model_name]]
  temp_df <- data.frame(
    Model = model_name,
    Sensitivity = roc_obj$sensitivities,
    Specificity = roc_obj$specificities
  )
  plot_data <- rbind(plot_data, temp_df)
}

# Merge AUC and CI information
plot_data <- merge(plot_data, ci_df, by = "Model")

# Metabolite name mapping
name_mapping <- c(
  "Metabolite Score" = "Metabolite Score",
  "Combined Model" = "Combined Model",
  "Com_35341_pos" = "Cytosine",
  "Com_8451_pos" = "DHEA",               
  "Com_550_pos" = "Betaine",
  "Com_22881_pos" = "7-Methylxanthine",
  "Com_27317_pos" = "Cholecalciferol",
  "Com_57_neg" = "LPC 18:2",
  "Com_45_neg" = "Elaidic acid"
)

plot_data$Model_Display <- name_mapping[plot_data$Model]

# Create labels containing AUC and 95% CI
plot_data$Model_Label <- sprintf("%s\nAUC: %.3f (%.3f-%.3f)", 
                                 plot_data$Model_Display,
                                 plot_data$AUC,
                                 plot_data$AUC_lower,
                                 plot_data$AUC_upper)

# Order by AUC
model_order <- unique(plot_data$Model_Label[order(-plot_data$AUC)])
plot_data$Model_Label <- factor(plot_data$Model_Label, levels = model_order)

# 8. Plot ROC curves with 95% CI -------------------------------------------

n_models <- length(unique(plot_data$Model_Label))
color_palette <- scales::hue_pal()(n_models)

box_plot <- ggplot(plot_data, aes(x = 1-Specificity, y = Sensitivity, 
                                  color = Model_Label, group = Model_Label)) +
  geom_line(size = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray", alpha = 0.7) +
  labs(title = "Metabolites for 6-month gross motor prediction",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = "Model with AUC (95% CI)") +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    # Adjust legend style - key modifications here
    legend.text = element_text(size = 8.5, 
                               margin = margin(t = 4, b = 4, unit = "pt")),  # Increase top and bottom margins
    legend.title = element_text(size = 10, face = "bold", 
                                margin = margin(b = 10, unit = "pt")),  # Increase title bottom margin
    legend.spacing.y = unit(8, "pt"),  # Increase vertical spacing between legend items
    legend.key.height = unit(16, "pt"),  # Increase legend key height
    legend.key = element_rect(fill = "white", color = "white"),  # White background
    legend.background = element_rect(fill = "white", color = NA)  # White legend background
  ) +
  scale_color_manual(values = color_palette,
                     labels = function(x) gsub("\\n", "\n", x)) +  # Ensure line breaks are parsed correctly
  coord_equal() +
  guides(color = guide_legend(
    override.aes = list(size = 1.5),
    ncol = 1,
    byrow = TRUE,
    keyheight = unit(18, "pt")  # Further increase key height
  ))

print(box_plot)






