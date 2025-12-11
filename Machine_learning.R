# Load necessary packages
library(mixOmics)
library(ggplot2)
library(plotly)
library(readr)

# Set random seed for reproducibility
set.seed(123)

data=read_csv("path.csv")


# Simplified metabolomics analysis pipeline (fixed errors and added 95% CI)
library(mixOmics)
library(glmnet)
library(pROC)
library(caret)
library(doParallel)

set.seed(123)

# 1. Data preparation and preprocessing ------------------------------------------------------
data_scaled <- na.omit(data)

# Metabolite column processing
metabolite_cols <- 1:447
data_scaled[, metabolite_cols] <- log2(data_scaled[, metabolite_cols] + 1)
data_scaled[, metabolite_cols] <- scale(data_scaled[, metabolite_cols])

# Categorical variable processing
data_scaled$Group <- factor(data_scaled$Group)
data_scaled$sex <- factor(data_scaled$sex)
data_scaled$zygosity <- factor(data_scaled$zygosity)
data_scaled$gestational_age <- factor(data_scaled$gestational_age)

# 2. PLS-DA initial screening ----------------------------------------------------------
plsda_model <- mixOmics::plsda(X = data_scaled[, metabolite_cols], 
                               Y = data_scaled$Group, ncomp = 2)

vip_scores <- mixOmics::vip(plsda_model)
vip_values <- vip_scores[, 1]
selected_vip <- names(vip_values)[vip_values > 1.5]
cat("PLS-DA selected metabolites:", length(selected_vip), "\n")

# 3. Non-parametric test screening -------------------------------------------------------
p_values <- sapply(selected_vip, function(x) wilcox.test(data_scaled[[x]] ~ data_scaled$Group)$p.value)
significant_metabolites <- names(p_values)[p.adjust(p_values, method = "BH") < 0.05] 
cat("Significant metabolites:", length(significant_metabolites), "\n")

# 4. LASSO stepwise screening -------------------------------------------------------
if(length(significant_metabolites) > 0) {
  
  # Prepare data matrix
  x_matrix <- cbind(
    data_scaled[, significant_metabolites, drop = FALSE],
    gestational_age = as.numeric(data_scaled$gestational_age) - 1,
    sex = as.numeric(data_scaled$sex) - 1,
    zygosity = as.numeric(data_scaled$zygosity) - 1
  )
  x_matrix <- as.matrix(x_matrix)
  y_factor <- data_scaled$Group
  
  covariate_names <- c("gestational_age", "sex", "zygosity")
  
  # Set penalty factor (covariates not penalized)
  penalty_factor <- c(rep(1, length(significant_metabolites)), 0, 0, 0)
  
  # Initial full model
  full_model <- cv.glmnet(x_matrix, y_factor, 
                          family = "binomial", alpha = 1,
                          type.measure = "auc", nfolds = 10,
                          penalty.factor = penalty_factor)
  
  # Stepwise elimination process
  current_metabs <- significant_metabolites
  minimal_metabs <- 2
  
  repeat {
    if (length(current_metabs) <= minimal_metabs) {
      final_metabolites <- current_metabs
      break
    }
    
    # Evaluate contribution of each metabolite
    temp_results <- list()
    for(metab in current_metabs) {
      temp_metabs <- setdiff(current_metabs, metab)
      temp_predictors <- c(temp_metabs, covariate_names)
      temp_penalty <- c(rep(1, length(temp_metabs)), 0, 0, 0)
      
      temp_model <- cv.glmnet(x_matrix[, temp_predictors, drop = FALSE], y_factor,
                              family = "binomial", alpha = 1,
                              type.measure = "auc", nfolds = 10,
                              penalty.factor = temp_penalty)
      
      pred_full <- predict(full_model, x_matrix, s = "lambda.min", type = "response")[,1]
      pred_temp <- predict(temp_model, x_matrix[, temp_predictors, drop = FALSE], 
                           s = "lambda.min", type = "response")[,1]
      
      roc_full <- roc(y_factor, pred_full, quiet = TRUE)
      roc_temp <- roc(y_factor, pred_temp, quiet = TRUE)
      
      delong_test <- roc.test(roc_full, roc_temp, method = "delong")
      
      temp_results[[metab]] <- data.frame(
        metabolite = metab,
        temp_auc = max(temp_model$cvm),
        p_value = delong_test$p.value
      )
    }
    
    # Combine results
    auc_results <- do.call(rbind, temp_results)
    auc_results$adj_p <- p.adjust(auc_results$p_value, method = "fdr")
    
    # Find the least important metabolite
    max_p_row <- which.max(auc_results$adj_p)
    max_p_metab <- auc_results$metabolite[max_p_row]
    max_p <- auc_results$adj_p[max_p_row]
    
    # Termination condition
    if (max_p < 0.05) {
      final_metabolites <- current_metabs
      break
    }
    
    # Remove this metabolite
    current_metabs <- setdiff(current_metabs, max_p_metab)
    cat("Removed:", max_p_metab, "| Adjusted p-value:", round(max_p, 3), 
        "| Remaining:", length(current_metabs), "\n")
  }
  
  cat("LASSO final selected metabolites:", length(final_metabolites), "\n")
  
  # Build final model
  final_predictors <- c(final_metabolites, covariate_names)
  final_penalty <- c(rep(1, length(final_metabolites)), 0, 0, 0)
  
  final_model <- cv.glmnet(x_matrix[, final_predictors, drop = FALSE], y_factor,
                           family = "binomial", alpha = 1,
                           type.measure = "auc", nfolds = 10,
                           penalty.factor = final_penalty)
  
  # Fix coefficient extraction error
  coefs <- coef(final_model, s = "lambda.min")
  coef_values <- as.numeric(coefs)  # Convert to numeric vector
  coef_names <- rownames(coefs)     # Get variable names
  
  # Extract non-zero coefficient metabolites (excluding intercept and covariates)
  non_zero_indices <- which(coef_values != 0)
  non_zero_coefs <- coef_values[non_zero_indices]
  non_zero_names <- coef_names[non_zero_indices]
  
  # Keep only metabolites (excluding intercept and covariates)
  metab_indices <- non_zero_names %in% final_metabolites
  final_coef_values <- non_zero_coefs[metab_indices]
  names(final_coef_values) <- non_zero_names[metab_indices]
  
  # Performance evaluation
  predictions <- predict(final_model, x_matrix[, final_predictors, drop = FALSE], 
                         s = "lambda.min", type = "response")[,1]
  roc_result <- roc(y_factor, predictions, quiet = TRUE)
  auc_ci <- ci.auc(roc_result)
  
  cat("LASSO model AUC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", 
                                  auc(roc_result), auc_ci[1], auc_ci[3]))
  
} else {
  stop("No significant metabolites entered LASSO analysis")
}

# 5. Calculate metabolite scores and model evaluation ------------------------------------------------
if(length(final_metabolites) > 0 && length(final_coef_values) > 0) {
  
  # Calculate metabolite scores using LASSO coefficients - fixed variable name error
  metabolite_scores <- as.matrix(data_scaled[, final_metabolites]) %*% final_coef_values[final_metabolites]
  data_final <- data.frame(
    Group = data_scaled$Group,
    Score = as.numeric(metabolite_scores),
    data_scaled[, final_metabolites]
  )
  
  # Set cross-validation parameters
  ctrl <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 5,
    summaryFunction = twoClassSummary,
    classProbs = TRUE,
    savePredictions = "final"
  )
  
  cat("\n=== Starting model cross-validation comparison ===\n")
  
  # Method 1: Multi-metabolite model evaluation
  multi_model <- train(
    Group ~ .,
    data = data_final[, c("Group", final_metabolites)],
    method = "glm",
    family = "binomial",
    trControl = ctrl,
    metric = "ROC"
  )
  
  # Method 2: Score model evaluation
  score_model <- train(
    Group ~ Score,
    data = data_final,
    method = "glm",
    family = "binomial",
    trControl = ctrl,
    metric = "ROC"
  )
  
  # Direct ROC analysis (using entire dataset)
  multi_pred <- predict(multi_model$finalModel, 
                        newdata = data_final[, final_metabolites, drop = FALSE],
                        type = "response")
  
  score_pred <- predict(score_model$finalModel,
                        newdata = data.frame(Score = data_final$Score),
                        type = "response")
  
  # Calculate ROC and 95% CI
  roc_multi <- roc(data_final$Group, as.numeric(multi_pred), quiet = TRUE)
  roc_score <- roc(data_final$Group, as.numeric(score_pred), quiet = TRUE)
  
  multi_auc_ci <- ci.auc(roc_multi)
  score_auc_ci <- ci.auc(roc_score)
  
  # DeLong test comparing two models
  delong_test <- roc.test(roc_multi, roc_score, method = "delong")
  
  # 6. Results output -----------------------------------------------------------
  cat("\n=== Final results ===\n")
  cat("Screening process:\n")
  cat("- PLS-DA initial screening:", length(selected_vip), "metabolites\n")
  cat("- Significant differences:", length(significant_metabolites), "metabolites\n") 
  cat("- LASSO final selection:", length(final_metabolites), "metabolites\n\n")
  
  cat("Model performance (95% CI):\n")
  cat("- Multi-metabolite model AUC:", sprintf("%.3f (%.3f-%.3f)\n", 
                                               auc(roc_multi), multi_auc_ci[1], multi_auc_ci[3]))
  cat("- Score model AUC:", sprintf("%.3f (%.3f-%.3f)\n", 
                                    auc(roc_score), score_auc_ci[1], score_auc_ci[3]))
  cat("- AUC difference p-value (DeLong test):", round(delong_test$p.value, 4), "\n\n")
  
  cat("Final metabolites:\n")
  result_df <- data.frame(
    metabolite = names(final_coef_values),
    coefficient = final_coef_values,
    stringsAsFactors = FALSE
  )
  result_df <- result_df[order(-abs(result_df$coefficient)), ]
  print(result_df)
  
  # 7. Visualization --------------------------------------------------------
  par(mfrow = c(1, 2))
  
  # ROC curve comparison
  plot(roc_multi, col = "blue", main = "ROC curve comparison (95% CI)")
  plot(roc_score, col = "red", add = TRUE)
  legend("bottomright", 
         legend = c(
           paste("Multi-metabolite AUC=", round(auc(roc_multi), 3), 
                 " (", round(multi_auc_ci[1], 3), "-", round(multi_auc_ci[3], 3), ")"),
           paste("Score model AUC=", round(auc(roc_score), 3), 
                 " (", round(score_auc_ci[1], 3), "-", round(score_auc_ci[3], 3), ")")
         ),
         col = c("blue", "red"), lwd = 2, cex = 0.8)
  
  # Metabolite coefficient plot
  if(nrow(result_df) > 0) {
    barplot(result_df$coefficient, names.arg = result_df$metabolite,
            horiz = TRUE, las = 1, cex.names = 0.7,
            main = "Final metabolite LASSO coefficients", xlab = "Coefficient value", col = "steelblue")
    abline(v = 0, lty = 2, col = "gray")
  }
  
  par(mfrow = c(1, 1))
  
} else {
  cat("LASSO did not select any metabolites, cannot perform score calculation\n")
} 



data_final <- data_final[,-2]
# 2. Optimize metabolite boxplot
plot_data <- reshape2::melt(data_final, id.vars = "Group")

#gross6#
metab_names <- c(
  "Com_35341_pos" = "Cytosine",
  "Com_8451_pos" = "DHEA",               
  "Com_550_pos" = "Betaine",
  "Com_22881_pos" = "7-Methylxanthine",
  "Com_27317_pos" = "Cholecalciferol",
  "Com_57_neg" = "LPC 18:2",
  "Com_45_neg" = "Elaidic acid"
  
)



# Key correction: Replace metabolite IDs with actual names
plot_data$variable <- metab_names[as.character(plot_data$variable)]


box_plot <- ggplot(plot_data, aes(x = Group, y = value, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, lwd = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, shape = 21) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  labs(y = "log2(Intensity)", 
       subtitle = paste(length(final_metabolites), "significant metabolites for 12-month gross motor"),
       fill = "Group") +  # Set legend title
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",  # Legend at bottom
    legend.title = element_text(size = 10),  # Legend title size
    legend.text = element_text(size = 9),    # Legend text size
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 15),
    panel.spacing = unit(0.8, "lines")
  ) +
  scale_fill_manual(
    values = c("#4E84C4", "#D16103"),  # Colors remain unchanged
    labels = c("Group1", "Group2")     # Clearly label group names
  )

print(box_plot)

