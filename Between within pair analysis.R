# ==============================================
# Twin Study Analysis: Metabolism and development
# Author: [liqin hu]
# Date: [2025/12/11]
# Description: Analysis of twin data examining relationships between
#              metabolites, birth parameters, and neurodevelopmental outcomes
# ==============================================

# ------------------------------
# 1. LOAD REQUIRED LIBRARIES
# ------------------------------
library(lme4)       # Linear Mixed Effects Models
library(tidyverse)  # Data manipulation and visualization
library(afex)       # Analysis of Factorial Experiments
library(dplyr)      # Data manipulation
library(readr)      # Fast data reading

# ------------------------------
# 2. DATA LOADING AND PREPROCESSING
# ------------------------------

# Load dataset 1
df4 <- read_csv("path.csv")

# Log-transform metabolite columns (columns 33-479)
df4[, 33:479] <- log2(df4[, 33:479])


# Convert to data frame
df4 <- data.frame(df4)

# ------------------------------------------------------------
# 3. Population-level analysis using mixed linear models
# ------------------------------------------------------------

# Initialize results container
result1 <- c()

# Loop through each metabolite column
for (i in 33:479) {
  # Fit mixed linear model: birth weight ~ metabolite + covariates + random effects
  lmer_model <- lmer(birth_weight ~ df4[, i] + gest_age1 + sex + zygosity + (1 | pair), 
                     data = df4)
  
  # Extract coefficients and p-values
  result1 <- rbind(result1, c(colnames(df4)[i], 
                              coef(summary(lmer_model))[2, c(1, 2, 4, 5)]))
}

# Convert results to data frame
result1 <- data.frame(result1)
colnames(result1) <- c("Metabolite", "Estimate", "Std_Error", "t_value", "p_value")

# Apply multiple testing correction (Benjamini-Hochberg)
result1$corrected_p_value <- p.adjust(as.numeric(result1$p_value), method = "BH")

# Save results
write.csv(result1, 
          "path",
          row.names = FALSE)



# ------------------------------
# 4. CALCULATE WITHIN-PAIR DIFFERENCES
# ------------------------------

# 4.1 Calculate differences for metabolite columns
# Get metabolite column names
probe_columns <- colnames(df4)[33:479]

# Convert to numeric
df4[probe_columns] <- lapply(df4[probe_columns], as.numeric)

# Calculate within-pair differences for each metabolite
for (col in probe_columns) {
  new_col_name <- paste0(col, "_dif")  # New column name for differences
  
  df4 <- df4 %>%
    arrange(cohort, pair) %>%          # Sort by cohort and pair
    group_by(cohort) %>%               # Group by cohort
    mutate(!!new_col_name := .data[[col]] - lag(.data[[col]]))  # Calculate difference
}

# Remove grouping
df4 <- df4 %>% ungroup()


# 4.2 Calculate differences for birth weight
birth_weight_col <- colnames(df4)[13]  # Assuming birth weight is in column 13

# Convert to numeric
df4[birth_weight_col] <- lapply(df4[birth_weight_col], as.numeric)

# Calculate within-pair difference for birth weight
new_col_name <- paste0(birth_weight_col, "_dif")
df4 <- df4 %>%
  arrange(cohort, pair) %>%
  group_by(cohort) %>%
  mutate(!!new_col_name := .data[[birth_weight_col]] - lag(.data[[birth_weight_col]]))

df4 <- df4 %>% ungroup()


# ------------------------------
# 5. ANALYZE MONOZYGOTIC (MZ) TWINS
# ------------------------------

# Create MZ twins subset (zygosity = 1 for MZ)
df4_MZ <- subset(df4, zygosity == 1)

# Initialize results data frame
result2 <- data.frame(Column = character(), 
                      Estimate = numeric(), 
                      StdError = numeric(), 
                      PValue = numeric(), 
                      stringsAsFactors = FALSE)

# Loop through difference columns (480:926 are the _dif columns)
for (i in 480:926) {
  # Fit GLM: birth weight difference ~ metabolite difference + covariates
  fit <- glm(birth_weight_dif ~ df4_MZ[[i]] + gest_age1 + sexc, 
             data = df4_MZ)
  
  # Extract results
  result3 <- data.frame(
    Column = colnames(df4_MZ)[i],
    Estimate = coef(summary(fit))[2, 1],
    StdError = coef(summary(fit))[2, 2],
    PValue = coef(summary(fit))[2, 4],
    stringsAsFactors = FALSE
  )
  
  # Append to results
  result2 <- rbind(result2, result3)
}

# Multiple testing correction
result2$corrected_p_value <- p.adjust(result2$PValue, method = "BH")

# Save MZ analysis results
write.csv(result2,
          "path",
          row.names = FALSE)

# ------------------------------
# 6. ANALYZE DIZYGOTIC (DZ) TWINS
# ------------------------------


# Create DZ twins subset (zygosity = 2 for DZ)
df4_DZ <- subset(df4, zygosity == 2)  # Note: variable name changed to df4_DZ for clarity

# Initialize results data frame
result4 <- data.frame(Column = character(), 
                      Estimate = numeric(), 
                      StdError = numeric(), 
                      PValue = numeric(), 
                      stringsAsFactors = FALSE)

# Loop through difference columns
for (i in 480:926) {
  # Fit GLM: birth weight difference ~ metabolite difference + covariates
  fit <- glm(birth_weight_dif ~ df4_DZ[[i]] + gest_age1 + sexc, 
             data = df4_DZ)
  
  # Extract results
  result5 <- data.frame(
    Column = colnames(df4_DZ)[i],
    Estimate = coef(summary(fit))[2, 1],
    StdError = coef(summary(fit))[2, 2],
    PValue = coef(summary(fit))[2, 4],
    stringsAsFactors = FALSE
  )
  
  # Append to results
  result4 <- rbind(result4, result5)
}

# Multiple testing correction
result4$corrected_p_value <- p.adjust(result4$PValue, method = "BH")

# Save DZ analysis results
# Note: Consider renaming output file to indicate DZ analysis
write.csv(result4,
          "path",
          row.names = FALSE)

# ==============================================
# END OF ANALYSIS SCRIPT
# ==============================================


















