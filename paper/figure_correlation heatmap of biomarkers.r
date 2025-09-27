rm(list = ls());gc();source(".Rprofile")


# Baseline
baseline_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>%
  dplyr::filter(age == earliest_age)



# Correlation heatmap for specified variables
library(corrplot)
library(ggplot2)
library(reshape2)

# Variables for correlation analysis
corr_vars <- c("bmi", "hba1c", "homa2b", "homa2ir", "ldlc", "hdlc", "tgl", "sbp", "dbp", "egfr_ckdepi_2021")


# Select only the variables we want and remove NA values
corr_data <- baseline_df %>%
  dplyr::select(all_of(corr_vars)) %>%
  na.omit()

# Calculate correlation matrix
corr_matrix <- cor(corr_data, use = "complete.obs")

# Create variable labels for better display
var_labels <- c(
  "bmi" = "BMI (kg/m²)",
  "hba1c" = "HbA1c (%)",
  "homa2b" = "HOMA2-B (%)",
  "homa2ir" = "HOMA2-IR",
  "ldlc" = "LDL (mg/dL)",
  "hdlc" = "HDL (mg/dL)",
  "tgl" = "Triglycerides (mg/dL)",
  "sbp" = "SBP (mmHg)",
  "dbp" = "DBP (mmHg)",
  "egfr_ckdepi_2021" = "eGFR (mL/min/1.73m²)"
)

# Apply labels to correlation matrix
rownames(corr_matrix) <- var_labels[rownames(corr_matrix)]
colnames(corr_matrix) <- var_labels[colnames(corr_matrix)]

# Create correlation heatmap using corrplot
png(paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/correlation_heatmap_corrplot.png"), 
    width = 800, height = 800, res = 120)
corrplot(corr_matrix, 
         method = "color", 
         type = "upper", 
         order = "hclust",
         addCoef.col = "black", 
         tl.col = "black", 
         tl.srt = 45,
         number.cex = 0.7,
         col = colorRampPalette(c("blue", "white", "red"))(200))
dev.off()

