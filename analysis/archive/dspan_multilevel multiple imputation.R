rm(list = ls());gc();source(".Rprofile")


analytic_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>% 
  mutate(original_joint_id = paste(study, study_id, sep = "_"),   # for later restoration
         mice_id = as.integer(as.factor(paste(study, study_id, sep = "_"))))  # for multilevel mice

colnames(analytic_df)


# detect outliers
library(purrr)

# detect variables all NA for some people
vars_to_check  <- c("age", "height","weight","bmi","wc","sbp", "dbp","hba1c", 
                     "ldlc","hdlc","glucosef","insulinf",
                     "tgl", "serumcreatinine","homa2b", "homa2ir")


problem_vars <- c()
for (var in vars_to_check) {
  n_all_na <- analytic_df %>%
    group_by(mice_id) %>%
    summarise(all_na = all(is.na(.data[[var]]))) %>%
    # summarise(all_na = all(is.na(.[,var])) %>%
    dplyr::filter(all_na) %>%
    nrow()
  
  if(n_all_na > 0){
    problem_vars = c(problem_vars,var)
  }
  
}

print(problem_vars)



multilevel_vars <- c("age", "height","bmi","sbp", "dbp","hba1c")

# proportion_vars <- c("female")
# 
# grouped_vars <- c("race")

# Moved dmagediag to an ID variable
id_vars <- c("study_id", "study", "mice_id","original_joint_id","cluster_study_id", 
             "cluster","newdm_event","dmagediag", "t", "earliest_age", "censored_age",
             # no NA
             "female", "race")


library(survey)
library(mice)

before_imputation <- analytic_df  %>% 
  dplyr::select(
    any_of(id_vars),
    any_of(problem_vars),
    any_of(multilevel_vars)
  ) 

# Get initial method and pred
mi_null <- mice(before_imputation, maxit = 0)
method = mi_null$method
pred = mi_null$predictorMatrix

# assign methods
method[problem_vars] <- "pmm" 
method[multilevel_vars ] <- "2l.norm"
# method[proportion_vars] <- "2l.bin"
# method[grouped_vars] <- "2l.pmm"
method[id_vars] <- ""

method["weight"] <- "~I(bmi*(height/100)^2)"

pred[c("homa2b","homa2ir"),] <- 0
pred[c("homa2b","homa2ir"),c("insulinf","glucosef")] <- 1



# Initialize predictor matrix
# pred[,] <- 1  # Start by allowing all variables to predict each other

# Clear ID variables from being predictors or outcomes
pred[id_vars,] <- 0
pred[,id_vars] <- 0

# Set cluster indicators (-2) for multilevel variables
for (v in multilevel_vars) {
  if(v %in% colnames(before_imputation)) {
    pred[v, "mice_id"] <- -2
  } 
}

# Handle special cases
if(all(c("homa2b", "homa2ir", "insulinf", "glucosef") %in% colnames(before_imputation))) {
  # HOMA variables should only use insulin and glucose as predictors
  pred[c("homa2b","homa2ir"),] <- 0
  pred[c("homa2b","homa2ir"),c("insulinf","glucosef")] <- 1
} 


# Print predictor matrix for debugging
print("Predictor Matrix:")
print(pred)
print("\nMethods:")
print(method)

mi_dfs <- mice(before_imputation,
               method = method,
               predictorMatrix = pred,
               m=10,maxit=50,seed=500)

df <- complete(mi_dfs, action = 1)

saveRDS(mi_dfs, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/mi_dfs_new.RDS"))





# Correlation heatmap for specified variables
library(corrplot)
library(ggplot2)
library(reshape2)

# Variables for correlation analysis
corr_vars <- c("bmi", "hba1c", "homa2b", "homa2ir", "ldlc", "hdlc", "sbp", "dbp", "egfr_ckdepi_2021")

# Calculate correlation matrix using complete data (first imputed dataset)
df_complete <- complete(mi_dfs, action = 1)

# Select only the variables we want and remove NA values
corr_data <- df_complete %>%
  dplyr::select(all_of(corr_vars)) %>%
  na.omit()

# Calculate correlation matrix
corr_matrix <- cor(corr_data, use = "complete.obs")

# Create correlation heatmap using corrplot
png(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/figures/correlation_heatmap_corrplot.png"), 
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

# Alternative: Create correlation heatmap using ggplot2
corr_melted <- melt(corr_matrix)

p_corr <- ggplot(corr_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(title = "Correlation Heatmap of Clinical Variables")

ggsave(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/figures/correlation_heatmap_ggplot.png"), 
       plot = p_corr, width = 10, height = 8, dpi = 300)

# Print correlation matrix
print("Correlation Matrix:")
print(round(corr_matrix, 3))

# Print summary statistics
cat("\nSample size for correlation analysis:", nrow(corr_data), "\n")
cat("Variables included:", paste(corr_vars, collapse = ", "), "\n")


