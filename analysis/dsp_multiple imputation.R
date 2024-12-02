rm(list = ls());gc();source(".Rprofile")

analytic_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01_analytic df.RDS")) %>% 
  mutate(race_eth = as.factor(race_eth))

colnames(analytic_df)

continuous_vars <- c("age", "dmagediag", "sbp", "dbp", "height", "wc", "bmi", "hba1c", "insulinf",
                     "glucosef", "glucose2h", "tgl", "hdlc", "ldlc", "serumcreatinine", "urinecreatinine",
                     "egfr", "apo_a", "apo_b", "uric_acid", "vldlc", "diagDays", "lab_StudyDays", 
                     "anthro_StudyDays", "hc", "triceps", "iliac", "abdominal", "medial", "ast", "alt",
                     "insulinf2", "glucosef2", "urinealbumin", "uacr", "weight", "homa2b", "homa2ir")

proportion_vars <- c("female")

grouped_vars <- c("race_eth")

id_vars <- c("study_id", "study", "visit", "year", "exam", "cluster_study_id", "newdm")

library(survey)
library(mice)

before_imputation <- analytic_df %>% 
  dplyr::select(
    any_of(id_vars),
    any_of(continuous_vars),
    any_of(proportion_vars),
    any_of(grouped_vars)
  )

mi_null <- mice(before_imputation, maxit = 0)

method = mi_null$method
method[proportion_vars] <- "logreg"

pred = mi_null$predictorMatrix

id_var_indices <- which(colnames(pred) %in% id_vars)
pred[id_var_indices] <- 0
pred[,id_var_indices] <- 0


mi_dfs <- mice(before_imputation,
               method = method,
               pred = pred,
               m=1,maxit=50,seed=500)

df <- complete(mi_dfs, action = 1)

saveRDS(df, "analysis/mi_dfs.RDS")
