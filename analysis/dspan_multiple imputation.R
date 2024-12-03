rm(list = ls());gc();source(".Rprofile")

analytic_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01_analytic df.RDS")) %>% 
  mutate(race_eth = as.factor(race_eth),
         # To avoid the warning that 'Imputation method logreg is for categorical data' -- we can convert it back later
         female = factor(female,levels=c(0,1)))

colnames(analytic_df)

# Why do you have both insulinf, glucosef and insulinf2, glucosef2?!
continuous_vars <- c("age", "sbp", "dbp", "height", "wc", "bmi", "hba1c", "insulinf",
                     "glucosef", "glucose2h", "tgl", "hdlc", "ldlc", "serumcreatinine", "urinecreatinine",
                     "egfr", "apo_a", "apo_b", "uric_acid", "vldlc", 
                     "hc", "triceps", "iliac", "abdominal", "medial", "ast", "alt",
                     "insulinf2", "glucosef2", "urinealbumin", "uacr", "weight", "homa2b", "homa2ir")

proportion_vars <- c("female")

grouped_vars <- c("study","race_eth")

# Moved dmagediag to an ID variable
id_vars <- c("study_id",  "visit", "year", "exam", "cluster_study_id", "newdm", "dmagediag","diagDays", "anthro_StudyDays", "lab_StudyDays")

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
# method[method == "pmm"] <- "rf" # Takes too long
method[proportion_vars] <- "logreg"
method[grouped_vars] <- "polyreg"
method[id_vars] <- ""

# Made glucosef and insulinf dependent on glucosef2 and insulinf2
method["glucosef"] <- "~I(glucosef2/0.0555)"
method["insulinf"] <- "~I(insulinf2/6)"
method["weight"] ~ "I(bmi*(height/100)^2)"

pred = mi_null$predictorMatrix

# Corrected below --------
pred[id_vars,] <- 0
pred[,id_vars] <- 0


mi_dfs <- mice(before_imputation,
               method = method,
               pred = pred,
               m=1,maxit=50,seed=500)

df <- complete(mi_dfs, action = 1)

saveRDS(df, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/mi_dfs.RDS"))
