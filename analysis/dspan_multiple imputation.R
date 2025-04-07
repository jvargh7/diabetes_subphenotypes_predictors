rm(list = ls());gc();source(".Rprofile")


analytic_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/8 cohorts/dsppre01_analytic df.RDS")) %>% 
  mutate(race_clean = as.factor(race_clean),
         # To avoid the warning that 'Imputation method logreg is for categorical data' -- we can convert it back later
         female = factor(female,levels=c(0,1)),
         med_dm_use = as.factor(med_dm_use),
         med_bp_use = as.factor(med_bp_use),
         med_chol_use = as.factor(med_chol_use)) %>% 
  select(-c("apo_a","apo_b","alt","sex"))

colnames(analytic_df)


# detect outliers
library(purrr)

continuous_vars <- c("age", "height","weight","bmi","wc","sbp", "dbp","hba1c", 
                     "totalc","ldlc","hdlc","vldlc","glucosef","insulinf","glucose2h",
                     "tgl", "serumcreatinine","urinecreatinine","urinealbumin",
                     "uacr","egfr", "homa2b", "homa2ir")

proportion_vars <- c("female","med_dm_use","med_bp_use","med_chol_use")

grouped_vars <- c("race_clean")

# Moved dmagediag to an ID variable
id_vars <- c("study_id", "study", "wave", "cluster_study_id", "cluster", "new_id", "ratio_th",
             "dmagediag", "dmduration", "dmfamilyhistory","available_labs", "available_anthro")


summary_stats <- analytic_df %>%
  summarise(across(all_of(continuous_vars), list(
    Min = ~min(., na.rm = TRUE),
    `1st Qu.` = ~quantile(., 0.25, na.rm = TRUE),
    Median = ~median(., na.rm = TRUE),
    Mean = ~mean(., na.rm = TRUE),
    `3rd Qu.` = ~quantile(., 0.75, na.rm = TRUE),
    Max = ~max(., na.rm = TRUE),
    NA_s = ~sum(is.na(.))
  ), .names = "{.col}_{.fn}")) %>%
  pivot_longer(cols = everything(), names_to = "name", values_to = "value") %>%
  separate(name, into = c("variable", "stat"), sep = "_") %>%
  pivot_wider(names_from = stat, values_from = value)


write.csv(summary_stats, "analysis/check outliers before mice.csv")


library(survey)
library(mice)

before_imputation <- analytic_df  %>% 
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

method["weight"] <- "~I(bmi*(height/100)^2)"

pred = mi_null$predictorMatrix

# Corrected below --------
pred[id_vars,] <- 0
pred[,id_vars] <- 0


mi_dfs <- mice(before_imputation,
               method = method,
               pred = pred,
               m=10,maxit=50,seed=500)

#df <- complete(mi_dfs, action = 1)

saveRDS(mi_dfs, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/8 cohorts/mi_dfs.RDS"))
