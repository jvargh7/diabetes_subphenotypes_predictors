rm(list = ls());gc();source(".Rprofile")


analytic_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>% 
  mutate(joint_id = as.integer(as.factor(joint_id)))
  # mutate(race = as.factor(race),
  #        # To avoid the warning that 'Imputation method logreg is for categorical data' -- we can convert it back later
  #        female = factor(female,levels=c(0,1))) 

colnames(analytic_df)


# detect outliers
library(purrr)

# detect variables all NA for some people
vars_to_check  <- c("age", "height","weight","bmi","wc","sbp", "dbp","hba1c", 
                     "ldlc","hdlc","glucosef","insulinf","glucose2h",
                     "tgl", "serumcreatinine","homa2b", "homa2ir")


problem_vars <- c()
for (var in vars_to_check) {
  n_all_na <- analytic_df %>%
    group_by(joint_id) %>%
    summarise(all_na = all(is.na(.data[[var]]))) %>%
    dplyr::filter(all_na) %>%
    nrow()
  if (n_all_na > 0) problem_vars <- c(problem_vars, var)
}
print(problem_vars)



multilevel_vars <- c("age", "height","weight","bmi","sbp", "dbp","hba1c")

# proportion_vars <- c("female")
# 
# grouped_vars <- c("race")

# Moved dmagediag to an ID variable
id_vars <- c("study_id", "study", "joint_id","cluster_study_id", "cluster","newdm_event",
             "dmagediag", "t", "earliest_age", "censored_age",
             # no NA
             "female", "race")


library(survey)
library(mice)

before_imputation <- analytic_df  %>% 
  dplyr::select(
    any_of(id_vars),
    any_of(problem_vars),
    any_of(multilevel_vars)
    # any_of(proportion_vars),
    # any_of(grouped_vars)
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



pred[,] <- 1         # Allow all variables as predictors by default

# Corrected below --------
pred[id_vars,] <- 0
pred[,id_vars] <- 0

# For multilevel, set cluster variable (joint_id) to -2 for the multilevel vars
for (v in vars_to_check) {
  if (v %in% problem_vars) {
    pred[v, "joint_id"] <- 0
  } else {
    pred[v, "joint_id"] <- -2
  }
}



mi_dfs <- mice(before_imputation,
               method = method,
               pred = pred,
               m=10,maxit=50,seed=500)

#df <- complete(mi_dfs, action = 1)

saveRDS(mi_dfs, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/mi_dfs.RDS"))


