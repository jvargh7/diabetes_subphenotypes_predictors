rm(list=ls());gc();source("dsphyc3_config.R")


analytic_df <- readRDS(paste0(path_dsp_hyperc3_folder,"/data/dspan01_analytic sample.RDS")) %>% 
  mutate(original_joint_id = paste(study, study_id, sep = "_"),   # for later restoration
         mice_id = as.integer(as.factor(paste(study, study_id, sep = "_"))))  # for multilevel mice

colnames(analytic_df)



problem_vars <- c()

# detect variables all NA for some people
vars_to_check  <- c("age", "height","weight","bmi","wc","sbp", "dbp","hba1c", 
                    "ldlc","hdlc","glucosef","insulinf",
                    "tgl", "serumcreatinine","homa2b", "homa2ir")

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


mi_dfs <- mice(before_imputation,
               method = method,
               predictorMatrix = pred,
               m=10,maxit=50,seed=500)

df <- complete(mi_dfs, action = 1)
cat(paste0("Number of rows: ",nrow(df)))

saveRDS(mi_dfs, paste0(path_dsp_hyperc3_folder,"/working/dsphyc301_mi_dfs.RDS"))




  
 print("Imputation completed")
