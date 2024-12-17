rm(list = ls());gc();source(".Rprofile")

# cases with at least one wave before diagnosis
g1wave_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dsppre01_analytic df.RDS")) %>% 
  

# cases with at least one wave before diagnosis + available cluster

# apply IPCW


# multiple imputation