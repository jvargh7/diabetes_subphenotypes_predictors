rm(list = ls());gc();source(".Rprofile")

library(mice)
library(purrr)
library(emmeans)
library(contrast)

mi_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/mi_dfs.RDS"))

for(i in 1:mi_dfs$m) {
  df <- complete(mi_dfs, action = 1) 
  
  fup15y <- df %>%
    group_by(study, study_id) %>%
    mutate(min_age = min(age)) %>%
    # Restrict to observations within 15 years of earliest age
    dplyr::filter(age <= (min_age + 15)) %>%
    mutate(event = case_when(# Individuals who are never diagnosed
      is.na(dmagediag) ~ 0,
      
      # Individuals who are diagnosed within 15 years of earliest wave
      dmagediag <= (min_age + 15) ~ 1,
      
      # Individuals who are diagnosed after 15 years of earliest wave
      TRUE ~ 0
      
    )) %>% 
    # without diagnosed T2D at baseline 
    dplyr::filter(event == 0 | ((age <= dmagediag) & (min_age < dmagediag))) %>% 
    ungroup()
  
  # cluster available at the age of T2D detection
  cluster_ava <- fup15y %>% 
    mutate(age_int = case_when(study == "dppos" ~ as.integer(age),
                               TRUE ~ age),
           dmagediag_int = case_when(study == "dppos" ~ as.integer(dmagediag),
                               TRUE ~ dmagediag)) %>% 
    group_by(study, study_id) %>%
    dplyr::filter(event == 0 | ((age_int == dmagediag_int) & !is.na(cluster))) %>% 
    ungroup()
  
  # T2D: >= 1 wave before diagnosis; no T2D: >= 2 wave from baseline
  wave_df <- cluster_ava %>% 
    group_by(study, study_id) %>%
    dplyr::filter(
      (event == 0 & any(age > min_age)) |
        (event == 1 & any(age > min_age & age < dmagediag))
    )
    ungroup()
  
  analytic_dfs[[i]] <- wave_df
}



# cases with at least one wave before diagnosis

  

# cases with at least one wave before diagnosis + available cluster

# apply IPCW


# multiple imputation