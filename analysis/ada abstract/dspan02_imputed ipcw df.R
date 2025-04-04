rm(list = ls());gc();source(".Rprofile")

library(mice)
library(purrr)
library(emmeans)
library(contrast)

source("functions/egfr_ckdepi_2021.R")
mi_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/mi_dfs.RDS"))


ipcw_dfs <- list()

for(i in 1:mi_dfs$m) {
  df <- complete(mi_dfs, action = i) %>% 
    dplyr::filter(!study %in% c("aric", "cardia"))
  
  fup15y <- df %>%
    mutate(egfr_ckdepi_2021 = egfr_ckdepi_2021(scr = serumcreatinine,female = female,age = age)) %>% 
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
    ungroup() %>% 
    mutate(joint_id = paste(study, study_id, sep = "_"))
  
  # T2D: >= 1 wave before diagnosis; no T2D: >= 2 wave from baseline
  wave_df <- fup15y %>% 
    group_by(study, study_id) %>%
    mutate(has_age_after_min = any(age > min_age),
           has_age_before_dmagediag = any(age < dmagediag)) %>% 
    dplyr::filter(
      (event == 0 & has_age_after_min) |
        (event == 1 & has_age_before_dmagediag)
    ) %>% 
    mutate(censored_age = case_when(is.na(dmagediag) ~ max(age),
                                    
                                    !is.na(dmagediag) & dmagediag <= (min_age + 15) ~ dmagediag,
                                    
                                    TRUE ~ max(age)),
           time_to_event = censored_age - age) %>% 
    select(-has_age_after_min, -has_age_before_dmagediag) %>% 
    ungroup()
  
  # cluster available id at the age of T2D detection
  cluster_avaid <- wave_df %>% 
    mutate(age_int = case_when(study == "dppos" ~ as.integer(age),
                               TRUE ~ age),
           dmagediag_int = case_when(study == "dppos" ~ as.integer(dmagediag),
                               TRUE ~ dmagediag)) %>% 
    group_by(study, study_id) %>%
    dplyr::filter(event == 0 | ((age_int == dmagediag_int) & !is.na(cluster))) %>% 
    ungroup() %>% 
    select(-age_int,-dmagediag_int)
  
  cluster_ava <- wave_df %>% 
    dplyr::filter((joint_id %in% cluster_avaid$joint_id))
  
  
  
  # apply IPCW for people without cluster
  dm_df <- wave_df %>% 
    dplyr::filter(event == 1) %>% 
    mutate(clu_available = case_when(joint_id %in% cluster_avaid$joint_id ~ 1,
                                     TRUE ~ 0))
  
  ltfu_equation <- clu_available ~ study + female + race + min_age
  ltfu_cluster_model <- glm(ltfu_equation, data = dm_df, family = "binomial")
  dm_df$prob_cluster_fup = predict(ltfu_cluster_model,newdata=dm_df,type="response")
  
  
  dm_weight <- dm_df %>%
    mutate(ipcw_cluster = case_when(clu_available == 1 ~ 1/prob_cluster_fup,
                                    TRUE ~ 1/(1-prob_cluster_fup))) %>% 
    distinct(joint_id, .keep_all = TRUE)
  
  ipcw_df <- cluster_ava %>%
    left_join(dm_weight %>% select(study,study_id,joint_id,ipcw_cluster), by = c("joint_id","study","study_id")) %>% 
    mutate(ipcw_cluster = case_when(is.na(ipcw_cluster) ~ 1,
                                    TRUE ~ ipcw_cluster))
  
  
  ipcw_dfs[[i]] <- ipcw_df
}

saveRDS(ipcw_dfs, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan02_ipcw dfs.RDS"))



