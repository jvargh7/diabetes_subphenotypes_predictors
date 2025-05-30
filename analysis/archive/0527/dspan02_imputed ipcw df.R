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
    # no available clusters in these studies
    dplyr::filter(!study %in% c("aric", "cardia")) 
  
  
  fup15y <- df %>%
    mutate(joint_id = paste(study, study_id, sep = "_")) %>% 
    mutate(
      egfr_ckdepi_2021 = egfr_ckdepi_2021(scr = serumcreatinine,female = female,age = age),
      # ratio_th = case_when(is.na(ratio_th) ~ tgl/hdlc,
      #                           TRUE ~ ratio_th),
      uacr = case_when(is.na(uacr) ~ urinealbumin/urinecreatinine,
                           TRUE ~ uacr)) %>% 
    group_by(study, study_id) %>%
    mutate(min_age = min(age)) %>%
    mutate(dmage_filter = case_when(# Individuals who are never diagnosed
      is.na(dmagediag) ~ 0,
      
      # Individuals who are diagnosed within 15 years of earliest wave
      dmagediag <= (min_age + 15) ~ 1,
      
      # Individuals who are diagnosed after 15 years of earliest wave
      TRUE ~ 0
      
    )) %>% 
    # Restrict to observations within 15 years of earliest age ------------ automatically
    # without diagnosed T2D at baseline 
    dplyr::filter((dmage_filter == 0 & age <= (min_age + 15)) | 
                    (dmage_filter == 1 & min_age < dmagediag)) %>% 
    ungroup() 
  
  
  # For OBS when age >= diagnosis age, pick the 1st obs
  exceed_df <- fup15y %>% 
    dplyr::filter(dmage_filter == 1 & age >= dmagediag) %>% 
    arrange(study_id, study, age) %>%
    group_by(study, study_id) %>% 
    slice(1) %>% 
    mutate(age = dmagediag)  # N = 2,708
  
  
  endpoint_df <- bind_rows(fup15y %>% 
                             group_by(study, study_id) %>% 
                             dplyr::filter(dmage_filter == 0 | age < dmagediag),
                           exceed_df) %>% 
    arrange(study,study_id,age) %>% 
    group_by(study,study_id) %>%
    mutate(event = case_when(
      dmage_filter == 0 ~ 0,
      dmage_filter == 1 & row_number() != n() ~ 0,  # For dmage_filter == 1, event is 0 unless it's the last wave
      dmage_filter == 1 & row_number() == n() ~ 1,  # For dmage_filter == 1, event is 1 for the last wave
      TRUE ~ NA_real_  # Handle unexpected cases
    )
  ) %>%
  ungroup() # 2,827 events

  
  # T2D: >= 1 wave before diagnosis; no T2D: >= 2 wave from baseline
  wave_df <- endpoint_df %>% 
    arrange(study,study_id,age) %>% 
    group_by(study, study_id) %>%
    mutate(has_age_after_min = any(age > min_age),
           has_age_before_dmagediag = any(age < dmagediag)) %>% 
    dplyr::filter(
      (dmage_filter == 0 & has_age_after_min) |
        (dmage_filter == 1 & has_age_before_dmagediag)
    ) %>% # N = 2,644 events
    mutate(censored_age = case_when(is.na(dmagediag) ~ max(age),
                                    
                                    !is.na(dmagediag) & dmagediag <= (min_age + 15) ~ dmagediag,
                                    
                                    TRUE ~ max(age)),
           time_to_event = censored_age - age) %>% 
    select(-has_age_after_min, -has_age_before_dmagediag) %>% 
    ungroup()
  
  # cluster available id at the age of T2D detection, N = 8,173
  cluster_avaid <- wave_df %>% 
    mutate(age_int = case_when(study == "dppos" ~ as.integer(age),
                               TRUE ~ age),
           dmagediag_int = case_when(study == "dppos" ~ as.integer(dmagediag),
                               TRUE ~ dmagediag)) %>% 
    group_by(study, study_id) %>%
    dplyr::filter(dmage_filter == 0 | ((age_int == dmagediag_int) & !is.na(cluster))) %>% 
    ungroup() %>% 
    select(-age_int,-dmagediag_int)
  # dppos: 2,529, jhs: 1,476, mesa: 5,358
  
  cluster_ava <- wave_df %>% 
    dplyr::filter((joint_id %in% cluster_avaid$joint_id))
  

  
  # apply IPCW for people without cluster; 1,910 available vs 717 unavailable
  dm_df <- wave_df %>% 
    dplyr::filter(event == 1) %>% 
    mutate(clu_available = case_when(joint_id %in% cluster_avaid$joint_id ~ 1,
                                     TRUE ~ 0))
  
  # Ensure all categorical variables are factors and have correct levels if new levels might be present
  dm_df$study <- factor(dm_df$study, levels = unique(dm_df$study))
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



