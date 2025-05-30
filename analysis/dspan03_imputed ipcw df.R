rm(list = ls());gc();source(".Rprofile")

library(mice)
library(purrr)
library(emmeans)
library(contrast)

source("functions/egfr_ckdepi_2021.R")
mi_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/mi_dfs.RDS"))


wave_dfs <- list()

for(i in 1:mi_dfs$m) {
  df <- complete(mi_dfs, action = i) %>% 
    # no available clusters in these studies
    dplyr::filter(!study %in% c("aric", "cardia")) %>%
    mutate(egfr_ckdepi_2021 = egfr_ckdepi_2021(scr = serumcreatinine,female = female,age = age)) %>% 
    dplyr::filter(newdm_event == 1) # no missing values in subtype
    

  # For OBS when age >= diagnosis age, pick the 1st obs
  exceed_df <- df %>% 
    dplyr::filter(age >= dmagediag) %>% 
    arrange(study_id, study, joint_id, age) %>%
    group_by(study, study_id, joint_id) %>% 
    slice(1) %>% 
    mutate(age = dmagediag)  # N = 1,743
  
  
  endpoint_df <- bind_rows(df %>% 
                             group_by(study, study_id, joint_id) %>% 
                             dplyr::filter(age < dmagediag),
                           exceed_df) %>% 
    arrange(study,study_id,joint_id,age) %>% 
    group_by(study,study_id,joint_id) %>%
    mutate(min_age = min(age)) %>% 
    mutate(event = case_when(
      row_number() != n() ~ 0,  # event is 0 unless it's the last wave
      row_number() == n() ~ 1,  # event is 1 for the last wave
      TRUE ~ NA_real_  # Handle unexpected cases
    )
  ) %>%
  ungroup() # 1,743 events

  
  wave_df <- endpoint_df %>% 
    mutate(censored_age = case_when(is.na(dmagediag) ~ max(age),
                                    
                                    !is.na(dmagediag) & dmagediag <= (min_age + 15) ~ dmagediag,
                                    
                                    TRUE ~ max(age)),
           time_to_event = censored_age - age) %>% 
    ungroup()
 
  
  wave_dfs[[i]] <- wave_df
}

saveRDS(wave_dfs, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan02_predictors analytic dfs.RDS"))



