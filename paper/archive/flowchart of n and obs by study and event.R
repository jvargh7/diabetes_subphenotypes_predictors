rm(list = ls());gc();source(".Rprofile")

# all participants without diagnosed T2D at baseline
all_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dsppre01_analytic df.RDS")) %>%
  group_by(study, study_id) %>%
  mutate(min_age = min(age)) %>%
  # Restrict to observations within 15 years of earliest age
  dplyr::filter(age <= (min_age + 15)) %>%
  # There are individuals whose dmagediag < minimum age
  # exclude people being diagnosed at baseline (1st visit)
  dplyr::filter(is.na(dmagediag) | ((age <= dmagediag) & (min_age != dmagediag))) %>% 
  
  mutate(event = case_when(# Individuals who are never diagnosed
    is.na(dmagediag) ~ 0,
    
    # Individuals who are diagnosed within 15 years of earliest wave
    dmagediag <= (min_age + 15) ~ 1,
    
    # Individuals who are diagnosed after 15 years of earliest wave
    TRUE ~ 0
    
  )) %>%
  ungroup()

# obs
table(all_df$study)
table(all_df$event, all_df$study)

# N
all_df %>%
  group_by(study) %>%
  summarise(total_unique_ids = n_distinct(study_id))

all_df %>%
  dplyr::filter(event %in% c(0, 1)) %>%
  group_by(study, event) %>%
  summarise(total_unique_ids = n_distinct(study_id), .groups = "drop")
  
#------------------------------------------------------------------------------------
# at least 1 visit time before diagnosis
mt1vst <- all_df %>% 
  dplyr::filter((event == 1 & dmagediag > age) | event == 0)
  
# obs
table(mt1vst$study)
table(mt1vst$event, mt1vst$study)

# N
mt1vst %>%
  group_by(study) %>%
  summarise(total_unique_ids = n_distinct(study_id))

mt1vst %>%
  dplyr::filter(event %in% c(0, 1)) %>%
  group_by(study, event) %>%
  summarise(total_unique_ids = n_distinct(study_id), .groups = "drop")

#------------------------------------------------------------------------------------
# with available cluster
clus_ava <- mt1vst %>% 
  dplyr::filter((event == 1 & !is.na(bmi) & !is.na(hba1c)) | event == 0)
  
# obs
table(clus_ava$study)
table(clus_ava$event, clus_ava$study)

# N
clus_ava %>%
  group_by(study) %>%
  summarise(total_unique_ids = n_distinct(study_id))

clus_ava %>%
  dplyr::filter(event %in% c(0, 1)) %>%
  group_by(study, event) %>%
  summarise(total_unique_ids = n_distinct(study_id), .groups = "drop") 
  
  
  
  