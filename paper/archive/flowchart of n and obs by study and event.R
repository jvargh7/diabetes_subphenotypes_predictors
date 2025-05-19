rm(list = ls());gc();source(".Rprofile")

# restrict to 15y follow-up time for all 
fup15y <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dsppre01_analytic df.RDS")) %>%
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
  ungroup() %>% 
  # create unique id, easier to count N; N = 30374
  mutate(joint_id = paste(study, study_id, sep = "_"))

#------------------------------------------------------------
# no T2D within 15y follow-up time - censored, N = 24310
nodm15y <- fup15y %>% 
  group_by(study, study_id) %>%
  dplyr::filter(event == 0) %>% 
  ungroup()

# obs
table(nodm15y$study)
table(nodm15y$event, nodm15y$study)

# N
nodm15y %>%
  group_by(study) %>%
  summarise(total_unique_ids = n_distinct(joint_id))

nodm15y %>%
  dplyr::filter(event %in% c(0, 1)) %>%
  group_by(study, event) %>%
  summarise(total_unique_ids = n_distinct(joint_id), .groups = "drop")


# at least have 2 wave
nodm2wv <- nodm15y %>% 
  group_by(study, study_id) %>%
  mutate(has_age_after_min = any(age > min_age)) %>% 
  dplyr::filter(has_age_after_min) %>% 
  ungroup()

# obs
table(nodm2wv$study)
table(nodm2wv$event, nodm2wv$study)

# N
nodm2wv %>%
  group_by(study) %>%
  summarise(total_unique_ids = n_distinct(joint_id))

nodm2wv %>%
  dplyr::filter(event %in% c(0, 1)) %>%
  group_by(study, event) %>%
  summarise(total_unique_ids = n_distinct(joint_id), .groups = "drop")

#------------------------------------------------------------
# T2D within 15y follow-up time

# not diagnosed at baseline, N = 4933
dm15y <- fup15y %>% 
  group_by(study, study_id) %>%
  dplyr::filter(event == 1) %>% 
  # There are individuals whose dmagediag < minimum age
  # exclude people being diagnosed at baseline (1st visit)
  dplyr::filter((age <= dmagediag) & (min_age < dmagediag)) %>% 
  ungroup()

# obs
table(dm15y$study)
table(dm15y$event, dm15y$study)

# N
dm15y %>%
  group_by(study) %>%
  summarise(total_unique_ids = n_distinct(joint_id))

dm15y %>%
  dplyr::filter(event %in% c(0, 1)) %>%
  group_by(study, event) %>%
  summarise(total_unique_ids = n_distinct(joint_id), .groups = "drop")


# available bmi & a1c at diagnosis age, N = 748
clus_avaid <- dm15y %>% 
  dplyr::filter(study != "dppos") %>% 
  group_by(study, study_id) %>%
  dplyr::filter((age == dmagediag) & !is.na(cluster)) %>% 
  ungroup()

# N = 1387
clus_avaid_dpp <- dm15y %>% 
  dplyr::filter(study == "dppos") %>% 
  mutate(age_int = as.integer(age),
         dmagediag_int = as.integer(dmagediag)) %>% 
  group_by(study, study_id) %>%
  dplyr::filter((age_int == dmagediag_int) & !is.na(cluster)) %>% 
  ungroup()

# N = 2135
clus_ava <- dm15y %>% 
  dplyr::filter((joint_id %in% clus_avaid$joint_id) | (joint_id %in% clus_avaid_dpp$joint_id))

# obs
table(clus_ava$study)
table(clus_ava$event, clus_ava$study)

# N
clus_ava %>%
  group_by(study) %>%
  summarise(total_unique_ids = n_distinct(joint_id))

clus_ava %>%
  dplyr::filter(event %in% c(0, 1)) %>%
  group_by(study, event) %>%
  summarise(total_unique_ids = n_distinct(joint_id), .groups = "drop")

#------------------------------------------------------------------------
# 64 + 190 = 254 ppl don't have visit age that at the same year of their diagnosis
df <- dm15y %>% 
  dplyr::filter(study != "dppos") %>% 
  group_by(study, study_id) %>%
  mutate(no_age_equals_dmagediag = all(age != dmagediag)) %>%
  # Filter to keep only groups where no age equals dmagediag
  dplyr::filter(no_age_equals_dmagediag) %>%
  # Remove the helper column as it's no longer needed
  select(study,study_id,joint_id,age,dmagediag,no_age_equals_dmagediag) %>%
  # Ungroup to prevent accidental grouping effects in later analysis
  ungroup()

df_dpp <- dm15y %>% 
  dplyr::filter(study == "dppos") %>% 
  mutate(age_int = as.integer(age),
         dmagediag_int = as.integer(dmagediag)) %>% 
  group_by(study, study_id) %>%
  mutate(no_age_equals_dmagediag = all(age_int != dmagediag_int)) %>%
  # Filter to keep only groups where no age equals dmagediag
  dplyr::filter(no_age_equals_dmagediag) %>%
  # Remove the helper column as it's no longer needed
  select(study,study_id,joint_id,age,dmagediag,no_age_equals_dmagediag) %>%
  # Ungroup to prevent accidental grouping effects in later analysis
  ungroup()

#------------------------------------------------------------------------
# at least have 1 wave before diagnosis, N = 1731
dm1wv <- clus_ava %>% 
  group_by(study, study_id) %>%
  mutate(has_age_before_dmagediag = any(age < dmagediag)) %>%
  # Filter groups that have at least one valid age record before dmagediag
  dplyr::filter(has_age_before_dmagediag) %>%
  # Remove the helper column if not needed further
  select(-has_age_before_dmagediag) %>%
  ungroup()

# obs
table(dm1wv$study)
table(dm1wv$event, dm1wv$study)

# N
dm1wv %>%
  group_by(study) %>%
  summarise(total_unique_ids = n_distinct(joint_id))

dm1wv %>%
  dplyr::filter(event %in% c(0, 1)) %>%
  group_by(study, event) %>%
  summarise(total_unique_ids = n_distinct(joint_id), .groups = "drop")








