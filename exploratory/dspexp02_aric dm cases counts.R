rm(list=ls());gc();source(".Rprofile")

#--------------------------------------------------------------------------------------------------------------
### ARIC ###
#--------------------------------------------------------------------------------------------------------------
# includes new diagnosed dm + undiagnosed
aric_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01a_aric.RDS"))

# Identify all DM (either dmagediag is not NA OR (dmagediag is NA, HbA1c >= 6.5))
# N = 4842
aric_dm_all <- aric_longitudinal %>%
  group_by(study_id) %>% 
  dplyr::filter((!is.na(dmagediag) | 
                   is.na(dmagediag) & hba1c >= 6.5)) %>%
  ungroup()

# Among diagnosed DM, duration <= 1 year, N = 4831
aric_dm_newdiag <- aric_dm_all %>% 
  group_by(study_id) %>% 
  dplyr::filter(!is.na(dmagediag) & !is.na(age) & 
                  (age - dmagediag) >= 0 & 
                  (age - dmagediag) <= 1) %>%
  ungroup()

# Identify undiagnosed DM based on A1c. Set agediagnosed_dm = current age, N = 11
aric_dm_undiag <- aric_longitudinal %>%
  group_by(study_id) %>% 
  dplyr::filter((is.na(dmagediag) & hba1c >= 6.5)) %>%
  ungroup() %>% 
  mutate(dmagediag = age)


# Exclude all DM to get no DM, N = 9454
aric_ndm <- aric_longitudinal %>% 
  dplyr::filter(!study_id %in% aric_dm_all$study_id) 

# distinct(study_id) %>%
# nrow()

### ARRIC Newly diagnosed dm: 4842 ###

#-------------------------------------------------------------------------
# Total sample (no T2D + new T2D), N = 14296, obs = 42158
aric_total <- bind_rows(aric_dm_newdiag,
                       aric_dm_undiag,
                       aric_ndm)

# N = 7500
aric_female <- aric_total %>%
  dplyr::filter(female == 1)  

# N = 3356
aric_racemin <- aric_total %>%
  dplyr::filter(!is.na(race)) %>% 
  dplyr::filter(race != "NH White") 

aric_age <- aric_total %>% 
  dplyr::filter(!is.na(age)) %>% 
  group_by(study_id) %>% 
  slice_max(order_by = visit, with_ties = FALSE) %>%
  ungroup()

aric_fuptime <- aric_total %>% 
  dplyr::filter(!is.na(age)) %>% 
  group_by(study_id) %>% 
  mutate(fuptime = max(age, na.rm = TRUE) - min(age, na.rm = TRUE)) %>%
  ungroup()

