rm(list=ls());gc();source(".Rprofile")

#--------------------------------------------------------------------------------------------------------------
### ARIC ###
#--------------------------------------------------------------------------------------------------------------
# includes new diagnosed dm + undiagnosed
# N =  13,817
aric_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01a_aric.RDS"))

# Identify all DM (either dmagediag is not NA OR (dmagediag is NA, HbA1c >= 6.5 | fasting glucose >= 126))
# N = 4,386
aric_dm_all <- aric_longitudinal %>%
  group_by(study_id) %>% 
  dplyr::filter((!is.na(dmagediag) | 
                   is.na(dmagediag) & (hba1c >= 6.5 | glucosef >= 126))) %>%
  ungroup()

# Among diagnosed DM, duration <= 1 year, N = 4,352
aric_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/aric_newdm.RDS")) 

aric_newdm_long <- aric_longitudinal %>% 
  dplyr::filter(study_id %in% aric_newdm$study_id) 

# Identify undiagnosed DM based on A1c. Set agediagnosed_dm = current age, N = 34
aric_dm_undiag <- aric_longitudinal %>%
  group_by(study_id) %>% 
  dplyr::filter((is.na(dmagediag) & (hba1c >= 6.5 | glucosef >= 126))) %>%
  ungroup() %>% 
  mutate(dmagediag = age)


# Exclude all DM to get no DM, N = 9,431
aric_ndm <- aric_longitudinal %>% 
  dplyr::filter(!study_id %in% aric_dm_all$study_id) 

# distinct(study_id) %>%
# nrow()

### ARRIC Newly diagnosed dm: 4386 ###

#-------------------------------------------------------------------------
# Total sample (no T2D + new T2D), N = 13,817, obs = 57,482
aric_total <- bind_rows(aric_newdm_long,
                       aric_dm_undiag,
                       aric_ndm) 

saveRDS(aric_total, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp02_aric new and no dm.RDS"))

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

