rm(list=ls());gc();source(".Rprofile")

#--------------------------------------------------------------------------------------------------------------
### cardia ###
#--------------------------------------------------------------------------------------------------------------
# includes new diagnosed dm + undiagnosed
# N = 5,054
cardia_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01b_cardia.RDS"))

# Identify all DM (either dmagediag is not NA OR (dmagediag is NA, HbA1c >= 6.5 | fasting glucose >= 126))
# N = 769
cardia_dm_all <- cardia_longitudinal %>%
  group_by(study_id) %>% 
  dplyr::filter((!is.na(dmagediag) | 
                   is.na(dmagediag) & (hba1c >= 6.5 | glucosef >= 126))) %>%
  ungroup()

# Among diagnosed DM, duration <= 1 year, N = 623
cardia_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/cardia_newdm.RDS")) 

cardia_newdm_long <- cardia_longitudinal %>% 
  dplyr::filter(study_id %in% cardia_newdm$study_id) 


# Identify undiagnosed DM based on A1c. Set agediagnosed_dm = current age, N = 146
cardia_dm_undiag <- cardia_longitudinal %>%
  group_by(study_id) %>% 
  dplyr::filter((is.na(dmagediag) & (hba1c >= 6.5 | glucosef >= 126))) %>%
  ungroup() %>% 
  mutate(dmagediag = age)


# Exclude all DM to get no DM, N = 4,285
cardia_ndm <- cardia_longitudinal %>% 
  dplyr::filter(!study_id %in% cardia_dm_all$study_id) 

# distinct(study_id) %>% nrow()

### CARDIA Newly diagnosed dm: 769 ###

#-------------------------------------------------------------------------
# Total sample (no T2D + new T2D), N = 5,054, obs = 35,365
cardia_total <- bind_rows(cardia_newdm_long,
                        cardia_dm_undiag,
                        cardia_ndm)

saveRDS(cardia_total, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp03_cardia new and no dm.RDS"))

# N = 2738
cardia_female <- cardia_total %>%
  dplyr::filter(female == 1)  

# N = 2592
cardia_racemin <- cardia_total %>%
  dplyr::filter(!is.na(race)) %>% 
  dplyr::filter(race != "NH White") 

cardia_age <- cardia_total %>% 
  dplyr::filter(!is.na(age)) %>% 
  group_by(study_id) %>% 
  slice_max(order_by = year, with_ties = FALSE) %>%
  ungroup()

cardia_fuptime <- cardia_total %>% 
  dplyr::filter(!is.na(age)) %>% 
  group_by(study_id) %>% 
  mutate(fuptime = max(age, na.rm = TRUE) - min(age, na.rm = TRUE)) %>%
  ungroup()

