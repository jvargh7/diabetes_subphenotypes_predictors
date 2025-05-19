rm(list=ls());gc();source(".Rprofile")

#--------------------------------------------------------------------------------------------------------------
### MESA ###
#--------------------------------------------------------------------------------------------------------------
# includes new diagnosed dm + undiagnosed
# N = 6,094
mesa_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01f_mesa.RDS"))

# Identify all DM (either dmagediag is not NA OR (dmagediag is NA, HbA1c >= 6.5 | glucosf2 >= 6.993))
# N = 989
mesa_dm_all <- mesa_longitudinal %>%
  group_by(study_id) %>% 
  dplyr::filter((!is.na(dmagediag) | 
                   is.na(dmagediag) & (hba1c >= 6.5 | glucosef >= 126))) %>%
  ungroup()

# Among diagnosed DM, duration <= 1 year, N = 989
mesa_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/mesa_newdm.RDS"))

mesa_newdm_long <- mesa_longitudinal %>% 
  dplyr::filter(study_id %in% mesa_newdm$study_id)

# Identify undiagnosed DM based on A1c. Set agediagnosed_dm = current age, N = 0
mesa_dm_undiag <- mesa_longitudinal %>%
  group_by(study_id) %>% 
  dplyr::filter((is.na(dmagediag) & (hba1c >= 6.5 | glucosef >= 126))) %>%
  ungroup() %>% 
  mutate(dmagediag = age)


# Exclude all DM to get no DM, N = 5,105
mesa_ndm <- mesa_longitudinal %>% 
  dplyr::filter(!study_id %in% mesa_dm_all$study_id) 

# %>% distinct(study_id) %>% nrow()

### ARRIC Newly diagnosed dm: 989 ###

#-------------------------------------------------------------------------
# Total sample (no T2D + new T2D), N = 6,094, obs = 26,618
mesa_total <- bind_rows(mesa_newdm_long,
                        mesa_dm_undiag,
                        mesa_ndm)

saveRDS(mesa_total, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp08_mesa new and no dm.RDS"))

# N = 3253
mesa_female <- mesa_total %>%
  dplyr::filter(female == 1)  

# N = 3612
mesa_racemin <- mesa_total %>%
  dplyr::filter(!is.na(race)) %>% 
  dplyr::filter(race != "NH White") 

mesa_age <- mesa_total %>% 
  dplyr::filter(!is.na(age)) %>% 
  group_by(study_id) %>% 
  slice_max(order_by = exam, with_ties = FALSE) %>%
  ungroup()

mesa_fuptime <- mesa_total %>% 
  dplyr::filter(!is.na(age)) %>% 
  group_by(study_id) %>% 
  mutate(fuptime = max(age, na.rm = TRUE) - min(age, na.rm = TRUE)) %>%
  ungroup()

