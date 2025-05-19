rm(list=ls());gc();source(".Rprofile")

#--------------------------------------------------------------------------------------------------------------
### DPP/DPPOS ###
#--------------------------------------------------------------------------------------------------------------
# includes new diagnosed dm + undiagnosed
# N = 3,589
dppos_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01c_dppos.RDS"))

# Identify all DM (either dmagediag is not NA OR (dmagediag is NA, HbA1c >= 6.5))
# N = 2,206
dppos_dm_all <- dppos_longitudinal %>%
  group_by(study_id) %>% 
  dplyr::filter((!is.na(dmagediag) | 
                   is.na(dmagediag) & (hba1c >= 6.5 | glucosef >= 126))) %>%
  ungroup()

# Among diagnosed DM, duration <= 1 year, N = 1,766
dpp_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/dpp_newdm.RDS")) # N = 800
dos_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/dos_newdm.RDS")) %>% # N = 973
  anti_join(dpp_newdm,
            by=c("study_id")) 

dppos_newdm <- dppos_longitudinal %>% 
  dplyr::filter(study_id %in% c(dpp_newdm$study_id, dos_newdm$study_id))

# Identify undiagnosed DM based on A1c. Set agediagnosed_dm = current age, N = 440
dppos_dm_undiag <- dppos_longitudinal %>%
  group_by(study_id) %>% 
  dplyr::filter((is.na(dmagediag) & (hba1c >= 6.5 | glucosef >= 126))) %>%
  ungroup() %>% 
  mutate(dmagediag = age)


# Exclude all DM to get no DM, N = 1,383
dppos_ndm <- dppos_longitudinal %>% 
  dplyr::filter(!study_id %in% dppos_dm_all$study_id) 

# %>% distinct(study_id) %>% nrow()

### ARRIC Newly diagnosed dm: 2206 ###

#-------------------------------------------------------------------------
# Total sample (no T2D + new T2D), N = 3,589, obs = 67,950
dppos_total <- bind_rows(dppos_newdm,
                        dppos_dm_undiag,
                        dppos_ndm)

saveRDS(dppos_total, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp04_dpp new and no dm.RDS"))

# N = 2389
dppos_female <- dppos_total %>%
  dplyr::filter(sex == 2)  

# N = 1509
dppos_racemin <- dppos_total %>%
  dplyr::filter(!is.na(race_eth)) %>% 
  dplyr::filter(race_eth != "NH White") 

dppos_age <- dppos_total %>% 
  dplyr::filter(!is.na(age)) %>% 
  group_by(study_id) %>% 
  slice_max(order_by = age, with_ties = FALSE) %>%
  ungroup()

dppos_fuptime <- dppos_total %>% 
  dplyr::filter(!is.na(age)) %>% 
  group_by(study_id) %>% 
  mutate(fuptime = max(age, na.rm = TRUE) - min(age, na.rm = TRUE)) %>%
  ungroup()

