rm(list=ls());gc();source(".Rprofile")

#--------------------------------------------------------------------------------------------------------------
### JHS ###
#--------------------------------------------------------------------------------------------------------------
# includes new diagnosed dm + undiagnosed
# N = 1,817
jhs_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01e_jhs.RDS"))

# Identify all DM (either dmagediag is not NA OR (dmagediag is NA, HbA1c >= 6.5))
# N = 282
jhs_dm_all <- jhs_longitudinal %>%
  group_by(study_id) %>% 
  dplyr::filter((!is.na(dmagediag) | 
                   is.na(dmagediag) & (hba1c >= 6.5 | glucosef >= 126))) %>%
  ungroup()

# Among diagnosed DM, duration <= 1 year, N = 268
jhs_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/jhs_newdm.RDS")) 

jhs_newdm_long <- jhs_longitudinal %>% 
  dplyr::filter(study_id %in% jhs_newdm$study_id) 

# Identify undiagnosed DM based on A1c. Set agediagnosed_dm = current age, N = 14
jhs_dm_undiag <- jhs_longitudinal %>%
  group_by(study_id) %>% 
  dplyr::filter((is.na(dmagediag) & (hba1c >= 6.5 | glucosef >= 126))) %>%
  ungroup() %>% 
  mutate(dmagediag = age)


# Exclude all DM to get no DM, N = 1,535
jhs_ndm <- jhs_longitudinal %>% 
  dplyr::filter(!study_id %in% jhs_dm_all$study_id) 

# %>% distinct(study_id) %>% nrow()

### JHS Newly diagnosed dm: 282 ###

#-------------------------------------------------------------------------
# Total sample (no T2D + new T2D), N = 1,817, obs = 4,621
jhs_total <- bind_rows(jhs_newdm_long,
                        jhs_dm_undiag,
                        jhs_ndm)

saveRDS(jhs_total, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp06_jhs new and no dm.RDS"))

# N = 1111
jhs_female <- jhs_total %>%
  dplyr::filter(female == 1)  

# N = 1817
jhs_racemin <- jhs_total %>%
  dplyr::filter(!is.na(race_eth)) %>% 
  dplyr::filter(race_eth != "NH White") 

jhs_age <- jhs_total %>% 
  dplyr::filter(!is.na(age)) %>% 
  group_by(study_id) %>% 
  slice_max(order_by = visit, with_ties = FALSE) %>%
  ungroup()

jhs_fuptime <- jhs_total %>% 
  dplyr::filter(!is.na(age)) %>% 
  group_by(study_id) %>% 
  mutate(fuptime = max(age, na.rm = TRUE) - min(age, na.rm = TRUE)) %>%
  ungroup()

