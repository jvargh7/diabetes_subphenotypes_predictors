rm(list=ls());gc();source(".Rprofile")

#--------------------------------------------------------------------------------------------------------------
### MESA ###
#--------------------------------------------------------------------------------------------------------------
# includes new diagnosed dm + undiagnosed
mesa_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01f_mesa.RDS"))

# Identify all DM (either dmagediag is not NA OR (dmagediag is NA, HbA1c >= 6.5 | glucosf2 >= 6.993))
# N = 989
mesa_dm_all <- mesa_longitudinal %>%
  group_by(study_id) %>% 
  dplyr::filter((!is.na(dmagediag) | 
                   is.na(dmagediag) & (hba1c >= 6.5 | glucosef2 >= 6.993))) %>%
  ungroup()

# Among diagnosed DM, duration <= 1 year, N = 989
mesa_dm_newdiag <- mesa_dm_all %>% 
  group_by(study_id) %>% 
  dplyr::filter(!is.na(dmagediag) & !is.na(age) & 
                  (age - dmagediag) >= 0 & 
                  (age - dmagediag) <= 1) %>%
  ungroup()

# Identify undiagnosed DM based on A1c. Set agediagnosed_dm = current age, N = 0
mesa_dm_undiag <- mesa_longitudinal %>%
  group_by(study_id) %>% 
  dplyr::filter((is.na(dmagediag) & (hba1c >= 6.5 | glucosef2 >= 6.993))) %>%
  ungroup() %>% 
  mutate(dmagediag = age)


# Exclude all DM to get no DM, N = 5105
mesa_ndm <- mesa_longitudinal %>% 
  dplyr::filter(!study_id %in% mesa_dm_all$study_id) 

# distinct(study_id) %>%
# nrow()

### ARRIC Newly diagnosed dm: 989 ###

#-------------------------------------------------------------------------
# Total sample (no T2D + new T2D), N = 6094, obs = 23334
mesa_total <- bind_rows(mesa_dm_newdiag,
                        mesa_dm_undiag,
                        mesa_ndm)

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

