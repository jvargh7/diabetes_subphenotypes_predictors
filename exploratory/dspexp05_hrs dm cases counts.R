rm(list=ls());gc();source(".Rprofile")

#--------------------------------------------------------------------------------------------------------------
### HRS ###
#--------------------------------------------------------------------------------------------------------------
# N = 42,951
hrs_long <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01i_hrs.RDS"))
#-------------------------------------------------------------------------------------------
# New-dm: Self-reported Duration <= 1 year OR HbA1c >= 6.5% OR FPG >=126 mg/dL.
# Existing DM: Self-reported Duration > 1 year
# No DM: Neither of the above
# We need only New-DM and No-DM. We need to exclude Existing DM.

# Identify all DM (either dmagediag is not NA OR (dmagediag is NA, HbA1c >= 6.5))
# N = 10,243
hrs_dm_all <- hrs_long %>%
  group_by(study_id) %>% 
  dplyr::filter((!is.na(dmagediag) | 
                   is.na(dmagediag) & hba1c >= 6.5)) %>%
  ungroup()

# Identify diagnosed DM, N = 7,808
hrs_dm_diag <- hrs_dm_all %>% 
  group_by(study_id) %>%
  dplyr::filter(diagnosed_dm == 1) %>%
  ungroup()

# Among diagnosed DM, duration <= 1 year, N = 2,293
hrs_dm_newdiag <- hrs_dm_diag %>% 
  dplyr::filter(!is.na(dmagediag) & !is.na(age)) %>%
  group_by(study_id) %>% 
  dplyr::filter(between(abs(age - dmagediag), 0, 1)) %>%
  ungroup()

# Identify undiagnosed DM based on A1c. Set dmagediag = curent age, N = 1,175
hrs_dm_undiag <- hrs_long %>%
  group_by(study_id) %>% 
  dplyr::filter((is.na(dmagediag) & hba1c >= 6.5)) %>%
  ungroup() %>% 
  mutate(dmagediag = age)


# Exclude all DM from all HRS to get no DM, N = 32,708
hrs_ndm <- hrs_long %>% 
  dplyr::filter(!study_id %in% hrs_dm_all$study_id) 

# %>% distinct(study_id) %>% nrow()

### HRS Newly diagnosed dm: 3454 ###

#-------------------------------------------------------------------------
# Total sample (no T2D + new T2D), N = 35,987, obs = 289,084
hrs_total <- bind_rows(hrs_dm_newdiag,
                       hrs_dm_undiag,
                       hrs_ndm) %>% 
  mutate(race = haven::as_factor(race),
         ethnicity = haven::as_factor(ethnicity)) %>% 
  mutate(race = case_when(race == "1.white/caucasian" ~ "White",
                          race == "2.black/african american" ~ "Black",
                          race == "3.other" ~ "Other",
                          TRUE ~ "Unknown"),
         ethnicity = case_when(ethnicity == "0.not hispanic" ~ "non-hispanic",
                               ethnicity == "1.hispanic" ~ "hispanic",
                               TRUE ~ "unknown")) %>% 
  mutate(race_clean = race_eth)

saveRDS(hrs_total, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp05_hrs new and no dm.RDS"))

# N = 20263
hrs_female <- hrs_total %>%
  dplyr::filter(gender == "Female")  

attributes(hrs_total$race)$label <- NULL

# 1: White
# 2: Black/African American
# 3: Other
# N = 9835
hrs_racemin <- hrs_total %>%
  dplyr::filter(!is.na(race)) %>% 
  dplyr::filter(race != 1) 
  
hrs_age <- hrs_total %>% 
  dplyr::filter(!is.na(age)) %>% 
  group_by(study_id) %>% 
  slice_max(order_by = wave, with_ties = FALSE) %>%
  ungroup()

hrs_fuptime <- hrs_total %>% 
  dplyr::filter(!is.na(age)) %>% 
  group_by(study_id) %>% 
  mutate(fuptime = max(age, na.rm = TRUE) - min(age, na.rm = TRUE)) %>%
  ungroup()


