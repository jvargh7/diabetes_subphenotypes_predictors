rm(list=ls());gc();source(".Rprofile")

#--------------------------------------------------------------------------------------------------------------
### HRS ###
#--------------------------------------------------------------------------------------------------------------
biomarkers <- readRDS(paste0(path_g2a_longitudinal_folder,"/working/hrs biomarkers.RDS"))

process_hrs_wave <- function(wave_num, path_folder) {
  male_file <- paste0(path_folder, "/working/hrs/G2A HRS Wave ", wave_num, " male.RDS")
  female_file <- paste0(path_folder, "/working/hrs/G2A HRS Wave ", wave_num, " female.RDS")

  bind_rows(
    readRDS(male_file), 
    readRDS(female_file)
  ) %>%
    dplyr::select(hhid,pn,hhidpn,age,gender,race,diagnosed_dm,agediagnosed_dm,medication_dm) %>%
    mutate(wave = wave_num)
}


hrs_long <- bind_rows(
  lapply(8:14, process_hrs_wave, path_folder = path_g2a_longitudinal_folder)
) %>% 
  left_join(biomarkers %>% 
              dplyr::select(hhid,pn,wave,year,a1c_adj),
            by = c("hhid","pn","wave"))
#-------------------------------------------------------------------------------------------
# New-dm: Self-reported Duration <= 1 year OR HbA1c >= 6.5% OR FPG >=126 mg/dL.
# Existing DM: Self-reported Duration > 1 year
# No DM: Neither of the above
# We need only New-DM and No-DM. We need to exclude Existing DM.

# Identify all DM (either agediagnosed_dm is not NA OR (agediagnosed_dm is NA, HbA1c >= 6.5))
# N = 10243
hrs_dm_all <- hrs_long %>%
  group_by(hhidpn) %>% 
  dplyr::filter((!is.na(agediagnosed_dm) | 
                   is.na(agediagnosed_dm) & a1c_adj >= 6.5)) %>%
  ungroup()

# Identify diagnosed DM, N = 7808
hrs_dm_diag <- hrs_dm_all %>% 
  group_by(hhidpn) %>%
  dplyr::filter(diagnosed_dm == 1) %>%
  ungroup()

# Among diagnosed DM, duration <= 1 year, N = 2279
hrs_dm_newdiag <- hrs_dm_diag %>% 
  group_by(hhidpn) %>% 
  dplyr::filter(!is.na(agediagnosed_dm) & !is.na(age) & 
                  (age - agediagnosed_dm) >= 0 & 
                  (age - agediagnosed_dm) <= 1) %>%
  ungroup()

# Identify undiagnosed DM based on A1c. Set agediagnosed_dm = current age, N = 1175
hrs_dm_undiag <- hrs_long %>%
  group_by(hhidpn) %>% 
  dplyr::filter((is.na(agediagnosed_dm) & a1c_adj >= 6.5)) %>%
  ungroup() %>% 
  mutate(agediagnosed_dm = age)


# Exclude all DM from all HRS to get no DM, N = 32708
hrs_ndm <- hrs_long %>% 
  dplyr::filter(!hhidpn %in% hrs_dm_all$hhidpn) 

# %>%
#   summarise(unique_hhidpn = n_distinct(hhidpn)) %>%
#   pull(unique_hhidpn)

### HRS Newly diagnosed dm: 3454 ###

#-------------------------------------------------------------------------
# Total sample (no T2D + new T2D), N = 35973, obs = 288954
hrs_total <- bind_rows(hrs_dm_newdiag,
                       hrs_dm_undiag,
                       hrs_ndm)

# N = 20263
hrs_female <- hrs_total %>%
  dplyr::filter(gender == "Female")  

attributes(hrs_total$race)$label <- NULL
# N = 9024
hrs_racemin <- hrs_total %>%
  dplyr::filter(!is.na(race)) %>% 
  dplyr::filter(.data$race != 1) 
  
hrs_age <- hrs_total %>% 
  dplyr::filter(!is.na(age)) %>% 
  group_by(hhidpn) %>% 
  slice_max(order_by = wave, with_ties = FALSE) %>%
  ungroup()

hrs_fuptime <- hrs_total %>% 
  dplyr::filter(!is.na(age)) %>% 
  group_by(hhidpn) %>% 
  mutate(fuptime = max(age, na.rm = TRUE) - min(age, na.rm = TRUE)) %>%
  ungroup()
































