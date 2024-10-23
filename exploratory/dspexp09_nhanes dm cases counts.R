rm(list=ls());gc();source(".Rprofile")

#------------------------------------------------------------------------------------------------------------------------------------------------
### Survey - NHANES ###
#--------------------------------------------------------------------------------------------------------------
years <- c("19992000", "20012002", "20032004", "20052006", "20072008", "20092010", "20112012", "20132014", "20152016", "20172018", "2017Mar2020","20212023")

process_nhanes_data <- function(years, path) {
  list_of_data <- list()
  
  for (year in years) {
    file_path <- paste0(path, "/working/cleaned/nhanes_", year, ".rds")
    data <- readRDS(file_path)
    data <- mutate(data, year = paste0(substr(year, 1, 4), "-", substr(year, 5, 8)))
    list_of_data[[year]] <- data
  }
  
  source_df <- bind_rows(list_of_data)
  
  return(source_df)
}

nhanes_data <- process_nhanes_data(years, path_nhanes_ckm_folder) %>% 
  mutate(dm = case_when(dm_doc_told == 1 | glycohemoglobin >= 6.5 | fasting_glucose >= 126 ~ 1,
                        TRUE ~ 0)) 

#-------------------------------------------------------------------------------------------
# Identify all DM (either dm_age is not NA OR (dm_age is NA, HbA1c >= 6.5))
# N = 10636
nhanes_dm_all <- nhanes_data %>% 
  group_by(respondentid) %>% 
  dplyr::filter((!is.na(dm_age) | 
                   is.na(dm_age) & glycohemoglobin >= 6.5)) %>%
  ungroup()

# Identify diagnosed DM, N = 8697
nhanes_dm_diag <- nhanes_dm_all %>% 
  group_by(respondentid) %>%
  dplyr::filter(dm_doc_told == 1) %>%
  ungroup()

# Among diagnosed DM, duration <= 1 year, N = 985
nhanes_dm_newdiag <- nhanes_dm_diag %>% 
  group_by(respondentid) %>% 
  dplyr::filter(!is.na(dm_age) & !is.na(age) & 
                  (age - dm_age) >= 0 & 
                  (age - dm_age) <= 1) %>%
  ungroup()

# Identify undiagnosed DM based on A1c. Set dm_age = current age, N = 1960
nhanes_dm_undiag <- nhanes_data %>%
  group_by(respondentid) %>% 
  dplyr::filter((is.na(dm_age) & glycohemoglobin >= 6.5)) %>%
  ungroup() %>% 
  mutate(dm_age = age)


# Exclude all DM from all HRS to get no DM, N = 93164
nhanes_ndm <- nhanes_data %>% 
  dplyr::filter(!respondentid %in% nhanes_dm_all$respondentid) 


# distinct(respondentid) %>%
# nrow()

### NHANES Newly diagnosed (duration <= 1y + undiagnosed) dm: 2717 ###

#-------------------------------------------------------------------------
# Total sample (no T2D + new T2D), N = 96109, obs = 96109
nhanes_total <- bind_rows(nhanes_dm_newdiag,
                          nhanes_dm_undiag,
                          nhanes_ndm)

# N = 49390
nhanes_female <- nhanes_total %>%
  dplyr::filter(gender == 2)  

# N = 58899
nhanes_racemin <- nhanes_total %>%
  dplyr::filter(race != 3 & !is.na(race)) 

nhanes_age <- nhanes_total %>% 
  dplyr::filter(!is.na(age)) %>% 
  group_by(respondentid) %>% 
  slice_max(order_by = age, with_ties = FALSE) %>%
  ungroup()

nhanes_fuptime <- nhanes_total %>% 
  dplyr::filter(!is.na(age)) %>% 
  group_by(respondentid) %>% 
  mutate(fuptime = max(age, na.rm = TRUE) - min(age, na.rm = TRUE)) %>%
  ungroup()
