rm(list=ls());gc();source(".Rprofile")
library(haven)
accord_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/accord_newdm.RDS")) %>% 
  rename(race = race_eth) %>% 
  dplyr::select(study_id, study, race, female) %>% 
  mutate(study = "accord")

aric_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp02_aric new and no dm.RDS")) %>% 
  mutate(study_id = as.numeric(study_id),
         study = "aric") %>% 
  dplyr::select(study_id, study, race, female)

cardia_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp03_cardia new and no dm.RDS")) %>% 
  mutate(study = "cardia") %>% 
  dplyr::select(study_id, study, race, female)

dppos_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp04_dpp new and no dm.RDS")) %>% 
  rename(race = race_eth) %>% 
  mutate(female = case_when(sex == 1 ~ 0,
                            sex == 2 ~ 1,
                            TRUE ~ NA_real_),
         study = "dpp/dppos") %>% 
  dplyr::select(study_id, study, race, female)

attributes(hrs_total$race)$label <- NULL
hrs_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp05_hrs new and no dm.RDS")) %>% 
  mutate(study_id = as.numeric(hhidpn),
         study = "hrs",
         race = as.numeric(race),
         race = case_when(race == 1 ~ "NH White",
                          race == 2 ~ "NH Black",
                          race == 3 ~ "NH Other",
                          TRUE ~ NA_character_))

jhs_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp06_jhs new and no dm.RDS")) %>% 
  rename(race = race_eth) %>%
  mutate(study = "jhs") %>% 
  dplyr::select(study_id, study, race, female)

la_newdm <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/la_newdm.RDS")) %>% 
  rename(race = race_eth) %>%
  mutate(study = "look ahead") %>% 
  dplyr::select(study_id, study, race, female)

mesa_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp08_mesa new and no dm.RDS")) %>% 
  mutate(race = case_when(race == "Other" ~ "NH Other",
                          TRUE ~ race),
         study = "mesa") %>% 
  dplyr::select(study_id, study, race, female)

nhanes_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp09_nhanes new and no dm.RDS"))  %>% 
  mutate(study_id = respondentid,
         study = "nhanes",
         race = case_when(race == 3 ~ "NH White",
                          race == 4 ~ "NH Black",
                          race == 6 ~ "Asian",
                          race == 1 | race == 2 ~ "Hispanic",
                          TRUE ~ "NH Other"),
         female = case_when(gender == 1 ~ 0,
                            gender == 2 ~ 1,
                            TRUE ~ NA_real_)) %>% 
  dplyr::select(study_id, study, race, female)

oneflorida <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp10_oneflorida new and no dm.RDS"))  %>% 
  mutate(study_id = as.numeric(ID),
         study = "oneflorida",
         race = case_when(nhwhite == 1 & nhblack == 0 & hispanic == 0 & nhother == 0 ~ "NH White",
                          nhwhite == 0 & nhblack == 1 & hispanic == 0 & nhother == 0 ~ "NH Black",
                          nhwhite == 0 & nhblack == 0 & hispanic == 1 & nhother == 0 ~ "Hispanic",
                          nhwhite == 0 & nhblack == 0 & hispanic == 0 & nhother == 1 ~ "NH Other",
                          TRUE ~ NA_character_)) %>% 
  dplyr::select(study_id, study, race, female)

search <- readRDS(paste0(path_diabetes_subphenotypes_youth_folder,"/working/search/search_etiologic.RDS")) %>% 
  dplyr::select(-race) %>% 
  rename(race = race_eth) %>% 
  mutate(study = "search") %>% 
  dplyr::select(study_id, study, race, female)

today <- readRDS(paste0(path_diabetes_subphenotypes_youth_folder,"/working/today/today_baseline.RDS")) %>% 
  dplyr::select(-race) %>% 
  mutate(study_id = as.numeric(study_id),
         study = "today") %>% 
  rename(race = race_eth) %>% 
  dplyr::select(study_id, study, race, female)


pooled_df <- bind_rows(
  accord_newdm,
  aric_total,
  cardia_total,
  dppos_total,
  hrs_total,
  jhs_total,
  la_newdm,
  mesa_total,
  nhanes_total,
  oneflorida,
  search,
  today
)




















