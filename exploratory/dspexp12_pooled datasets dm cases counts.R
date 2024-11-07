rm(list=ls());gc();source(".Rprofile")

accord_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/accord_newdm.RDS")) %>% 
  dplyr::select(study_id, race_eth, female) %>% 
  mutate(study = "accord",
         ethnicity = case_when(race_eth == "Hispanic" ~ "hispanic",
                               race_eth == "NH Black" | race_eth == "NH White" | race_eth == "NH Other" ~ "non-hispanic",
                               TRUE ~ "unknown"))

aric_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp02_aric new and no dm.RDS")) %>% 
  dplyr::select(study_id, race, female) %>% 
  mutate(study_id = as.numeric(study_id),
         study = "aric",
         ethnicity = case_when(race == "NH Black" | race == "NH White" ~ "non-hispanic",
                               TRUE ~ "unknown"),
         race_eth = race) 

cardia_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp03_cardia new and no dm.RDS")) %>% 
  dplyr::select(study_id, race, female) %>% 
  mutate(study = "cardia",
         ethnicity = case_when(race == "NH Black" | race == "NH White" ~ "non-hispanic",
                               TRUE ~ "unknown"),
         race_eth = race)
  

dppos_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp04_dpp new and no dm.RDS")) %>% 
  dplyr::select(study_id, race_eth, sex) %>% 
  mutate(female = case_when(sex == 1 ~ 0,
                            sex == 2 ~ 1,
                            TRUE ~ NA_real_),
         study = "dpp/dppos",
         ethnicity = case_when(race_eth == "Hispanic" ~ "hispanic",
                               race_eth == "NH Black" | race_eth == "NH White" | race_eth == "NH Other" ~ "non-hispanic",
                               TRUE ~ "unknown"))
  

hrs_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp05_hrs new and no dm.RDS")) %>%
  dplyr::select(hhidpn, race, gender) %>% 
  mutate(study_id = as.numeric(hhidpn),
         study = "hrs",
         female = case_when(gender == "Female" ~ 1,
                            gender == "Male" ~ 0,
                            TRUE ~ NA_real_),
         ethnicity = case_when(race == 1 | race == 2 | race == 3 ~ "non-hispanic",
                               TRUE ~ "unknown"),
         race_eth = case_when(race == 1 ~ "NH White",
                          race == 2 ~ "NH Black",
                          race == 3 ~ "NH Other",
                          TRUE ~ "Unknown"))

jhs_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp06_jhs new and no dm.RDS")) %>% 
  dplyr::select(study_id, race_eth, female) %>% 
  mutate(study = "jhs",
         ethnicity = case_when(race_eth == "NH Black" ~ "non-hispanic",
                               TRUE ~ "unknown")) 

la_newdm <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/la_newdm.RDS")) %>% 
  dplyr::select(study_id, race_eth, female) %>%
  mutate(study = "look ahead",
         ethnicity = case_when(race_eth == "Hispanic" ~ "hispanic",
                               race_eth == "NH Black" | race_eth == "NH White" | race_eth == "NH Other" ~ "non-hispanic",
                               TRUE ~ "unknown"))

mesa_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp08_mesa new and no dm.RDS")) %>% 
  dplyr::select(study_id, race, female) %>% 
  mutate(ethnicity = case_when(race == "Hispanic" ~ "hispanic",
                               race == "NH Black" | race == "NH White" | race == "Other" ~ "non-hispanic",
                               TRUE ~ "unknown"),
         race_eth = case_when(race == "Other" ~ "NH Other",
                          TRUE ~ race),
         study = "mesa")
  

nhanes_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp09_nhanes new and no dm.RDS")) %>% 
  dplyr::select(respondentid, race, gender) %>% 
  mutate(study_id = respondentid,
         study = "nhanes",
         ethnicity = case_when(race == 1 | race == 2 ~ "hispanic",
                               race == 3 | race == 4 | race == 6 | race == 7 ~ "non-hispanic",
                               TRUE ~ "unknown"),
         race_eth = case_when(race == 3 ~ "NH White",
                          race == 4 ~ "NH Black",
                          race == 6 ~ "Asian",
                          race == 7 ~ "NH Other",
                          race == 1 | race == 2 ~ "Hispanic",
                          TRUE ~ "Unknown"),
         female = case_when(gender == 1 ~ 0,
                            gender == 2 ~ 1,
                            TRUE ~ NA_real_)) 
  

oneflorida <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp10_oneflorida new and no dm.RDS")) %>% 
  dplyr::select(ID, nhwhite, nhblack, hispanic, nhother, female) %>% 
  mutate(study_id = as.numeric(ID),
         study = "oneflorida",
         race_eth = case_when(nhwhite == 1 & nhblack == 0 & hispanic == 0 & nhother == 0 ~ "NH White",
                              nhwhite == 0 & nhblack == 1 & hispanic == 0 & nhother == 0 ~ "NH Black",
                              nhwhite == 0 & nhblack == 0 & hispanic == 1 & nhother == 0 ~ "Hispanic",
                              nhwhite == 0 & nhblack == 0 & hispanic == 0 & nhother == 1 ~ "NH Other",
                              TRUE ~ "Unknown"),
         ethnicity = case_when(race_eth == "Hispanic" ~ "hispanic",
                               race_eth == "NH Black" | race_eth == "NH White" | race_eth == "NH Other" ~ "non-hispanic",
                               TRUE ~ "unknown")) 

search <- readRDS(paste0(path_diabetes_subphenotypes_youth_folder,"/working/search/search_etiologic.RDS")) %>% 
  dplyr::select(study_id, race_eth, ethnicity, female) %>% 
  mutate(study = "search",
         ethnicity = case_when(race_eth == "Hispanic" ~ "hispanic",
                               TRUE ~ ethnicity))
  

today <- readRDS(paste0(path_diabetes_subphenotypes_youth_folder,"/working/today/today_baseline.RDS")) %>% 
  dplyr::select(study_id, ethnicity, race_eth, female) %>% 
  mutate(study_id = as.numeric(study_id),
         study = "today") 


#-------------------------------------------------------------------------------------------------------------
dataset_list <- list(
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


prepare_data <- function(df) {
  df %>% 
    select(study_id, study, female, race_eth, ethnicity) %>%
    mutate(
      race_eth = case_when(is.na(race_eth) | race_eth == "Unknown" ~ "Unknown",
                           TRUE ~ race_eth),
      ethnicity = case_when(is.na(ethnicity) | ethnicity == "unknown" ~ "unknown",
                            TRUE ~ ethnicity),
      sex = case_when(female == 1 ~ "Female",
                      female == 0 ~ "Male",
                      TRUE ~ "Unknown")
    )
}

pooled_df <- bind_rows(lapply(dataset_list, prepare_data))

saveRDS(pooled_df, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp12_pooled dataset new and no dm.RDS"))

#-------------------------------------------------------------------------------------------------------------

df <- pooled_df %>% 
  mutate(eth_sex = case_when(sex == "Female" & ethnicity == "hispanic" ~ "hispanic female",
                             sex == "Male" & ethnicity == "hispanic" ~ "hispanic Male",
                             sex == "Unknown" & ethnicity == "hispanic" ~ "hispanic unknown",
                             sex == "Female" & ethnicity == "non-hispanic" ~ "non-hispanic female",
                             sex == "Male" & ethnicity == "non-hispanic" ~ "non-hispanic Male",
                             sex == "Unknown" & ethnicity == "non-hispanic" ~ "non-hispanic unknown",
                             sex == "Female" & ethnicity == "unknown" ~ "unknown female",
                             sex == "Male" & ethnicity == "unknown" ~ "unknown male",
                             TRUE ~ "unknown unknown"))

df %>%
  group_by(eth_sex) %>%
  summarise(
    unique_study_ids = n_distinct(study_id) 
  )

df %>%
  group_by(race_eth) %>%
  summarise(
    unique_study_ids = n_distinct(study_id) 
  )


df1 <- df %>% 
  dplyr::filter(ethnicity == "unknown" & sex == "Male")

df1 %>%
  group_by(race_eth) %>%
  summarise(
    unique_study_ids = n_distinct(study_id) 
  )









