rm(list=ls());gc();source(".Rprofile")

# cohort: 
# 1. only race, in answers, no "NH"; if they have "Hispanic" in their answers, we can assume others are NH
# 2. race has "NH" information

# N = 601
# raw data: raceclass (White, Black, Hispanic, Other); no separate enthnicty
accord_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/accord_newdm.RDS")) %>% 
  dplyr::select(study_id, race_eth, female) %>% 
  mutate(study = "accord",
         ethnicity = case_when(race_eth == "Hispanic" ~ "hispanic",
                               race_eth == "NH Black" | race_eth == "NH White" | race_eth == "NH Other" ~ "non-hispanic",
                               TRUE ~ "unknown"),
         race_cleaned = race_eth)
# N = 13,817
# raw data: no ethnicity; race (W, B); should use race_rev
# NH: to be consistent
aric_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp02_aric new and no dm.RDS")) %>% 
  dplyr::select(study_id, race_rev, female) %>% 
  mutate(study_id = as.integer(sub("C", "", study_id)),
         study = "aric",
         ethnicity = "unknown",
         race_cleaned = case_when(race_rev == "White" ~ "White",
                                  race_rev == "AA" ~ "Black",
                                  TRUE ~ NA_character_)) 
# N = 5054
# raw data: codebook has ethnicity; in extracted raw data, race only has nh white (4) and nh black (5)
cardia_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp03_cardia new and no dm.RDS")) %>% 
  dplyr::select(study_id, race, female) %>% 
  mutate(study = "cardia",
         ethnicity = case_when(race == "NH Black" | race == "NH White" ~ "non-hispanic",
                               TRUE ~ "unknown"),
         race_cleaned = race)
  
# N = 3589
# raw data: no ethnicity; race_eth in raw data (white, black, other, hispanic, assume NH)
dppos_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp04_dpp new and no dm.RDS")) %>% 
  dplyr::select(study_id, race_eth, sex) %>% 
  mutate(female = case_when(sex == 1 ~ 0,
                            sex == 2 ~ 1,
                            TRUE ~ NA_real_),
         study = "dpp/dppos",
         ethnicity = case_when(race_eth == "Hispanic" ~ "hispanic",
                               race_eth == "NH Black" | race_eth == "NH White" | race_eth == "NH Other" ~ "non-hispanic",
                               TRUE ~ "unknown"),
         race_cleaned = race_eth)
  
# N = 35,973
hrs_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp05_hrs new and no dm.RDS")) %>%
  dplyr::select(hhidpn, race, ethnicity, gender) %>% 
  mutate(study_id = as.numeric(hhidpn),
         study = "hrs",
         female = case_when(gender == "Female" ~ 1,
                            gender == "Male" ~ 0,
                            TRUE ~ NA_real_),
         race_cleaned = race)
# N = 1817
# ask Jithin about race_eth, not in codebook now
# participants: African Americans
jhs_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp06_jhs new and no dm.RDS")) %>% 
  dplyr::select(study_id, female) %>% 
  mutate(study = "jhs",
         ethnicity = "unknown",
         race_cleaned = "Black") 
# N = 877
# black, white, other
la_newdm <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/la_newdm.RDS")) %>% 
  dplyr::select(study_id, race_eth, female) %>%
  mutate(study = "look ahead",
         ethnicity = case_when(race_eth == "Hispanic" ~ "hispanic",
                               race_eth == "NH Black" | race_eth == "NH White" | race_eth == "NH Other" ~ "non-hispanic",
                               TRUE ~ "unknown"),
         race_cleaned = race_eth)
# N = 6094
# raw data: race (1,2,3,4...); white, black, chinese, hispanic; have ethnicity
mesa_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp08_mesa new and no dm.RDS")) %>% 
  dplyr::select(study_id, race, ethnicity, female) %>% 
  mutate(study = "mesa",
         ethnicity = case_when(ethnicity == 1 ~ "hispanic",
                               ethnicity == 0 ~ "non-hispanic",
                               TRUE ~ "unknown"),
         race_cleaned = race)
  
# N = 96,109
nhanes_total <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp09_nhanes new and no dm.RDS")) %>% 
  dplyr::select(respondentid, race, gender) %>% 
  mutate(study_id = respondentid,
         study = "nhanes",
         ethnicity = case_when(race == 1 | race == 2 ~ "hispanic",
                               race == 3 | race == 4 | race == 6 | race == 7 ~ "non-hispanic",
                               TRUE ~ "unknown"),
         race_cleaned = case_when(race == 3 ~ "NH White",
                          race == 4 ~ "NH Black",
                          race == 6 ~ "Asian",
                          race == 7 ~ "NH Other",
                          race == 1 | race == 2 ~ "Hispanic",
                          TRUE ~ "Unknown"),
         female = case_when(gender == 1 ~ 0,
                            gender == 2 ~ 1,
                            TRUE ~ NA_real_)) 
  

# N = 329,545
oneflorida <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp10_oneflorida new and no dm.RDS")) %>% 
  dplyr::select(ID, nhwhite, nhblack, hispanic, nhother, female) %>% 
  group_by(ID) %>%
  mutate(study_id = cur_group_id()) %>%
  ungroup() %>% 
  mutate(study = "oneflorida",
         race_cleaned = case_when(nhwhite == 1 & nhblack == 0 & hispanic == 0 & nhother == 0 ~ "NH White",
                                  nhwhite == 0 & nhblack == 1 & hispanic == 0 & nhother == 0 ~ "NH Black",
                                  nhwhite == 0 & nhblack == 0 & hispanic == 1 & nhother == 0 ~ "Hispanic",
                                  nhwhite == 0 & nhblack == 0 & hispanic == 0 & nhother == 1 ~ "NH Other",
                                  TRUE ~ "Unknown"),
         ethnicity = case_when(hispanic == 1 ~ "hispanic",
                               hispanic == 0 ~ "non-hispanic",
                               TRUE ~ "unknown")) 
# N = 641
youth <- read.csv(paste0(path_diabetes_subphenotypes_youth_folder,"/working/cleaned/etiologic/setdy01a_analytic sample.csv")) %>% 
  dplyr::select(study_id, study, race, race_eth, ethnicity, female) %>% 
  mutate(study_id = as.numeric(gsub("-", "", study_id)),
         race_cleaned = race)


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
  youth
)


prepare_data <- function(df) {
  df %>% 
    select(study_id, study, female, race_cleaned, ethnicity) %>%
    mutate(
      race_cleaned = case_when(is.na(race_cleaned) | race_cleaned == "Unknown" ~ "Unknown",
                       TRUE ~ race_cleaned),
      ethnicity = case_when(is.na(ethnicity) | ethnicity == "unknown" ~ "unknown",
                            TRUE ~ ethnicity),
      sex = case_when(female == 1 ~ "female",
                      female == 0 ~ "male",
                      TRUE ~ "unknown")
    )
}
# N = 494,117
pooled_df <- bind_rows(lapply(dataset_list, prepare_data)) %>% 
  mutate(new_id = paste(study, study_id, sep = "_"))

saveRDS(pooled_df, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp12_pooled dataset new and no dm.RDS"))

#-------------------------------------------------------------------------------------------------------------

df <- pooled_df %>% 
  mutate(eth_sex = case_when(sex == "female" & ethnicity == "hispanic" ~ "hispanic female",
                             sex == "male" & ethnicity == "hispanic" ~ "hispanic male",
                             sex == "unknown" & ethnicity == "hispanic" ~ "hispanic unknown",
                             sex == "female" & ethnicity == "non-hispanic" ~ "non-hispanic female",
                             sex == "male" & ethnicity == "non-hispanic" ~ "non-hispanic male",
                             sex == "unknown" & ethnicity == "non-hispanic" ~ "non-hispanic unknown",
                             sex == "female" & ethnicity == "unknown" ~ "unknown female",
                             sex == "male" & ethnicity == "unknown" ~ "unknown male",
                             TRUE ~ "unknown unknown"))

df %>%
  group_by(eth_sex) %>%
  summarise(
    unique_ids = n_distinct(new_id) 
  )

df %>%
  group_by(race_cleaned) %>%
  summarise(
    unique_ids = n_distinct(new_id) 
  )


df1 <- df %>% 
  dplyr::filter(ethnicity == "unknown" & sex == "female")

df1 %>%
  group_by(race_cleaned) %>%
  summarise(
    unique_ids = n_distinct(new_id) 
  )









