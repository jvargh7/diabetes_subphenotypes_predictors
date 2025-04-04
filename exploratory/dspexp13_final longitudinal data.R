rm(list=ls());gc();source(".Rprofile")

# Longitudinal dataset: new + no T2D; 8 cohorts; historical + follow-up

accord_long = readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp01_accord new dm.RDS")) %>% 
  mutate(study = "accord") %>% 
  mutate(bmi = weight/((height/100)^2),
         ratio_th=tgl/hdlc) %>% 
  rename(age = bsage,
         wave = visit)

la_long <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp07_look ahead new dm.RDS")) %>% 
  mutate(study = "look ahead") %>% 
  mutate(ratio_th=tgl/hdlc) %>% 
  rename(age = bsage,
         wave = visit) 
 
aric_long <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp02_aric new and no dm.RDS")) %>% 
  mutate(study_id = as.numeric(str_replace(study_id,"C","")),
         study = "aric") %>% 
  mutate(race_eth = case_when(
    race_rev == "White" ~ "White",
    race_rev == "AA" ~ "Black",
    TRUE ~ NA_character_
  ),
  smoking = case_when(
    smk_cur %in% c("1", "T", "Y") ~ "Present",
    smk_evr %in% c("1", "T", "Y") ~ "Past",
    (smk_evr %in% c("0", "N")) & (smk_cur %in% c("0", "N")) ~ "Never",
    TRUE ~ "Missing"
  )) %>%
  mutate(wave = as.character(visit)) %>% 
  select(-race,-race_rev) 

cardia_long <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp03_cardia new and no dm.RDS")) %>% 
  mutate(study = "cardia") %>% 
  mutate(race_eth = case_when(
    race_rev == "White" ~ "White",
    race_rev == "AA" ~ "Black",
    TRUE ~ NA_character_
  )) %>% 
  mutate(wave = as.character(year)) %>% 
  select(-race_rev)

jhs_long <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp06_jhs new and no dm.RDS")) %>% 
  mutate(study = "jhs") %>% 
  mutate(ratio_th=tgl/hdlc) %>% 
  mutate(wave = as.character(visit))

mesa_long <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp08_mesa new and no dm.RDS")) %>% 
  mutate(study = "mesa") %>% 
  mutate(insulinf = insulinf2/6,
         glucosef = glucosef2/0.0555) %>% 
  mutate(wave = as.character(exam)) %>% 
  rename(race_eth = race)

dppos_long <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp04_dpp new and no dm.RDS")) %>% 
  mutate(study = "dppos") %>% 
  mutate(ratio_th=tgl/hdlc) %>% 
  mutate(wave = as.character(wave))

hrs_long <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp05_hrs new and no dm.RDS")) %>% 
  mutate(study = "hrs") %>% 
  mutate(wave = as.character(wave))



#-------------------------------------------------------------------------------------------------------------
dataset_list <- list(
  accord_long,
  aric_long,
  cardia_long,
  dppos_long,
  hrs_long,
  jhs_long,
  la_long,
  mesa_long
)

var_list <- c("study_id", "study", "wave", "female", "race_eth", "age", "dmagediag", "dmduration", "dmfamilyhistory",
              "insulinf", "glucosef", "glucose2h", "med_dm_use", "med_bp_use", "med_chol_use",
              "height", "weight", "bmi", "wc", "sbp", "dbp", "hba1c", "totalc", "ldlc", "hdlc", "vldlc", "tgl", "wc", "ratio_th",
              "serumcreatinine", "urinealbumin", "urinecreatinine", "uacr", "egfr", "alt", "apo_a", "apo_b",
              "available_labs", "available_anthro")

prepare_data <- function(df) {
  df %>% 
    select(any_of(var_list)) %>% 
    mutate(sex = case_when(female == 1 ~ "female", 
                           female == 0 ~ "male",
                           TRUE ~ "unknown")
    ) 
}
# N = 494,117
pooled_df <- bind_rows(lapply(dataset_list, prepare_data)) %>% 
  mutate(new_id = paste(study, study_id, sep = "_"), 
         # assign "Other race" to 3 ppl with multiple race_eth
         race_eth = case_when(new_id == "mesa_1076" | new_id == "mesa_1359" | new_id == "mesa_2614" ~ "Other",
                                  TRUE ~ race_eth))

saveRDS(pooled_df, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp13_final longitudinal data.RDS"))
