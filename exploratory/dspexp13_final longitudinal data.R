rm(list=ls());gc();source(".Rprofile")

library(stringr)

# Longitudinal dataset: new + no T2D; 8 cohorts; historical + follow-up

accord_long = readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp01_accord new dm.RDS")) %>% 
  mutate(study = "accord") %>% 
  mutate(wave = visit,
         rece_clean = race_eth)
 

la_long <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp07_look ahead new dm.RDS")) %>% 
  mutate(study = "look ahead") %>% 
  mutate(wave = visit) 
 
aric_long <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp02_aric new and no dm.RDS")) %>% 
  mutate(study_id = as.numeric(str_replace(study_id,"C","")),
         study = "aric") %>% 
  mutate(wave = as.character(visit))

cardia_long <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp03_cardia new and no dm.RDS")) %>% 
  mutate(study = "cardia") %>% 
  mutate(race_clean = race) %>% 
  mutate(wave = as.character(year)) 

jhs_long <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp06_jhs new and no dm.RDS")) %>% 
  mutate(study = "jhs") %>% 
  mutate(wave = as.character(visit)) %>% 
  mutate(uacr = urinealbumin/urinecreatinine)

mesa_long <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp08_mesa new and no dm.RDS")) %>% 
  mutate(study = "mesa") %>% 
  mutate(wave = as.character(exam)) 

dppos_long <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp04_dpp new and no dm.RDS")) %>% 
  mutate(study = "dppos") %>% 
  mutate(ratio_th=tgl/hdlc) %>% 
  mutate(wave = as.character(wave))

hrs_long <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp05_hrs new and no dm.RDS")) %>% 
  mutate(study = "hrs") %>% 
  mutate(wave = as.character(wave)) %>% 
  mutate(bmi = weight/(height^2),
         height = height*100)



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

var_list <- c("study_id", "study", "wave", "female", "race_clean", "age", "dmagediag", "dmduration", "dmfamilyhistory",
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
         race_clean = case_when(new_id == "mesa_1076" | new_id == "mesa_1359" | new_id == "mesa_2614" ~ "Other",
                                  TRUE ~ race_clean))

saveRDS(pooled_df, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp13_final longitudinal data.RDS"))


