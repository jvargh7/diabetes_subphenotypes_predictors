rm(list = ls());gc();source(".Rprofile")
#  ,"hc","triceps","iliac","abdominal","medial"--> not there in jhs
anthro_vars <- c("sbp","dbp","height","weight","wc","bmi")
# "glucose2h","apo_a","apo_b","vldlc","ast","alt","uric_acid" --> not there in jhs
lab_vars <- c("hba1c","insulinf","glucosef","tgl","hdlc","ldlc","totalc","ratio_th",
              "serumcreatinine","urinecreatinine","urinealbumin","egfr")
med_vars <- c("med_dm_use")

jhs_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/jhs_newdm.RDS")) 
jhs_baselinedm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/jhs_baseline_dm.RDS")) 

jhs_analysis <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/jhspre01_jhs_analysis.RDS")) %>% 
  dplyr::filter(aric == 0) %>% 
  dplyr::filter(!study_id %in% jhs_baselinedm$study_id) %>%
  arrange(study_id,visit) %>% 
  mutate(race_eth = "NH Black",
         female = case_when(female == "Female" ~ 1,
                            female == "Male" ~ 0,
                            TRUE ~ NA_real_)) %>% 
  mutate(ratio_th=tgl/hdlc,
         rece_clean = race_eth)

# jhs_longitudinal --------------

# Haven't filtered out observations where age >= dmagediag

jhs_longitudinal = jhs_analysis %>% 
  dplyr::select(-dmagediag) %>% 
  # Bringing the updated dmagediag from aric_events
  left_join(jhs_newdm %>% 
              dplyr::select(study_id,dmagediag),
            by=c("study_id")) %>% 
  mutate(med_dm_use = case_when(
    (dmmeds==1 | dmmedsoral==1 | dmmedsins==1) ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(
    available_labs = rowSums(!is.na(.[,lab_vars])),
    available_anthro = rowSums(!is.na(.[,anthro_vars]))) %>% 
  dplyr::select(study_id,visit,age,dmagediag,female,race_eth,rece_clean,
                alcohol,smoking,
                available_labs,available_anthro,one_of(anthro_vars),one_of(lab_vars),one_of(med_vars))



saveRDS(jhs_longitudinal,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01e_jhs.RDS"))
jhs_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01e_jhs.RDS"))


