rm(list = ls());gc();source(".Rprofile")
# ,"hc","triceps","iliac","abdominal","medial" --> not there in mesa
anthro_vars <- c("sbp","dbp","height","weight","wc","bmi")
# "glucose2h","vldlc","ast","alt","apo_a","apo_b","uric_acid" --> not there in mesa
lab_vars <- c("hba1c","insulinf","glucosef","tgl","hdlc","ldlc","totalc",
              "serumcreatinine","urinealbumin","urinecreatinine","uacr","egfr")
med_vars <- c("med_bp_use","med_chol_use","med_dm_use")


mesa_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/mesa_newdm.RDS")) 
mesa_baselinedm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/mesa_baseline_dm.RDS")) 

mesa_dat_all <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/mesa_dat_all.RDS")) %>% 
  dplyr::filter(!study_id %in% mesa_baselinedm$study_id) %>% 
  arrange(study_id,exam) %>% 
  dplyr::mutate(ratio_th=tgl/hdlc,
                glucosef2=glucosef*0.0555,
                insulinf2=insulinr*6,
                female = 1 - female,
                race_rev = case_when(
                  race == 1 ~ "White",
                  race == 3 ~ "AA",
                  race == 2 | 4 ~ "Other",
                  TRUE ~ NA_character_  
                ),
                race = case_when(
                  race == 1 ~ "NH White",
                  race == 3 ~ "NH Black",
                  race == 2 ~ "Other",
                  race == 4 ~ "Hispanic",
                  TRUE ~ NA_character_  
                )) %>% 
  dplyr::filter(!is.na(age)) %>% 
  mutate(race_clean = race,
         insulinf = insulinr)

mesa_longitudinal = mesa_dat_all %>% 
  arrange(study_id,exam) %>% 
  # Bringing the updated dmagediag from aric_events
  left_join(mesa_newdm %>% 
              dplyr::select(study_id,dmagediag),
            by=c("study_id")) %>% 
  mutate(
    med_bp_use = case_when(med_bp == 1 ~ 1,
                           TRUE ~ 0),
    med_chol_use = case_when(med_lipid == 1 ~ 1,
                           TRUE ~ 0),
    med_dm_use = case_when(
      (dia_med==2 | !is.na(dia_med_type) | dia_med_ins==1 | dia_med_ins_oh==1 | dia_med_ins_1st) ~ 1,
      TRUE ~ 0
    )) %>% 
  mutate(
    available_labs = rowSums(!is.na(.[,lab_vars])),
    available_anthro = rowSums(!is.na(.[,anthro_vars]))) %>% 
  dplyr::select(study_id,exam,age,dmagediag,female,race,race_clean,ethnicity,dmfamilyhistory,
                ratio_th,smk_cig,smk_pipe,smk_tob,
                available_labs,available_anthro,
                one_of(anthro_vars),one_of(lab_vars),one_of(med_vars))



saveRDS(mesa_longitudinal,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01f_mesa.RDS"))

mesa_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01f_mesa.RDS"))
write_csv(mesa_longitudinal,paste0(path_prediabetes_subphenotypes_folder,"/working/longitudinal/mesa.csv"))
