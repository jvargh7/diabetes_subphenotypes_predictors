rm(list = ls());gc();source(".Rprofile")
# ,"hc","triceps","iliac","abdominal","medial" --> not there in mesa
anthro_vars <- c("sbp","dbp","height","wc","bmi")
# "glucose2h","vldlc","ast","alt","apo_a","apo_b","uric_acid" --> not there in mesa
lab_vars <- c("hba1c","insulinf2","glucosef2","tgl","hdlc","ldlc",
              "serumcreatinine","urinecreatinine","egfr")

mesa_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/processed/final_dataset_temp.RDS")) %>% 
  dplyr::filter(study == "mesa") %>% 
  mutate(original_study_id = as.numeric(original_study_id),
         dmagediag_cluster = round(dmagediag,2))

mesa_dat_all <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/mesa_dat_all.RDS")) %>% 
  mutate(history_dm = case_when(exam > 1 ~ NA_real_,
                                dia_sr == 1 ~ 1,
                                TRUE ~ 0),
         
         diab_e1 = case_when(exam > 1 ~ NA_real_,
                             dia_sr == 1  ~ 1,
                             TRUE ~ 0),
         
         diab_new_e1 = case_when(exam > 1 ~ NA_real_,
                             dia_ada == 2 | (glucosef >= 126 & !is.na(glucosef)) ~ 1,
                             TRUE ~ 0),
         
         diab_new_e2 = case_when(exam != 2 ~ NA_real_,
                                 dia_ada == 2 | dia_urine == 1 | (glucosef >= 126 & !is.na(glucosef)) | (hba1c >= 6.5 & !is.na(hba1c)) ~ 1,
                                 TRUE ~ 0),
         diab_new_e3 = case_when(exam != 3 ~ NA_real_,
                                 dia_ada == 2 | (glucosef >= 126 & !is.na(glucosef)) ~ 1,
                                 TRUE ~ 0),
         
         diab_new_e4 = case_when(exam !=4 ~ NA_real_,
                                 dia_ada == 2 | (glucosef >= 126 & !is.na(glucosef)) ~ 1,
                                 TRUE ~ 0),
         
         diab_new_e5 = case_when(exam !=5 ~ NA_real_,
                                 dia_ada == 2 | dia_dx_lv ==1 | (glucosef >= 126 & !is.na(glucosef))| (hba1c >= 6.5 & !is.na(hba1c)) ~ 1,
                                 TRUE ~ 0)
         ) %>% 
  arrange(study_id,exam) %>% 
  dplyr::filter(!is.na(age)) %>% 
  group_by(study_id) %>% 
  mutate(across(matches("(history_dm|diab_e1|diab_new_e)"),function(x) zoo::na.locf(x,na.rm=FALSE))) %>% 
  mutate(earliest_age = min(age,na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(dmdiagexam = case_when(
                                 diab_e1 == 1 | diab_new_e1 == 1 ~ 1,
                                 diab_new_e2 == 1 ~ 2,
                                 diab_new_e3 == 1 ~ 3,
                                 diab_new_e4 == 1 ~ 4,
                                 diab_new_e5 == 1 ~ 5,
                                 TRUE ~ NA_real_)) %>% 
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
                ))


mesa_events = mesa_dat_all %>% 
  dplyr::filter(exam == dmdiagexam) %>% 
  dplyr::select(study_id,earliest_age,history_dm,age,exam) %>% 
  mutate(dmagediag = case_when(history_dm == 1 ~ -1,
                               TRUE ~ age)) 

with(mesa_events,table((age-dmagediag) <= 1))

write_csv(mesa_events,paste0(path_diabetes_subphenotypes_adults_folder,"/working/qc/dsppre01f_mesa_events from diabetes_subphenotypes_predictors.csv"))


# mesa_longitudinal --------------

mesa_longitudinal = mesa_dat_all %>% 
  # Bringing the updated dmagediag from mesa_events
  left_join(mesa_events %>% 
              dplyr::select(-age,-earliest_age,-history_dm) %>% 
              dplyr::select(-exam),
            by=c("study_id")) %>% 
  mutate(
    available_labs = rowSums(!is.na(.[,lab_vars])),
    available_anthro = rowSums(!is.na(.[,anthro_vars])))

# Before dmagediag
before_dmagediag = join_by(study_id == original_study_id,
                           age < dmagediag_cluster)


mesa_newdm_not_in_longitudinal = mesa_newdm %>% 
  anti_join(mesa_longitudinal %>% 
              dplyr::select(study_id),
            by=c("original_study_id" = "study_id"))

mesa_longitudinal_newdm = mesa_longitudinal %>% 
  inner_join(mesa_newdm %>% 
               dplyr::select(original_study_id,dmagediag) %>% 
               rename(dmagediag_cluster = dmagediag),
             by = before_dmagediag) %>% 
  # !is.na(hba1c) -- not collected for most participants
  dplyr::filter(!is.na(hba1c)|!is.na(glucosef2),!is.na(bmi)) %>% 
  mutate(diff_dmagediag = dmagediag - age)

mesa_longitudinal_neverdm = mesa_longitudinal %>% 
  dplyr::filter(is.na(dmagediag)) %>% 
  # !is.na(hba1c) -- not collected for most participants
  dplyr::filter(!is.na(hba1c)|!is.na(glucosef2),!is.na(bmi)) %>% 
  # Among study waves where fasting glucose, HbA1c and BMI are measured, get the penultimate (second-to-last) wave
  group_by(study_id) %>% 
  mutate(wave = 1:n()) %>% 
  mutate(diff_next = case_when(wave == n() ~ NA_real_,
                               TRUE ~ dplyr::lead(age,1) - age)) %>% 
  ungroup() 

mesa_selected = bind_rows(
  mesa_longitudinal_newdm %>% 
    # dplyr::filter(diff_dmagediag >= 9/12,diff_dmagediag <= 15/12) %>% 
    mutate(type = "pre_1y_newdm") %>% 
    rename(diff_age = diff_dmagediag),
  mesa_longitudinal_neverdm %>% 
    # dplyr::filter(!is.na(diff_next),diff_next >= 9/12,diff_next <= 15/12) %>% 
    mutate(type = "pre_1y_neverdm") %>% 
    rename(diff_age = diff_next)
  
) %>% 
  dplyr::select(type,study_id,diff_age,age,available_labs,available_anthro,one_of(anthro_vars),one_of(lab_vars)) %>% 
  arrange(type,study_id,age)

saveRDS(mesa_selected,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01f_mesa.RDS"))
mesa_selected <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01f_mesa.RDS"))
