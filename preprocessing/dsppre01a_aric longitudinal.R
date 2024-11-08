rm(list = ls());gc();source(".Rprofile")
# ,"hc","triceps","iliac","abdominal","medial" --> not there in aric
anthro_vars <- c("sbp","dbp","height","wc","bmi")
# "vldlc","ast","alt" --> not there in aric
lab_vars <- c("hba1c","insulinf","glucosef","glucose2h","tgl","hdlc","ldlc",
              "serumcreatinine","urinecreatinine","egfr","apo_a","apo_b","uric_acid")

aric_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/aric_newdm.RDS")) 
aric_baselinedm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/aric_baseline_dm.RDS")) 


aric_analysis <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/aric_analysis.RDS")) %>%
  dplyr::filter(!study_id %in% aric_baselinedm$study_id) %>%
  dplyr::filter(!is.na(age)) %>% 
  arrange(study_id,visit) %>% 
  dplyr::mutate(ratio_th=tgl/hdlc,
                glucosef2=glucosef*0.0555,
                insulinf2=insulinf*6,
                urinealbumin = urinealbumin/10,
                race_rev = case_when(
                  race == "W" ~ "White",
                  race == "B" ~ "AA",
                  TRUE ~ NA_character_  
                ),
                female = case_when(
                  female == "M" ~ 0,  
                  female == "F" ~ 1,  
                  TRUE ~ NA_integer_
                ),
                race = case_when(
                  race == "W" ~ "NH White",
                  race == "B" ~ "NH Black",
                  TRUE ~ NA_character_  
                ),
  ) %>% 
  group_by(study_id) %>% 
  mutate(across(one_of("female","race"),~zoo::na.locf(.))) %>% 
  ungroup()

# aric_longitudinal --------------

# Haven't filtered out observations where age >= dmagediag

aric_longitudinal = aric_analysis %>% 
  dplyr::select(-dmagediag) %>% 
  # Bringing the updated dmagediag from aric_events
  left_join(aric_newdm %>% 
              dplyr::select(study_id,dmagediag),
            by=c("study_id")) %>% 
  mutate(
         available_labs = rowSums(!is.na(.[,lab_vars])),
         available_anthro = rowSums(!is.na(.[,anthro_vars]))) %>% 
  dplyr::select(study_id,visit,age,dmagediag,female,available_labs,available_anthro,one_of(anthro_vars),one_of(lab_vars),race_rev,race)


saveRDS(aric_longitudinal,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01a_aric.RDS"))

# Writing to prediabetes phenotypes folder -------
aric_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01a_aric.RDS"))
write_csv(aric_longitudinal,paste0(path_prediabetes_subphenotypes_folder,"/working/longitudinal/aric.csv"))

