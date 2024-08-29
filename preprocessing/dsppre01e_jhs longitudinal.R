rm(list = ls());gc();source(".Rprofile")
#  ,"hc","triceps","iliac","abdominal","medial"--> not there in jhs
anthro_vars <- c("sbp","dbp","height","wc","bmi")
# "glucose2h","apo_a","apo_b","vldlc","ast","alt","uric_acid" --> not there in jhs
lab_vars <- c("hba1c","insulinf","glucosef","tgl","hdlc","ldlc",
              "serumcreatinine","urinecreatinine","egfr")

jhs_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/processed/final_dataset_temp.RDS")) %>% 
  dplyr::filter(study == "jhs") %>% 
  mutate(dmagediag_cluster = round(dmagediag,2),
         original_study_id = as.numeric(original_study_id))



jhs_analysis <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/jhspre01_jhs_analysis.RDS")) %>% 
  dplyr::filter(aric == 0) %>%
  arrange(study_id,visit)  %>% 
  group_by(study_id) %>% 
  # Self-reported diabetes carried forward as 1 if 1
  mutate(diabetes = zoo::na.locf(diabetes,na.rm=FALSE),
         dmagediag = zoo::na.locf(dmagediag,na.rm=FALSE)) %>% 
  mutate(dmagediag_ever = case_when(min(dmagediag,na.rm=TRUE) == Inf ~ NA_real_,
                                    TRUE ~ min(dmagediag,na.rm=TRUE))) %>% 
  ungroup() %>% 
  mutate(diab_v1 = case_when(visit > 1 ~ NA_real_,
                             dmagediag_ever < age ~ 1,
                             diabetes == 1 ~ 1,
                             TRUE ~ 0),
         diab_new_v1 = case_when(visit > 1 ~ NA_real_,
                             dmagediag_ever == age ~ 1,
                             (glucosef >=126 | hba1c >= 6.5) ~ 1,
                             TRUE ~ 0),
         diab_new_v2 = case_when(visit !=2 ~ NA_real_,
                                 dmagediag_ever <= age ~ 1,
                                 diabetes == 1 | (glucosef >=126 | hba1c >= 6.5) ~ 1,
                                 TRUE ~ 0
                                 ),
         diab_new_v3 = case_when(visit !=3 ~ NA_real_,
                                 dmagediag_ever <= age ~ 1,
                                 diabetes == 1 | (glucosef >=126 | hba1c >= 6.5) ~ 1,
                                 TRUE ~ 0)
         ) %>% 
  group_by(study_id) %>% 
  mutate(across(matches("(diab_v1|diab_new_v)"),function(x) zoo::na.locf(x,na.rm=FALSE))) %>% 
  mutate(earliest_age = min(age,na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(dmdiagvisit = case_when(
    diab_v1 == 1 | diab_new_v1 == 1 ~ 1,
    diab_new_v2 == 1 ~ 2,
    diab_new_v3 == 1 ~ 3,
    TRUE ~ NA_real_)) %>% 
  
  mutate(
    race_eth = "NH Black",
    female = case_when(female == "Female" ~ 1,
                       female == "Male" ~ 0,
                       TRUE ~ NA_real_)
    
  ) 

write_csv(jhs_analysis,paste0(path_diabetes_subphenotypes_adults_folder,"/working/qc/dsppre01e_jhs_analysis from diabetes_subphenotypes_predictors.csv"))

jhs_events = jhs_analysis %>% 
  dplyr::filter(visit == dmdiagvisit) %>% 
  dplyr::select(study_id,dmagediag_ever,earliest_age,age,visit) %>% 
  mutate(dmagediag = case_when(age < dmagediag_ever ~ age,
                               !is.na(dmagediag_ever) ~ dmagediag_ever,
                               TRUE ~ age)) 

with(jhs_events,table((age-dmagediag) <= 1))

write_csv(jhs_events,paste0(path_diabetes_subphenotypes_adults_folder,"/working/qc/dsppre01e_jhs_events from diabetes_subphenotypes_predictors.csv"))


# jhs_longitudinal --------------

jhs_longitudinal = jhs_analysis %>% 
  dplyr::select(-dmagediag) %>% 
  # Bringing the updated dmagediag from jhs_events
  left_join(jhs_events %>% 
              dplyr::select(-age,-earliest_age,-dmagediag_ever) %>% 
              dplyr::select(-visit),
            by=c("study_id")) %>% 
  mutate(
    available_labs = rowSums(!is.na(.[,lab_vars])),
    available_anthro = rowSums(!is.na(.[,anthro_vars])))

# Before dmagediag
before_dmagediag = join_by(study_id == original_study_id,
                           age < dmagediag_cluster)


jhs_newdm_not_in_longitudinal = jhs_newdm %>% 
  anti_join(jhs_longitudinal %>% 
              dplyr::select(study_id),
            by=c("original_study_id" = "study_id"))

jhs_longitudinal_newdm = jhs_longitudinal %>% 
  inner_join(jhs_newdm %>% 
               dplyr::select(original_study_id,dmagediag) %>% 
               rename(dmagediag_cluster = dmagediag),
             by = before_dmagediag) %>% 
  # !is.na(hba1c) -- not collected for most participants
  dplyr::filter(!is.na(hba1c)|!is.na(glucosef),!is.na(bmi)) %>% 
  mutate(diff_dmagediag = dmagediag - age)

jhs_longitudinal_neverdm = jhs_longitudinal %>% 
  dplyr::filter(is.na(dmagediag)) %>% 
  # !is.na(hba1c) -- not collected for most participants
  dplyr::filter(!is.na(hba1c)|!is.na(glucosef),!is.na(bmi)) %>% 
  # Among study waves where fasting glucose, HbA1c and BMI are measured, get the penultimate (second-to-last) wave
  group_by(study_id) %>% 
  mutate(wave = 1:n()) %>% 
  mutate(diff_next = case_when(wave == n() ~ NA_real_,
                               TRUE ~ dplyr::lead(age,1) - age)) %>% 
  ungroup() 

jhs_selected = bind_rows(
  jhs_longitudinal_newdm %>% 
    # dplyr::filter(diff_dmagediag >= 9/12,diff_dmagediag <= 15/12) %>% 
    mutate(type = "pre_1y_newdm") %>% 
    rename(diff_age = diff_dmagediag),
  jhs_longitudinal_neverdm %>% 
    # dplyr::filter(!is.na(diff_next),diff_next >= 9/12,diff_next <= 15/12) %>% 
    mutate(type = "pre_1y_neverdm") %>% 
    rename(diff_age = diff_next)
  
) %>% 
  dplyr::select(type,study_id,diff_age,age,available_labs,available_anthro,one_of(anthro_vars),one_of(lab_vars)) %>% 
  arrange(type,study_id,age)

saveRDS(jhs_selected,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01e_jhs.RDS"))
jhs_selected <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01e_jhs.RDS"))


