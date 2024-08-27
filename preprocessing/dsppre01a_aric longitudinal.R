rm(list = ls());gc();source(".Rprofile")
# ,"hc","triceps","iliac","abdominal","medial" --> not there in aric
anthro_vars <- c("sbp","dbp","height","wc","bmi")
# "vldlc","ast","alt" --> not there in aric
lab_vars <- c("hba1c","insulinf","glucosef","glucose2h","tgl","hdlc","ldlc",
              "serumcreatinine","urinecreatinine","egfr","apo_a","apo_b","uric_acid")

aric_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/processed/final_dataset_temp.RDS")) %>% 
  dplyr::filter(study == "aric") %>% 
  mutate(dmagediag = round(dmagediag,2))

aric_analysis <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/aric_analysis.RDS")) %>%
  mutate(across(contains("evr"),function(x) case_when(x=="N" ~ 0,
                                                      x == "Y" ~ 1,
                                                      TRUE ~ NA_real_))) %>% 
  rename(dmagediag_orig = dmagediag) %>% 
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
  mutate(diab_126_fast = case_when(
    diab_126_fast=="1"|diab_126_fast=="T"~1,
    diab_126_fast=="0"~0,
    TRUE ~ NA_real_),
    diab_evr = case_when(
      diab_evr=="Y"~1,
      diab_evr=="N"~0,
      TRUE ~ NA_real_),
    diab_v1 = case_when(
      (visit==1)&(diab_126_fast==1|diab_evr==1|diab_ind==1)~1,
      visit > 1 ~ NA_real_,
      TRUE~0)
  ) %>% 
  # group_by(study_id) %>% 
  # mutate(diab_v1 = zoo::na.locf(diab_v1)) %>% 
  # ungroup() %>% 
  # dplyr::filter(diab_v1==0) %>% 
  mutate(
    
    # VISIT 2 ------
    diab_doc = case_when(
      diab_doc=="Y"~1,
      diab_doc=="N"~0, #there are other letters with unknown meanings in v3,v4 and v5,code to NA for now 
      TRUE ~ NA_real_),
    diab_new_v2 = case_when(
      (visit==2)&(diab_126_fast==1|diab_doc==1|diab_ind==1) ~ 1,
      visit!=2~NA_real_,
      TRUE~0),
    # VISIT 3 -------
    diab_126 = case_when(
      diab_126=="1"|diab_126=="T"~1,
      diab_126=="0"~0,
      TRUE ~ NA_real_),
    diab_140 = case_when(
      diab_140=="1"|diab_140=="T"~1,
      diab_140=="0"~0,
      TRUE ~ NA_real_),
    diab_med_2w = case_when(
      diab_med_2w=="Y"~1,
      diab_med_2w=="N"~0,
      TRUE ~ NA_real_),
    diab_new_v3 = case_when(
      (visit==3)&(diab_126==1|diab_doc==1) ~ 1,
      visit!=3~NA_real_,
      TRUE~0),
   
    # VISIT 4 ---------
    diab_trt= case_when(
      diab_trt == "Y"~1,
      diab_trt=="N"~0,
      TRUE ~ NA_real_),
    diab_med_any = case_when(
      diab_med_any=="Y"~1,
      diab_med_any=="N"~0,
      TRUE ~ NA_real_),
    diab_new_v4 = case_when(
      (visit==4)&(diab_126==1|diab_doc==1|(glucosef >= 126 & !is.na(glucosef))|(glucose2h >= 200 & !is.na(glucose2h))) ~ 1,
      visit!=4~NA_real_,
      TRUE~0),
    
    # VISIT 5 -------
    diab_a1c65= case_when(
      diab_a1c65 == "1"|diab_a1c65 =="T"~1,
      diab_a1c65=="0"~0,
      TRUE ~ NA_real_),
    diab_med_4w= case_when(
      diab_med_4w == "1"|diab_med_4w =="T"~1,
      diab_med_4w=="0"~0,
      TRUE ~ NA_real_),
    diab_new_v5 = case_when(
      (visit==5)&(diab_126==1|diab_a1c65==1|diab_doc==1|(glucosef >= 126 & !is.na(glucosef))|(hba1c >= 6.5 & !is.na(hba1c))) ~ 1,
      visit!=5~NA_real_,
      TRUE~0),
    
    # VISIT 6 --------
    diab_new_v6 = case_when(
      (visit==6)&(diab_126==1|diab_a1c65==1|(glucosef >= 126 & !is.na(glucosef))) ~ 1,
      visit!=6~NA_real_,
      TRUE~0)) %>% 
  dplyr::filter(!is.na(age)) %>% 
  distinct(study_id,visit,bmi,glucosef,.keep_all =TRUE) %>% 
  arrange(study_id,visit) %>% 
  group_by(study_id) %>% 
  mutate(across(matches("(diab_v1|diab_new_v)"),function(x) zoo::na.locf(x,na.rm=FALSE))) %>% 
  mutate(earliest_age = min(age),
         dmagediag_orig = min(dmagediag_orig,na.rm=TRUE)) %>% 
  # mutate(across(matches("(diab_new_v)"),function(x) zoo::na.locf(x,fromLast = TRUE,na.rm=FALSE))) %>% 
  ungroup() %>% 
  mutate(dmagediag_orig = case_when(dmagediag_orig == Inf ~ NA_real_,
                                    TRUE ~ dmagediag_orig)) %>% 
  mutate(dmdiagvisit = case_when(dmagediag_orig < earliest_age ~ 1,
                                 diab_v1 == 1 ~ 1,
                                 diab_new_v2 == 1 ~ 2,
                                 diab_new_v3 == 1 ~ 3,
                                 diab_new_v4 == 1 ~ 4,
                                 diab_new_v5 == 1 ~ 5,
                                 diab_new_v6 == 1 ~ 6,
                                 TRUE ~ NA_real_))

aric_events = aric_analysis %>% 
  dplyr::filter(visit == dmdiagvisit) %>% 
  dplyr::select(study_id,dmagediag_orig,age,visit) %>% 
  mutate(dmagediag = case_when(age < dmagediag_orig ~ age,
                               !is.na(dmagediag_orig) ~ dmagediag_orig,
                               TRUE ~ age)) 
write_csv(aric_events,paste0(path_diabetes_subphenotypes_adults_folder,"/working/qc/dsppre01a_aric_events from diabetes_subphenotypes_predictors.csv"))

table(aric_events$visit > 1) # 3,621

aric_longitudinal = aric_analysis %>% 
  # Bringing the updated dmagediag from aric_events
  left_join(aric_events %>% 
              dplyr::select(-age,-dmagediag_orig) %>% 
              dplyr::select(-visit),
            by=c("study_id")) %>% 
  mutate(
         available_labs = rowSums(!is.na(.[,lab_vars])),
         available_anthro = rowSums(!is.na(.[,anthro_vars])))

# Before dmagediag
before_dmagediag = join_by(study_id == original_study_id,
                           age < dmagediag_cluster)


aric_newdm_not_in_longitudinal = aric_newdm %>% 
  anti_join(aric_longitudinal %>% 
              dplyr::select(study_id),
            by=c("original_study_id" = "study_id"))

aric_longitudinal_newdm = aric_longitudinal %>% 
  inner_join(aric_newdm %>% 
               dplyr::select(original_study_id,dmagediag) %>% 
               rename(dmagediag_cluster = dmagediag),
             by = before_dmagediag) %>% 
  # !is.na(hba1c) -- not collected for most participants
  dplyr::filter(!is.na(glucosef),!is.na(bmi)) %>% 
  mutate(diff_dmagediag = dmagediag - age)

aric_longitudinal_neverdm = aric_longitudinal %>% 
  dplyr::filter(is.na(dmagediag)) %>% 
  # !is.na(hba1c) -- not collected for most participants
  dplyr::filter(!is.na(glucosef),!is.na(bmi)) %>% 
  # Among study waves where fasting glucose, HbA1c and BMI are measured, get the penultimate (second-to-last) wave
  group_by(study_id) %>% 
  mutate(wave = 1:n()) %>% 
  mutate(diff_next = case_when(wave == n() ~ NA_real_,
                               TRUE ~ dplyr::lead(age,1) - age)) %>% 
  ungroup() 

aric_selected = bind_rows(
  aric_longitudinal_newdm %>% 
    # dplyr::filter(diff_dmagediag >= 9/12,diff_dmagediag <= 15/12) %>% 
    mutate(type = "pre_1y_newdm") %>% 
    rename(diff_age = diff_dmagediag),
  aric_longitudinal_neverdm %>% 
    # dplyr::filter(!is.na(diff_next),diff_next >= 9/12,diff_next <= 15/12) %>% 
    mutate(type = "pre_1y_neverdm") %>% 
    rename(diff_age = diff_next)
  
) %>% 
  dplyr::select(type,study_id,diff_age,age,available_labs,available_anthro,one_of(anthro_vars),one_of(lab_vars)) %>% 
  arrange(type,study_id,age)

saveRDS(aric_selected,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01a_aric.RDS"))
aric_selected <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01a_aric.RDS"))

# QC ----------

## Who are missing in aric_events? ----
aric_events_newdm = aric_events %>% 
  dplyr::filter(visit > 1)

# From final_dataset_temp.RDS....
missing_aric_events = aric_newdm %>% 
  anti_join(aric_events_newdm,
            by=c("original_study_id"="study_id"))

# Get rows for study_id that are missing in aric_events
missing_aric_longitudinal = aric_longitudinal %>% 
       inner_join(missing_aric_events %>% 
                    dplyr::select(original_study_id,dmagediag) %>% 
                    rename(dmagediag_cluster = dmagediag),
                  by=c("study_id"="original_study_id")) %>% 
  dplyr::select(study_id,dmagediag_cluster,visit,dmdiagvisit,dmagediag,dmagediag_orig,contains("age"),contains("diab"),"glucosef","hba1c")

write_csv(missing_aric_longitudinal,paste0(path_diabetes_subphenotypes_adults_folder,"/working/qc/dsppre01a_missing_aric_longitudinal from diabetes_subphenotypes_predictors.csv"))
