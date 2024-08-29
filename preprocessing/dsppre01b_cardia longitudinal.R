rm(list = ls());gc();source(".Rprofile")
# "hc","triceps","iliac","abdominal","medial" --> not there in cardia
anthro_vars <- c("sbp","dbp","height","wc","bmi")
#  --> not there in aric
# ,"ast","alt","uric_acid" --> not there in cardia
lab_vars <- c("hba1c","insulinf","glucosef","glucose2h","tgl","hdlc","ldlc",
              "serumcreatinine","urinecreatinine","egfr","apo_a","apo_b","vldlc")

cardia_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/processed/final_dataset_temp.RDS")) %>% 
  dplyr::filter(study == "cardia") %>% 
  mutate(dmagediag_cluster = round(dmagediag,2),
         original_study_id = as.numeric(original_study_id))

cardia_dat_all = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/cardia_dat_all.RDS")) %>% 
  arrange(study_id,year) %>% 
  mutate(diab_y0 = case_when(year > 0 ~ NA_real_,
                             glucosef>=126&!is.na(glucosef) ~ 1,
                             TRUE ~ 0),
         diab_new_y2 = case_when(year != 2 ~ NA_real_,
                                 diab_ind==2 ~ 1,
                                 TRUE ~ 0),
         diab_new_y5 = case_when(year !=5 ~ NA_real_,
                                 diab_ind==2 ~ 1,
                                 TRUE ~ 0),
         
         diab_new_y7 = case_when(year != 7 ~ NA_real_,
                                 ((dmagediag - age) >= 0 & (dmagediag - age) <= 1) ~ 1,
                                 (glucosef >= 126 & !is.na(glucosef)) ~ 1,
                                  TRUE ~ 0),

         
         diab_new_y10 = case_when(year != 10 ~ NA_real_,
                                  ((dmagediag - age) >= 0 & (dmagediag - age) <= 1) ~ 1,
                                 (glucosef >= 126 & !is.na(glucosef)) ~ 1,
                                 TRUE ~ 0),
         
         diab_new_y15 = case_when(year != 15 ~ NA_real_,
                                  ((dmagediag - age) >= 0 & (dmagediag - age) <= 1) ~ 1,
                                 (glucosef >= 126 & !is.na(glucosef)) ~ 1,
                                 TRUE ~ 0),
         
         diab_new_y20 = case_when(year != 20 ~ NA_real_,
                                  ((dmagediag - age) >= 0 & (dmagediag - age) <= 1) ~ 1,
                                 (glucosef >= 126 & !is.na(glucosef)) ~ 1,
                                 (hba1c >=6.5 & !is.na(hba1c)) ~ 1,
                                 TRUE ~ 0),
         
         diab_new_y25 = case_when(year != 25 ~ NA_real_,
                                  ((dmagediag - age) >= 0 & (dmagediag - age) <= 1) ~ 1,
                                 (glucosef >= 126 & !is.na(glucosef)) ~ 1,
                                 (hba1c >=6.5 & !is.na(hba1c)) ~ 1,
                                 TRUE ~ 0),
         
         diab_new_y30 = case_when(year != 30 ~ NA_real_,
                                  ((dmagediag - age) >= 0 & (dmagediag - age) <= 1) ~ 1,
                                 (glucosef >= 126 & !is.na(glucosef)) ~ 1,
                                 (hba1c >=6.5 & !is.na(hba1c)) ~ 1,
                                 TRUE ~ 0)
                                 
         ) %>% 
  group_by(study_id) %>% 
  mutate(across(matches("(diab_y|diab_new_y)"),function(x) zoo::na.locf(x,na.rm=FALSE))) %>% 
  dplyr::filter(!is.na(age)) %>%
  mutate(earliest_age = min(age,na.rm=TRUE),
         dmagediag_locf = zoo::na.locf(dmagediag,na.rm=FALSE)) %>% 
  # mutate(across(matches("(diab_new_v)"),function(x) zoo::na.locf(x,fromLast = TRUE,na.rm=FALSE))) %>% 
  ungroup() %>% 
  distinct(study_id,age,bmi,.keep_all = TRUE) %>% 
  # CHECK: Do we use this? ----------
  mutate(dmdiagyear = case_when(# !is.na(dmagediag) & dmagediag <= earliest_age ~ 0,
                                 diab_y0 == 1 ~ 0,
                                 diab_new_y2 == 1 ~ 2,
                                 diab_new_y5 == 1 ~ 5,
                                 diab_new_y7 == 1 ~ 7,
                                 diab_new_y10 == 1 ~ 10,
                                 diab_new_y15 == 1 ~ 15,
                                 diab_new_y20 == 1 ~ 20,
                                 diab_new_y25 == 1 ~ 25,
                                 diab_new_y30 == 1 ~ 30,
                                 TRUE ~ NA_real_))


  # cardia_dat_all %>% 
  #   dplyr::filter(diab_y0 == 1 | diab_new_y2 == 1) %>% 
  #   distinct(study_id) %>% 
  #   nrow() --> 70: confirming with id_sel2


cardia_dat_all %>%
  # dplyr::filter(diab_y0 == 1 | diab_new_y2 == 1) %>%
  # dplyr::filter(diab_y0 == 1 | diab_new_y2 == 1 | diab_new_y5 == 1) %>%
  # dplyr::filter(diab_y0 == 1 | diab_new_y2 == 1 | diab_new_y5 == 1|diab_new_y7 == 1) %>%
  # dplyr::filter(diab_y0 == 1 | diab_new_y2 == 1 | diab_new_y5 == 1|diab_new_y7 == 1 | diab_new_y10 == 1) %>%
  # dplyr::filter(diab_y0 == 1 | diab_new_y2 == 1 | diab_new_y5 == 1|diab_new_y7 == 1 | diab_new_y10 == 1 | diab_new_y15 == 1) %>%
  dplyr::filter(diab_y0 == 1 | diab_new_y2 == 1 | diab_new_y5 == 1|diab_new_y7 == 1 | diab_new_y10 == 1 | diab_new_y15 == 1 |
                  diab_new_y20 == 1 | diab_new_y25 == 1 | diab_new_y30 == 1 ) %>%
  distinct(study_id) %>%
  nrow() 
# --> 70: confirming with id_sel2
#-->121: confirming with id_sel5
# --> 152: confirming with id_sel7
# --> 199: confirming with id_sel10
# --> 244: confirming with id_sel15
# --> 666: confirming with id_sel30


cardia_events = cardia_dat_all %>% 
  dplyr::filter(year == dmdiagyear) %>% 
  dplyr::select(study_id,dmagediag_locf,dmdiagyear,earliest_age,age,year) %>% 
  mutate(dmagediag = case_when(age <= dmagediag_locf ~ age,
                               !is.na(dmagediag_locf) ~ dmagediag_locf,
                               TRUE ~ age)) 

with(cardia_events,table((age-dmagediag) <= 1))


write_csv(cardia_events,paste0(path_diabetes_subphenotypes_adults_folder,"/working/qc/dsppre01b_cardia_events from diabetes_subphenotypes_predictors.csv"))


# cardia_longitudinal -----------

cardia_longitudinal = cardia_dat_all %>% 
  dplyr::select(-dmagediag) %>% 
  # Bringing the updated dmagediag from cardia_events
  left_join(cardia_events %>% 
              dplyr::select(-age,-earliest_age,-dmdiagyear,-dmagediag_locf) %>% 
              dplyr::select(-year),
            by=c("study_id")) %>% 
  mutate(
    available_labs = rowSums(!is.na(.[,lab_vars])),
    available_anthro = rowSums(!is.na(.[,anthro_vars])))



# Before dmagediag
before_dmagediag = join_by(study_id == original_study_id,
                           age < dmagediag_cluster)


cardia_newdm_not_in_longitudinal = cardia_newdm %>% 
  anti_join(cardia_longitudinal %>% 
              dplyr::select(study_id),
            by=c("original_study_id" = "study_id"))

cardia_longitudinal_newdm = cardia_longitudinal %>% 
  inner_join(cardia_newdm %>% 
               dplyr::select(original_study_id,dmagediag) %>% 
               rename(dmagediag_cluster = dmagediag),
             by = before_dmagediag) %>% 
  # !is.na(hba1c) -- not collected for most participants
  dplyr::filter(!is.na(hba1c)|!is.na(glucosef),!is.na(bmi)) %>% 
  mutate(diff_dmagediag = dmagediag - age)

cardia_longitudinal_neverdm = cardia_longitudinal %>% 
  dplyr::filter(is.na(dmagediag)) %>% 
  # !is.na(hba1c) -- not collected for 25% participants
  dplyr::filter(!is.na(hba1c)|!is.na(glucosef),!is.na(bmi)) %>% 
  # Among study waves where fasting glucose, HbA1c and BMI are measured, get the penultimate (second-to-last) wave
  group_by(study_id) %>% 
  mutate(wave = 1:n()) %>% 
  mutate(diff_next = case_when(wave == n() ~ NA_real_,
                               TRUE ~ dplyr::lead(age,1) - age)) %>% 
  ungroup() 

cardia_selected = bind_rows(
  cardia_longitudinal_newdm %>% 
    # dplyr::filter(diff_dmagediag >= 9/12,diff_dmagediag <= 15/12) %>% 
    mutate(type = "pre_1y_newdm") %>% 
    rename(diff_age = diff_dmagediag),
  cardia_longitudinal_neverdm %>% 
    # dplyr::filter(!is.na(diff_next),diff_next >= 9/12,diff_next <= 15/12) %>% 
    mutate(type = "pre_1y_neverdm") %>% 
    rename(diff_age = diff_next)
  
) %>% 
  dplyr::select(type,study_id,diff_age,age,available_labs,available_anthro,one_of(anthro_vars),one_of(lab_vars)) %>% 
  arrange(type,study_id,age)

saveRDS(cardia_selected,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01b_cardia.RDS"))
cardia_selected <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01b_cardia.RDS"))

