
rm(list=ls());gc();source(".Rprofile")

anthro_vars <- c("sbp","dbp","height","wc","hc","triceps","iliac","abdominal","medial","bmi")
# "totalc" -- not there in dpp
lab_vars <- c("hba1c","insulinf","glucosef","glucose2h","vldlc","tgl","hdlc","ldlc",
              "serumcreatinine","ast","alt")

dpp_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/processed/final_dataset_temp.RDS")) %>% 
  dplyr::filter(study == "dpp") %>% 
  mutate(original_study_id = as.numeric(original_study_id))  %>% 
  mutate(dmagediag = round(dmagediag,2))

# Each 'release' contributed a row
dpp_demographics <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/dpppre01_demographic.RDS")) %>% 
  distinct(study_id,.keep_all=TRUE)

# There were some cases when the same 'lab_StudyDays' had multiple  values for the same parameter
dpp_lab <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/dpppre02_labs.RDS")) %>% 
  rename(lab_StudyDays = StudyDays) %>% 
  # -release -- not there in dpp; only in dppos
  dplyr::select(-quarter,-semi,-quarter,-visit) %>% 
  group_by(study_id,lab_StudyDays) %>% 
  summarize(across(everything(),~mean(.,na.rm=TRUE))) %>% 
  ungroup()

dpp_events <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/dpppre03_events.RDS")) %>% 
  mutate(diagDays = round(diabt*365.25)) %>% 
  distinct(study_id,diagDays,diabf) %>% 
  group_by(study_id,diagDays) %>% 
  # There were some cases when diabf was both 0 and 1 --> checked a few against actual labs of glucosef and glucose2h 
  dplyr::filter(diabf == max(diabf)) %>% 
  ungroup()

dpp_anthro <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/dpppre04_anthro.RDS")) %>% 
  rename(anthro_StudyDays = StudyDays) %>% 
  dplyr::select(-visit,-semi,-quarter,-visit) %>% 
  group_by(study_id,anthro_StudyDays) %>% 
  summarize(across(everything(),~mean(.,na.rm=TRUE))) %>% 
  ungroup()

anthro_lab = join_by(study_id==study_id,
                     closest(lab_StudyDays >= anthro_StudyDays),
                     closest(lab_StudyDays <= anthro_StudyDays_plus90))

dpp_longitudinal = dpp_demographics  %>% 
  left_join(dpp_lab,
            by=c("study_id")) %>% 
  left_join(dpp_events, 
            by=c("study_id","lab_StudyDays"="diagDays")) %>%
  left_join(dpp_anthro %>% 
              mutate(anthro_StudyDays_plus90 = anthro_StudyDays + 90),
            # Decided to use a rolling join defined as 'anthro_lab'
            # by=c("study_id","lab_StudyDays"="anthro_StudyDays"))   %>% 
            by=anthro_lab)   %>%
  dplyr::select(-anthro_StudyDays_plus90) %>% 
  mutate(age = case_when(agegroup == 1 ~ 37 + lab_StudyDays/365, # Less than 40
                         agegroup == 2 ~ 42 + lab_StudyDays/365,
                         agegroup == 3 ~ 47 + lab_StudyDays/365,
                         agegroup == 4 ~ 52 + lab_StudyDays/365,
                         agegroup == 5 ~ 57 + lab_StudyDays/365,
                         agegroup == 6 ~ 62 + lab_StudyDays/365,
                         agegroup == 7 ~ 67 + lab_StudyDays/365)) %>% 
  mutate(age = round(age,2),
         available_labs = rowSums(!is.na(.[,lab_vars])),
         available_anthro = rowSums(!is.na(.[,anthro_vars])))


# Before dmagediag
before_dmagediag = join_by(study_id == original_study_id,
                           age < dmagediag)
dpp_longitudinal_newdm = dpp_longitudinal %>% 
  inner_join(dpp_newdm %>% 
               dplyr::select(original_study_id,dmagediag),
             by = before_dmagediag) %>% 
  dplyr::filter(!is.na(hba1c)|!is.na(glucosef),!is.na(bmi)) %>% 
  mutate(diff_dmagediag = dmagediag - age)

dpp_longitudinal_neverdm = dpp_longitudinal %>% 
  group_by(study_id) %>% 
  mutate(ever_diabf = sum(diabf,na.rm=TRUE)) %>% 
  ungroup() %>% 
  dplyr::filter(ever_diabf == 0) %>% 
  dplyr::filter(!is.na(hba1c)|!is.na(glucosef),!is.na(bmi)) %>% 
  # Among study waves where fasting glucose, HbA1c and BMI are measured, get the penultimate (second-to-last) wave
  group_by(study_id) %>% 
  mutate(wave = 1:n()) %>% 
  mutate(diff_next = case_when(wave == n() ~ NA_real_,
                               TRUE ~ dplyr::lead(age,1) - age)) %>% 
  ungroup() 


dpp_selected = bind_rows(
  dpp_longitudinal_newdm %>% 
    # dplyr::filter(diff_dmagediag >= 9/12,diff_dmagediag <= 15/12) %>% 
    mutate(type = "pre_1y_newdm") %>% 
    rename(diff_age = diff_dmagediag),
  dpp_longitudinal_neverdm %>% 
    # dplyr::filter(!is.na(diff_next),diff_next >= 9/12,diff_next <= 15/12) %>% 
    mutate(type = "pre_1y_neverdm") %>% 
    rename(diff_age = diff_next)
  
) %>% 
  dplyr::select(type,study_id,diff_age,age,available_labs,available_anthro,one_of(anthro_vars),one_of(lab_vars)) %>% 
  arrange(type,study_id,age)

saveRDS(dpp_selected,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01c_dpp.RDS"))
dpp_selected <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01c_dpp.RDS"))
