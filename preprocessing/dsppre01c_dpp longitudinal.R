
rm(list=ls());gc();source(".Rprofile")

anthro_vars <- c("sbp","dbp","height","weight","wc","hc","triceps","iliac","abdominal","medial","bmi")
# "totalc" -- not there in dpp
lab_vars <- c("hba1c","insulinf","glucosef","glucose2h","vldlc","tgl","hdlc","ldlc","totalc",
              "serumcreatinine","ast","alt")
med_vars <- c("med_dm_use")

# Each 'release' contributed a row
dpp_demographics <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/dpppre01_demographic.RDS")) %>% 
  distinct(study_id,.keep_all=TRUE)

# Each 'release' contributed a row
dos_demographics <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/dospre01_demographic.RDS")) %>% 
  distinct(study_id,.keep_all=TRUE)

# Diabetes outcomes ------------
# Among intensive lifestyle, metformin, and placebo # participants, 849 had been diagnosed as having diabetes as of September, 2002, 
# and another 503 participants developed diabetes during the first phase of DPPOS
dpp_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/dpp_newdm.RDS")) 

dos_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/dos_newdm.RDS")) %>% 
  anti_join(dpp_newdm,
            by=c("study_id")) 

table(dos_newdm$study_id %in% dpp_newdm$study_id)
table(dpp_newdm$study_id %in% dos_newdm$study_id)

dos_nodm <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/dospre03_events.RDS")) %>% 
  anti_join(dpp_newdm,
            by=c("study_id")) %>% 
  anti_join(dos_newdm,
            by=c("study_id")) %>% 
  mutate(diagDays = round(diabt*365.25)) %>% 
  distinct(study_id,diagDays,diabf) %>% 
  group_by(study_id,diagDays) %>% 
  summarize(diabf = max(diabf),
            n = n()) %>% 
  ungroup() %>% 
  # There were some cases when diabf was both 0 and 1 --> checked a few against actual labs of glucosef and glucose2h 
  dplyr::select(-n) %>% 
  dplyr::filter(diabf == 0) %>% 
  group_by(study_id) %>% 
  dplyr::filter(diagDays == max(diagDays)) %>% 
  ungroup()

dpp_nodm <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/dpppre03_events.RDS")) %>% 
  anti_join(dpp_newdm,
            by=c("study_id")) %>% 
  anti_join(dos_newdm,
            by=c("study_id")) %>% 
  anti_join(dos_nodm,
            by=c("study_id")) %>% 
  mutate(diagDays = round(diabt*365.25)) %>% 
  distinct(study_id,diagDays,diabf) %>% 
  group_by(study_id,diagDays) %>% 
  summarize(diabf = max(diabf),
            n = n()) %>% 
  ungroup() %>% 
  # There were some cases when diabf was both 0 and 1 --> checked a few against actual labs of glucosef and glucose2h 
  dplyr::select(-n) %>% 
  dplyr::filter(diabf == 0)   %>% 
  group_by(study_id) %>% 
  dplyr::filter(diagDays == max(diagDays)) %>% 
  ungroup()




dppos_max_diagDays = bind_rows(
  dpp_newdm %>% mutate(dpp = 1, newdm = 1),
  dos_newdm %>% mutate(dpp = 0, newdm = 1),
  dpp_nodm %>% mutate(dpp = 1, newdm = 0),
  dos_nodm %>% mutate(dpp = 0, newdm = 0),
) %>% 
  dplyr::select(study_id,dpp,newdm,diagDays)

# saveRDS(dppos_max_diagDays,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01c_dpp_max_diagDays.RDS"))




# There were some cases when the same 'lab_StudyDays' had multiple  values for the same parameter
dpp_lab <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/dpppre02_labs.RDS")) %>% 
  rename(lab_StudyDays = StudyDays) %>% 
  # -release -- not there in dpp; only in dppos
  dplyr::select(-quarter,-semi,-quarter,-visit) %>% 
  group_by(study_id,lab_StudyDays) %>% 
  summarize(across(everything(),~mean(.,na.rm=TRUE))) %>% 
  ungroup() 

dos_lab <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/dospre02_labs.RDS")) %>% 
  rename(lab_StudyDays = StudyDays) %>% 
  # -release -- not there in dpp; only in dppos
  dplyr::select(-quarter,-semi,-quarter,-release,-visit) %>% 
  group_by(study_id,lab_StudyDays) %>% 
  summarize(across(everything(),~mean(.,na.rm=TRUE))) %>% 
  ungroup()



dpp_anthro <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/dpppre04_anthro.RDS")) %>% 
  rename(anthro_StudyDays = StudyDays) %>% 
  dplyr::select(-visit,-semi,-quarter,-visit) %>% 
  group_by(study_id,anthro_StudyDays) %>% 
  summarize(across(everything(),~mean(.,na.rm=TRUE))) %>% 
  ungroup()


dos_anthro <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/dospre04_anthro.RDS")) %>% 
  rename(anthro_StudyDays = StudyDays) %>% 
  dplyr::select(-visit,-semi,-quarter,-visit) %>% 
  group_by(study_id,anthro_StudyDays) %>% 
  summarize(across(everything(),~mean(.,na.rm=TRUE))) %>% 
  ungroup()


lab = bind_rows(dpp_lab,
                dos_lab) %>% 
  left_join(dppos_max_diagDays,
            by=c("study_id")) %>% 
  # Removed the below to get later visits also
  # dplyr::filter(lab_StudyDays <= diagDays) %>% 
  
  # Keep all observations for 'nodm' and only pre-diagnosis matches for 'newdm'
  # This is because for newdm == 1,  diagDays - lab_StudyDays are %in% c(0,365)
  dplyr::filter(newdm == 0 | (newdm == 1 & lab_StudyDays < diagDays) | (newdm == 1 & lab_StudyDays > (diagDays + 365)))

anthro = bind_rows(dpp_anthro,
                dos_anthro) %>% 
  left_join(dppos_max_diagDays,
            by=c("study_id")) %>% 
  
  # Removed the below to get later visits also
  # dplyr::filter(anthro_StudyDays <= diagDays) %>%  
  
  # Keep all observations for 'nodm' and only pre-diagnosis matches for 'newdm'
  # This is because for newdm == 1,  diagDays - anthro_StudyDays are %in% c(0,365)
  
  dplyr::filter(newdm == 0 | (newdm == 1 & anthro_StudyDays < diagDays) | (newdm == 1 & anthro_StudyDays > (diagDays + 365)))


anthro_lab = join_by(study_id==study_id,
                     closest(lab_StudyDays >= anthro_StudyDays),
                    # Decided an arbitrary cutoff of 90 days given the proximity in values and follow-up every 6 months
                     closest(lab_StudyDays <= anthro_StudyDays_plus90))


dppos_longitudinal = lab %>%
  left_join(anthro %>% 
              dplyr::select(-dpp,-newdm,-diagDays) %>% 
              mutate(anthro_StudyDays_plus90 = anthro_StudyDays + 90),
            # Decided to use a rolling join defined as 'anthro_lab'
            # by=c("study_id","lab_StudyDays"="anthro_StudyDays"))   %>% 
            by=anthro_lab)   %>%
  dplyr::select(-anthro_StudyDays_plus90) %>% 
  inner_join(dpp_demographics,
             by=c("study_id")) %>% 
  left_join(
    bind_rows(dpp_newdm ,
              dos_newdm) %>% 
      dplyr::select(study_id,dmagediag),
    by=c("study_id")
    
  ) %>% 
  mutate(age = case_when(agegroup == 1 ~ 37 + lab_StudyDays/365, # Less than 40
                         agegroup == 2 ~ 42 + lab_StudyDays/365,
                         agegroup == 3 ~ 47 + lab_StudyDays/365,
                         agegroup == 4 ~ 52 + lab_StudyDays/365,
                         agegroup == 5 ~ 57 + lab_StudyDays/365,
                         agegroup == 6 ~ 62 + lab_StudyDays/365,
                         agegroup == 7 ~ 67 + lab_StudyDays/365)) %>% 
  dplyr::filter(!is.na(bmi)) %>% 
  mutate(female = case_when(sex == 1 ~ 0,
                            sex == 2 ~ 1,
                            TRUE ~ NA_real_)) %>% 
  mutate(med_dm_use = case_when(
    treatment %in% c("Troglitazone","Metformin") ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(age = round(age,2),
         available_labs = rowSums(!is.na(.[,lab_vars])),
         available_anthro = rowSums(!is.na(.[,anthro_vars]))) %>% 
  mutate(race_clean = race_eth) %>% 
  dplyr::select(study_id,dpp,newdm,age,dmagediag,diagDays,race_eth,race_clean,female,
                lab_StudyDays,anthro_StudyDays,available_labs,available_anthro,
                one_of(anthro_vars),one_of(lab_vars),one_of(med_vars)) %>% 
  group_by(study_id) %>%  
  arrange(study_id,lab_StudyDays,age) %>% 
  mutate(wave = row_number()) %>% 
  ungroup()  



saveRDS(dppos_longitudinal,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01c_dppos.RDS"))

dppos_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01c_dppos.RDS"))
write_csv(dppos_longitudinal,paste0(path_prediabetes_subphenotypes_folder,"/working/longitudinal/dpp and dppos.csv"))
