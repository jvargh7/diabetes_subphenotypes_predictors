rm(list = ls());gc();source(".Rprofile")

source("functions/egfr_ckdepi_2021.R")

# 1. Use "joint_id" to check N
# 2. New T2D total: N = 7,623
# 3. T2D last visit defined as within the 1y after T2D diagnosis


# new dm with non-missing age at diagnosis, N = 7,623
final_dataset_temp = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/final_dataset_temp.RDS")) %>%
  mutate(joint_id = paste(study, original_study_id, sep = "_"))

# new T2D (n = 7,623) + No T2D (n = 22,373)
analytic_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dsppre01_analytic df.RDS")) %>% 
  mutate(egfr_ckdepi_2021 = egfr_ckdepi_2021(scr = serumcreatinine,female = female,age = age),
         joint_id = paste(study, study_id, sep = "_"),
         # newly diagnosed T2D or no T2D
         newdm_event = case_when(joint_id %in% final_dataset_temp$joint_id ~ 1,
                                 is.na(dmagediag) ~ 0,
                                 TRUE ~ NA)) %>%
  dplyr::filter(!is.na(newdm_event))

# Inclusion - Exclusion (N) --------------------------------------------

# new T2D, N = 7,623
newt2d_df <- analytic_df %>% 
  dplyr::filter(newdm_event == 1) 
  
# ID with prior A1c, N = 2,216
prior_a1c_newt2d <- newt2d_df %>% 
  dplyr::filter(age < dmagediag) %>%
  dplyr::filter(!is.na(hba1c)) %>% 
  arrange(joint_id, age) %>%                      
  group_by(joint_id) %>%                          
  slice_head(n = 1) %>%                          
  ungroup() %>% 
  mutate(earliest_age_newt2d = age)

# ID without prior A1c, N = 5,407
noprior_a1c_newt2d <- newt2d_df %>% 
  dplyr::filter(!joint_id %in% prior_a1c_newt2d$joint_id) %>% 
  distinct(joint_id)

# Prior A1c ---------------------------

## T2D diagnosis ## N = 2,216
t2d_diag_a1c <- newt2d_df %>% 
  dplyr::filter(joint_id %in% prior_a1c_newt2d$joint_id) %>% 
  # To address decimal values
  mutate(diff = age - dmagediag) %>%
  dplyr::filter(diff >= 0, diff <= 1) %>%
  group_by(joint_id) %>%
  slice_min(order_by = diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(max_age_newt2d = age) 

# 1. A1c + BMI, N = 1,570
a1c_bmi_newt2d <- t2d_diag_a1c %>% 
  dplyr::filter(!is.na(hba1c) & !is.na(bmi)) %>% 
  distinct(joint_id, max_age_newt2d)

# 2. A1c + no BMI, N = 22
a1c_nobmi <- t2d_diag_a1c %>% 
  dplyr::filter(!is.na(hba1c) & is.na(bmi)) %>% 
  distinct(joint_id)

# 3. no A1c + BMI, N = 579
noa1c_bmi <- t2d_diag_a1c %>% 
  dplyr::filter(is.na(hba1c) & !is.na(bmi)) %>% 
  distinct(joint_id)

# 4. no A1c + no BMI, N = 45
noa1c_nobmi <- t2d_diag_a1c %>% 
  dplyr::filter(is.na(hba1c) & is.na(bmi)) %>% 
  distinct(joint_id)


# No Prior A1c ---------------------------

## T2D diagnosis ## N = 5,407
t2d_diag <- newt2d_df %>% 
  dplyr::filter(joint_id %in% noprior_a1c_newt2d$joint_id) %>% 
  # To address decimal values
  mutate(diff = age - dmagediag) %>%
  dplyr::filter(diff >= 0, diff <= 1) %>%
  group_by(joint_id) %>%
  slice_min(order_by = diff, n = 1, with_ties = FALSE) %>%
  ungroup()

# 1. A1c + BMI, N = 1,819
a1c_bmi <- t2d_diag %>% 
  dplyr::filter(!is.na(hba1c) & !is.na(bmi)) %>% 
  distinct(joint_id)

# 2. A1c + no BMI, N = 124
a1c_nobmi <- t2d_diag %>% 
  dplyr::filter(!is.na(hba1c) & is.na(bmi)) %>% 
  distinct(joint_id)

# 3. no A1c + BMI, N = 3,404
noa1c_bmi <- t2d_diag %>% 
  dplyr::filter(is.na(hba1c) & !is.na(bmi)) %>% 
  distinct(joint_id)

# 4. no A1c + no BMI, N = 60
noa1c_nobmi <- t2d_diag %>% 
  dplyr::filter(is.na(hba1c) & is.na(bmi)) %>% 
  distinct(joint_id)





# No T2D, N = 22,373 -----------------------------------
not2d_df <- analytic_df %>% 
  dplyr::filter(newdm_event == 0) 

# last follow-up with A1c, N = 14,598
last_a1c <- not2d_df %>% 
  dplyr::filter(!is.na(hba1c)) %>% 
  arrange(joint_id, age) %>%                      
  group_by(joint_id) %>%                          
  slice_tail(n = 1) %>%                          
  ungroup() %>% 
  mutate(max_age_not2d = age) %>% 
  distinct(joint_id,max_age_not2d,hba1c,bmi)

# ID last follow-up with A1c, at least 1 prior visit N = 14,232
last_a1c_select <- not2d_df %>% 
  dplyr::filter(joint_id %in% last_a1c$joint_id) %>% 
  left_join(last_a1c, by = "joint_id") %>% 
  group_by(joint_id) %>% 
  mutate(has_age_before_max = any(age < max_age_not2d)) %>% 
  ungroup() %>% 
  dplyr::filter(has_age_before_max) %>% 
  select(-has_age_before_max)

# A1c available at only 1 visit, N = 366 -- assign to No prior A1c
a1c_1visit <- last_a1c %>% 
  dplyr::filter(!joint_id %in% last_a1c_select$joint_id) %>% 
  dplyr::filter(!is.na(hba1c) & !is.na(bmi)) # N = 366
  # dplyr::filter(!is.na(hba1c) & is.na(bmi)) # N = 0

# ID with prior A1c, N = 8,775
prior_a1c_not2d <- last_a1c_select %>% 
  dplyr::filter(age < max_age_not2d) %>%
  dplyr::filter(!is.na(hba1c)) %>% 
  arrange(joint_id, age) %>%                      
  group_by(joint_id) %>%                          
  slice_head(n = 1) %>%                          
  ungroup() %>% 
  mutate(earliest_age_not2d = age)

# ID without prior A1c, N = 5,457
noprior_a1c <- last_a1c_select %>% 
  dplyr::filter(!joint_id %in% prior_a1c_not2d$joint_id) %>% 
  distinct(joint_id)

# Prior A1c ---------------------------

## last follow-up ## N = 8,775
last_fup <- last_a1c_select %>% 
  dplyr::filter(joint_id %in% prior_a1c_not2d$joint_id) %>% 
  dplyr::filter(age == max_age_not2d)

# 1. A1c + BMI, N = 8,747
a1c_bmi_not2d <- last_fup %>% 
  dplyr::filter(!is.na(hba1c) & !is.na(bmi)) %>% 
  distinct(joint_id,max_age_not2d)

# 2. A1c + no BMI, N = 28
a1c_nobmi <- last_fup %>% 
  dplyr::filter(!is.na(hba1c) & is.na(bmi)) %>% 
  distinct(joint_id)

# 3. no A1c + BMI, N = 0
noa1c_bmi <- last_fup %>% 
  dplyr::filter(is.na(hba1c) & !is.na(bmi)) %>% 
  distinct(joint_id)

# 4. no A1c + no BMI, N = 0
noa1c_nobmi <- last_fup %>% 
  dplyr::filter(is.na(hba1c) & is.na(bmi)) %>% 
  distinct(joint_id)


# No Prior A1c ---------------------------

## T2D diagnosis ## N = 5,457
last_fup <- last_a1c_select %>% 
  dplyr::filter(joint_id %in% noprior_a1c$joint_id) %>% 
  dplyr::filter(age == max_age_not2d)

# 1. A1c + BMI, N = 5,420
a1c_bmi <- last_fup %>% 
  dplyr::filter(!is.na(hba1c) & !is.na(bmi)) %>% 
  distinct(joint_id)

# 2. A1c + no BMI, N = 37
a1c_nobmi <- last_fup %>% 
  dplyr::filter(!is.na(hba1c) & is.na(bmi)) %>% 
  distinct(joint_id)

# 3. no A1c + BMI, N = 0
noa1c_bmi <- last_fup %>% 
  dplyr::filter(is.na(hba1c) & !is.na(bmi)) %>% 
  distinct(joint_id)

# 4. no A1c + no BMI, N = 0
noa1c_nobmi <- last_fup %>% 
  dplyr::filter(is.na(hba1c) & is.na(bmi)) %>% 
  distinct(joint_id)





##### final analytic sample ######

# T2D: 1,570, OBS = 16,556; NoT2D: 8,747, OBS = 57,089
selected_df <- analytic_df %>% 
  dplyr::filter(joint_id %in% c(a1c_bmi_newt2d$joint_id, a1c_bmi_not2d$joint_id)) %>% 
  left_join(prior_a1c_newt2d %>% 
              select(joint_id, earliest_age_newt2d),
            by = "joint_id") %>% 
  left_join(prior_a1c_not2d %>% 
              select(joint_id, earliest_age_not2d),
            by = "joint_id") %>% 
  left_join(a1c_bmi_newt2d %>% 
              select(joint_id, max_age_newt2d),
            by = c("joint_id")) %>%
  left_join(a1c_bmi_not2d %>% 
              select(joint_id, max_age_not2d),
            by = c("joint_id")) %>%
  mutate(earliest_age = case_when(newdm_event == 1 ~ earliest_age_newt2d,
                                  TRUE ~ earliest_age_not2d),
         censored_age = case_when(newdm_event == 1 ~ max_age_newt2d, # defined as the age within 0-1y after T2D diagnosis
                                  TRUE ~ max_age_not2d)) %>% 
  dplyr::filter(age >= earliest_age, age <= (censored_age)) %>% 
  select(-earliest_age_newt2d,-earliest_age_not2d,-max_age_not2d,-max_age_newt2d) %>% 
  group_by(study,study_id,joint_id) %>% 
  mutate(t = age - censored_age) %>% 
  ungroup()

# redefine outliers
clean_df <- selected_df %>% 
  mutate(
    weight = case_when(weight>150 ~ 150,
                       TRUE ~ weight),
    bmi = case_when(bmi<15 ~ 15,
                    bmi>65 ~ 65,
                    TRUE ~ bmi),
    wc = case_when(wc<50 ~ 50,
                   wc>150 ~ 150,
                   TRUE ~ wc),
    sbp = case_when(sbp<70 ~ 70,
                    sbp>220 ~ 220,
                    TRUE ~ sbp),
    dbp = case_when(dbp<40 ~ 40,
                    dbp>130 ~ 130,
                    TRUE ~ dbp),
    hba1c = case_when(hba1c<4 ~ 4,
                      hba1c>14 ~ 14,
                      TRUE ~ hba1c),
    ldlc = case_when(ldlc<30 ~ 30,
                     ldlc>250 ~ 250,
                     TRUE ~ ldlc),
    hdlc = case_when(hdlc<10 ~ 10,
                     hdlc>100 ~ 100,
                     TRUE ~ hdlc),
    vldlc = case_when(vldlc<5 ~ 5,
                      vldlc>100 ~ 100,
                      TRUE ~ vldlc),
    glucosef = case_when(glucosef<50 ~ 50,
                         glucosef>150 ~ 150,
                         TRUE ~ glucosef),
    insulinf = case_when(insulinf<2 ~ 2,
                         insulinf>100 ~ 100,
                         TRUE ~ insulinf),
    glucose2h = case_when(glucose2h<40 ~ 40,
                          glucose2h>300 ~ 300,
                          TRUE ~ glucose2h),
    tgl = case_when(tgl<30 ~ 30,
                    tgl>300 ~ 300,
                    TRUE ~ tgl),
    egfr = case_when(egfr<5 ~ 5,
                     egfr>150 ~ 150,
                     TRUE ~ egfr),
    homa2b = case_when(homa2b<10 ~ 10,
                       homa2b>300 ~ 300,
                       TRUE ~ homa2b),
    homa2ir = case_when(homa2ir>10 ~ 10,
                        TRUE ~ homa2ir)
  )


saveRDS(clean_df, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS"))

