rm(list = ls());gc();source(".Rprofile")

library(mice)
library(survival)
library(survminer)
library(ggsurvfit)
library(broom)

source("functions/egfr_ckdepi_2021.R")
mi_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/mi_dfs.RDS"))

# TDCM - longitudinal data, hazards time-varying, HR constant
analytic_dfs <- list()
overall_tdcm <- list()
mard_tdcm <- list()
mod_tdcm <- list()
sidd_tdcm <- list()
sird_tdcm <- list()


for(i in 1:mi_dfs$m) {
  df <- complete(mi_dfs, action = i) %>% 
    mutate(egfr_ckdepi_2021 = egfr_ckdepi_2021(scr = serumcreatinine,female = female,age = age),
           time_to_event = censored_age - age) 
  
  
  analytic_df <- df %>% 
    arrange(study,study_id,joint_id,age) %>% 
    group_by(study,study_id,joint_id) %>%
    mutate(event = case_when(
      newdm_event == 1 & (age == censored_age) ~ 1,  # event is 1 for the last wave
      TRUE ~ 0 
    )
    ) %>%
    ungroup() 
  
  analytic_dfs[[i]] <- analytic_df
}


for (i in 1:length(analytic_dfs)) {
  df <- analytic_dfs[[i]]   
  
  cluster_df <- analytic_df %>% 
    mutate(mard = case_when(event == 1 & (cluster == "MARD") ~ 1,
                            TRUE ~ 0),
           mod = case_when(event == 1 & (cluster == "MOD") ~ 1,
                           TRUE ~ 0),
           sidd = case_when(event == 1 & (cluster == "SIDD") ~ 1,
                            TRUE ~ 0),
           sird = case_when(event == 1 & (cluster == "SIRD") ~ 1,
                            TRUE ~ 0))
  
  tdcm_df <- cluster_df %>%
    arrange(study, study_id, joint_id, age) %>%
    group_by(study, study_id, joint_id) %>% 
    mutate(
      tstart = case_when(row_number() == 1 ~ age, 
                         TRUE ~ dplyr::lag(age, n = 1)), 
      tstop = age
    ) %>% 
    ungroup() %>% 
    dplyr::filter((tstart < tstop) & (tstop <= censored_age)) %>% 
    mutate(across(c(bmi, hba1c, homa2b, homa2ir, ldlc, sbp, egfr_ckdepi_2021), ~replace(., is.infinite(.), NA))) %>% 
    # error due to 0 ppl in NH Other (sidd == 1), ignore this category
    mutate(race = case_when(race == "NH Other" ~ "Other", 
                            TRUE ~ race)) %>% 
    # scaling
    mutate(sbp_scaled = sbp/10,
           ldlc_scaled = ldlc/10,
           homa2b_scaled = homa2b/10,
           egfr_ckdepi_2021_scaled = egfr_ckdepi_2021/10)
  

  overall_tdcm[[i]] <- coxph(Surv(tstart, tstop, event) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                             + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled, 
                             data = tdcm_df)
  
  mard_tdcm[[i]] <- coxph(Surv(tstart, tstop, mard) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                          + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled, 
                          data = tdcm_df)
  
  mod_tdcm[[i]] <- coxph(Surv(tstart, tstop, mod) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                         + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled, 
                         data = tdcm_df)
  
  sidd_tdcm[[i]] <- coxph(Surv(tstart, tstop, sidd) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                          + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled, 
                          data = tdcm_df)
  
  # CARDIA has 0 people from SIRD 
  tdcm_sird <- tdcm_df %>% dplyr::filter(!(study == "cardia"))
  sird_tdcm[[i]] <- coxph(Surv(tstart, tstop, sird) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                          + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled, 
                          data = tdcm_sird)


}



source("functions/pool_results.R")

tdcm_results <- bind_rows(
  pool_results(overall_tdcm) %>% mutate(model = "Overall"),
  pool_results(mard_tdcm) %>% mutate(model = "MARD"),
  pool_results(mod_tdcm) %>% mutate(model = "MOD"),
  pool_results(sidd_tdcm) %>% mutate(model = "SIDD"),
  pool_results(sird_tdcm) %>% mutate(model = "SIRD")
  ) %>% 
  write_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan03_pooled tdcm results.csv"))

