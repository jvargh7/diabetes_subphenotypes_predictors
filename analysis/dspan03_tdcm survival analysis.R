rm(list = ls());gc();source(".Rprofile")

library(mice)
library(survival)
library(survminer)
library(ggsurvfit)
library(broom)

source("functions/egfr_ckdepi_2021.R")
mi_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dsphyc301_mi_dfs.RDS"))

clean_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>% 
  dplyr::filter(dpp_intervention == 1) %>% 
  distinct(study,study_id,joint_id,dpp_intervention) # n = 1,772


# TDCM - longitudinal data, hazards time-varying, HR constant
analytic_dfs <- list()
overall_tdcm <- list()
mard_tdcm <- list()
mod_tdcm <- list()
sidd_tdcm <- list()
sird_tdcm <- list()


for(i in 1:mi_dfs$m) {
  df <- complete(mi_dfs, action = i) %>% 
    rename(joint_id = original_joint_id) %>% 
    mutate(egfr_ckdepi_2021 = egfr_ckdepi_2021(scr = serumcreatinine,female = female,age = age),
           time_to_event = censored_age - age) %>% 
    left_join(clean_df,
              by = c("study","study_id","joint_id")) %>% 
    mutate(dpp_intervention = case_when(
      dpp_intervention == 1 ~ 1,
      TRUE ~ 0
    ))
  
  
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
  
  
  tdcm_df <- df %>%
    arrange(study, study_id, joint_id, age) %>%
    group_by(study, study_id, joint_id) %>% 
    mutate(
      tstart = age, 
      tstop = dplyr::lead(age),
      event_true = dplyr::lead(event)
    ) %>% 
    ungroup() %>% 
    dplyr::filter(age < censored_age) %>% 
    mutate(mard = case_when(event_true == 1 & (cluster == "MARD") ~ 1,
                            TRUE ~ 0),
           mod = case_when(event_true == 1 & (cluster == "MOD") ~ 1,
                           TRUE ~ 0),
           sidd = case_when(event_true == 1 & (cluster == "SIDD") ~ 1,
                            TRUE ~ 0),
           sird = case_when(event_true == 1 & (cluster == "SIRD") ~ 1,
                            TRUE ~ 0)) %>% 
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
  

  overall_tdcm[[i]] <- coxph(Surv(tstart, tstop, event_true) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                             + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention, 
                             data = tdcm_df, cluster = joint_id)
  
  mard_tdcm[[i]] <- coxph(Surv(tstart, tstop, mard) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                          + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention, 
                          data = tdcm_df, cluster = joint_id)
  
  mod_tdcm[[i]] <- coxph(Surv(tstart, tstop, mod) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                         + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention, 
                         data = tdcm_df, cluster = joint_id)
  
  sidd_tdcm[[i]] <- coxph(Surv(tstart, tstop, sidd) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                          + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention, 
                          data = tdcm_df, cluster = joint_id)
  
  # CARDIA has 0 people from SIRD 
  tdcm_sird <- tdcm_df %>% dplyr::filter(!study %in% c("cardia"))
  sird_tdcm[[i]] <- coxph(Surv(tstart, tstop, sird) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                          + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled, 
                          data = tdcm_df, cluster = joint_id)
  
  
  
  


}


unpooled_resuts <- bind_rows(
  imap_dfr(overall_tdcm,function(m,index){
    broom::tidy(m) %>% mutate(model = "Overall",iteration = index) %>% 
      return(.)
  }),
  imap_dfr(mard_tdcm,function(m,index){
    broom::tidy(m) %>% mutate(model = "MARD",iteration = index) %>% 
      return(.)
  }),
  imap_dfr(mod_tdcm,function(m,index){
    broom::tidy(m) %>% mutate(model = "MOD",iteration = index) %>% 
      return(.)
  }),
  imap_dfr(sidd_tdcm,function(m,index){
    broom::tidy(m) %>% mutate(model = "SIDD",iteration = index) %>% 
      return(.)
  }),
  imap_dfr(sird_tdcm,function(m,index){
    broom::tidy(m) %>% mutate(model = "SIRD",iteration = index) %>% 
      return(.)
  })
  
  
)

write_csv(unpooled_resuts, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan03_unpooled tdcm results.csv"))

source("functions/pool_results.R")

tdcm_results <- bind_rows(
  pool_results(overall_tdcm) %>% mutate(model = "Overall"),
  pool_results(mard_tdcm) %>% mutate(model = "MARD"),
  pool_results(mod_tdcm) %>% mutate(model = "MOD"),
  pool_results(sidd_tdcm) %>% mutate(model = "SIDD"),
  pool_results(sird_tdcm) %>% mutate(model = "SIRD")
  ) %>% 
  write_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan03_pooled tdcm results.csv"))

