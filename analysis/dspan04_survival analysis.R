rm(list = ls());gc();source(".Rprofile")

library(survival)
library(survminer)
library(ggsurvfit)
library(broom)

analytic_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan02_predictors analytic dfs.RDS"))

# TDCM - longitudinal data, hazards time-varying, HR constant

overall_tdcm <- list()
mard_tdcm <- list()
mod_tdcm <- list()
sidd_tdcm <- list()
sird_tdcm <- list()
not2d_tdcm <- list()

for (i in 1:length(analytic_dfs)) {
  df <- analytic_dfs[[i]]  
  
  cluster_df <- df %>% 
    mutate(mard = case_when(subtype == "MARD" ~ 1,
                            TRUE ~ 0),
           mod = case_when(subtype == "MOD" ~ 1,
                           TRUE ~ 0),
           sidd = case_when(subtype == "SIDD" ~ 1,
                            TRUE ~ 0),
           sird = case_when(subtype == "SIRD" ~ 1,
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
    mutate(across(c(bmi, hba1c, homa2b, homa2ir, ldlc, sbp, egfr_ckdepi_2021), ~replace(., is.infinite(.), NA))) %>% 
    dplyr::filter((tstart < tstop) & (tstop <= censored_age)) %>% 
    # error due to 0 ppl in NH Other (sidd == 1), ignore this category
    mutate(race = case_when(race == "NH Other" ~ "Other", 
                            TRUE ~ race)) %>% 
    # scaling
    mutate(sbp_scaled = sbp/10,
           ldlc_scaled = ldlc/10,
           egfr_ckdepi_2021_scaled = egfr_ckdepi_2021/10)
  

  overall_tdcm[[i]] <- coxph(Surv(tstart, tstop, event) ~ study + female + race + min_age + bmi + hba1c + homa2b 
                             + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled, 
                             data = tdcm_df)
  
  mard_tdcm[[i]] <- coxph(Surv(tstart, tstop, mard) ~ study + female + race + min_age + bmi + hba1c + homa2b 
                          + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled, 
                          data = tdcm_df)
  
  mod_tdcm[[i]] <- coxph(Surv(tstart, tstop, mod) ~ study + female + race + min_age + bmi + hba1c + homa2b 
                         + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled, 
                         data = tdcm_df)
  
  sidd_tdcm[[i]] <- coxph(Surv(tstart, tstop, sidd) ~ study + female + race + min_age + bmi + hba1c + homa2b 
                          + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled, 
                          data = tdcm_df)
  
  sird_tdcm[[i]] <- coxph(Surv(tstart, tstop, sird) ~ study + female + race + min_age + bmi + hba1c + homa2b 
                          + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled, 
                          data = tdcm_df)

}



source("functions/pool_results.R")

tdcm_results <- bind_rows(
  pool_results(overall_tdcm) %>% mutate(model = "Overall"),
  pool_results(mard_tdcm) %>% mutate(model = "MARD"),
  pool_results(mod_tdcm) %>% mutate(model = "MOD"),
  pool_results(sidd_tdcm) %>% mutate(model = "SIDD"),
  pool_results(sird_tdcm) %>% mutate(model = "SIRD")
  ) %>% 
  write_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan03_tdcm pooled results with multiple imputation.csv"))


