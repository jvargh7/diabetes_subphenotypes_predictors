rm(list = ls());gc();source(".Rprofile")

library(survival)
library(survminer)
library(ggsurvfit)
library(broom)

ipcw_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan02_ipcw dfs.RDS"))

# TDCM - longitudinal data, hazards time-varying, HR constant

overall_tdcm <- list()
mard_tdcm <- list()
mod_tdcm <- list()
sidd_tdcm <- list()
sird_tdcm <- list()


for (i in 1:length(ipcw_dfs)) {
  df <- ipcw_dfs[[i]]  
  
  cluster_df <- df %>% 
    mutate(mard = case_when(cluster == "MARD" ~ 1,
                            TRUE ~ 0),
           mod = case_when(cluster == "MOD" ~ 1,
                           TRUE ~ 0),
           sidd = case_when(cluster == "SIDD" ~ 1,
                            TRUE ~ 0),
           sird = case_when(cluster == "SIRD" ~ 1,
                            TRUE ~ 0))
  
  tdcm_df <- cluster_df %>%
    arrange(study, study_id, age) %>%
    group_by(study, study_id) %>%
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
                            TRUE ~ race))

  overall_tdcm[[i]] <- coxph(Surv(tstart, tstop, event) ~ study + female + race + min_age + bmi + hba1c + homa2b 
                             + homa2ir + ldlc + sbp + egfr_ckdepi_2021, 
                             data = tdcm_df, weights = ipcw_cluster)
  
  mard_tdcm[[i]] <- coxph(Surv(tstart, tstop, mard) ~ study + female + race + min_age + bmi + hba1c + homa2b 
                          + homa2ir + ldlc + sbp + egfr_ckdepi_2021, 
                          data = tdcm_df, weights = ipcw_cluster)
  
  mod_tdcm[[i]] <- coxph(Surv(tstart, tstop, mod) ~ study + female + race + min_age + bmi + hba1c + homa2b 
                         + homa2ir + ldlc + sbp + egfr_ckdepi_2021, 
                         data = tdcm_df, weights = ipcw_cluster)
  
  sidd_tdcm[[i]] <- coxph(Surv(tstart, tstop, sidd) ~ study + female + race + min_age + bmi + hba1c + homa2b 
                          + homa2ir + ldlc + sbp + egfr_ckdepi_2021, 
                          data = tdcm_df, weights = ipcw_cluster)
  
  sird_tdcm[[i]] <- coxph(Surv(tstart, tstop, sird) ~ study + female + race + min_age + bmi + hba1c + homa2b 
                          + homa2ir + ldlc + sbp + egfr_ckdepi_2021, 
                          data = tdcm_df, weights = ipcw_cluster)
}



source("functions/pool_results.R")

tdcm_results <- bind_rows(
  pool_results(overall_tdcm) %>% mutate(model = "Overall"),
  pool_results(mard_tdcm) %>% mutate(model = "MARD"),
  pool_results(mod_tdcm) %>% mutate(model = "MOD"),
  pool_results(sidd_tdcm) %>% mutate(model = "SIDD"),
  pool_results(sird_tdcm) %>% mutate(model = "SIRD")) %>% 
  write_csv(.,"analysis/dspan03_tdcm pooled results with multiple imputation.csv")


