rm(list = ls());gc();source(".Rprofile")

library(survival)
library(survminer)
library(ggsurvfit)
library(broom)
library(tidyr)
library(lme4)

ipcw_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/ipcw_dfs.RDS"))

overall_cp <- list()
mard_cp <- list()
mod_cp <- list()
sidd_cp <- list()
sird_cp <- list()

# Cox PH model
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
  
  # Cox PH - baseline data
  cross_df <- cluster_df %>% 
    group_by(study_id,study) %>% 
    dplyr::filter(age == min(age)) %>% 
    ungroup()
  
  overall_cp[[i]] <- coxph(Surv(time_to_event, event) ~ study + race + female + age + bmi + hba1c + homa2b + homa2ir 
                           + ldlc + sbp + egfr_ckdepi_2021, 
                           data = cross_df, weights = ipcw_cluster)
  
  mard_cp[[i]] <- coxph(Surv(time_to_event, mard) ~ study + race + female + age + bmi + hba1c + homa2b + homa2ir 
                        + ldlc + sbp + egfr_ckdepi_2021, 
                        data = cross_df, weights = ipcw_cluster)
  
  mod_cp[[i]] <- coxph(Surv(time_to_event, mod) ~ study + race + female + age + bmi + hba1c + homa2b + homa2ir 
                       + ldlc + sbp + egfr_ckdepi_2021, 
                       data = cross_df, weights = ipcw_cluster)
  
  sidd_cp[[i]] <- coxph(Surv(time_to_event, sidd) ~ study + race + female + age + bmi + hba1c + homa2b + homa2ir 
                        + ldlc + sbp + egfr_ckdepi_2021, 
                        data = cross_df, weights = ipcw_cluster)
  
  sird_cp[[i]] <- coxph(Surv(time_to_event, sird) ~ study + race + female + age + bmi + hba1c + homa2b + homa2ir 
                        + ldlc + sbp + egfr_ckdepi_2021, 
                        data = cross_df, weights = ipcw_cluster)
  
  
}


# Pooling coefficients ------------
source("functions/clean_mi_contrasts.R")

overall_cp_out = clean_mi_contrasts(overall_cp,link="coxph")











#--------------------------------------------------------------------------------------------------------------------
# TDCM - longitudinal data

overall_tdcm <- list()
mard_tdcm <- list()
mod_tdcm <- list()
sidd_tdcm <- list()
sird_tdcm <- list()

for (i in 1:length(ipcw_dfs)) {
  df <- ipcw_dfs[[1]]  
  
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
    # dplyr::filter(tstart < tstop)
    dplyr::filter((tstart < tstop) & (tstop <= censored_age))

  
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



#--------------------------------------------------------------------------------------------------------------------
# mixed effect model

overall_mix <- list()
mard_mix <- list()
mod_mix <- list()
sidd_mix <- list()
sird_mix <- list()

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
  
  
  overall_mix[[i]] <- glmer(event ~ time_to_event + (1|study) + female + race + min_age + bmi + hba1c + homa2b 
                             + homa2ir + ldlc + sbp + egfr_ckdepi_2021, 
                             data = cluster_df, weights = ipcw_cluster, family = binomial(link = "logit"))
  
  mard_mix[[i]] <- glmer(mard ~ time_to_event + (1|study) + female + race + min_age + bmi + hba1c + homa2b 
                         + homa2ir + ldlc + sbp + egfr_ckdepi_2021, 
                         data = cluster_df, weights = ipcw_cluster, family = binomial(link = "logit"))
  
  mod_mix[[i]] <- glmer(mod ~ time_to_event + (1|study) + female + race + min_age + bmi + hba1c + homa2b 
                        + homa2ir + ldlc + sbp + egfr_ckdepi_2021, 
                        data = cluster_df, weights = ipcw_cluster, family = binomial(link = "logit"))
  
  sidd_mix[[i]] <- glmer(sidd ~ time_to_event + (1|study) + female + race + min_age + bmi + hba1c + homa2b 
                         + homa2ir + ldlc + sbp + egfr_ckdepi_2021, 
                         data = cluster_df, weights = ipcw_cluster, family = binomial(link = "logit"))
  
  sird_mix[[i]] <- glmer(sird ~ time_to_event + (1|study) + female + race + min_age + bmi + hba1c + homa2b 
                         + homa2ir + ldlc + sbp + egfr_ckdepi_2021, 
                         data = cluster_df, weights = ipcw_cluster, family = binomial(link = "logit"))
  
  
}



overall_mix_out = clean_mi_contrasts(overall_mix,link="glmer log")




















