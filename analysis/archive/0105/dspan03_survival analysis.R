rm(list = ls());gc();source(".Rprofile")

library(survival)
library(survminer)
library(ggsurvfit)
library(broom)

ipcw_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan02_ipcw dfs.RDS"))

overall_cp <- list()
mard_cp <- list()
mod_cp <- list()
sidd_cp <- list()
sird_cp <- list()
cp_results <- list()
cp_output <- list()
df0 <- list()
df1 <- list()
df2 <- list()
df3 <- list()
df4 <- list()
coxph_output <- list()

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
    ungroup() %>% 
    mutate(across(c(bmi, hba1c, homa2b, homa2ir, ldlc, sbp, egfr_ckdepi_2021), ~replace(., is.infinite(.), NA))) %>% 
    mutate(race = case_when(race == "NH Other" ~ "Other", 
                            TRUE ~ race))
  
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
  
  
  cp_results[[i]] <- bind_rows(
    broom::tidy(overall_cp[[i]]) %>% mutate(model = "Overall"),
    broom::tidy(mard_cp[[i]]) %>% mutate(model = "MARD"),
    broom::tidy(mod_cp[[i]]) %>% mutate(model = "MOD"),
    broom::tidy(sidd_cp[[i]]) %>% mutate(model = "SIDD"),
    broom::tidy(sird_cp[[i]]) %>% mutate(model = "SIRD")) 
  
  # covert to Hazard Ratio
  cp_output[[i]] <- cp_results[[i]] %>% 
    mutate(HR = exp(estimate),
           lci = exp(estimate - 1.96 * std.error),
           uci = exp(estimate + 1.96 * std.error)) %>% 
    mutate(coef_ci = paste0(round(HR, 2), " (", round(lci, 2), ", ", round(uci, 2), ")")) %>% 
    pivot_wider(names_from = model, values_from = coef_ci) %>% 
    dplyr::select(term, Overall, MARD, MOD, SIDD, SIRD) 
  
  df0[[i]] <- cp_output[[i]] %>% 
    dplyr::select(term, Overall) 
  df1[[i]] <- cp_output[[i]] %>% 
    dplyr::select(term, MARD) 
  df2[[i]] <- cp_output[[i]] %>% 
    dplyr::select(term, MOD)
  df3[[i]] <- cp_output[[i]] %>% 
    dplyr::select(term, SIDD)
  df4[[i]] <- cp_output[[i]] %>% 
    dplyr::select(term, SIRD)
  
  coxph_output[[i]] <- na.omit(df0[[i]]) %>% 
    left_join(na.omit(df1[[i]]), by = "term") %>%
    left_join(na.omit(df2[[i]]), by = "term") %>% 
    left_join(na.omit(df3[[i]]), by = "term") %>% 
    left_join(na.omit(df4[[i]]), by = "term") %>% 
    mutate(model = paste0("m", i))
  
}

coxph_output_results <- bind_rows(coxph_output) %>% 
  write_csv(.,"analysis/dspan03_cox ph with multiple imputation.csv")


#--------------------------------------------------------------------------------------------------------------------
# TDCM - longitudinal data

overall_tdcm <- list()
mard_tdcm <- list()
mod_tdcm <- list()
sidd_tdcm <- list()
sird_tdcm <- list()


D = length(ipcw_dfs)

for (i in 1:D) {
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
    # dplyr::filter(tstart < tstop)
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


