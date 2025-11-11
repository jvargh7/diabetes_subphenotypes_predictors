rm(list = ls());gc();source(".Rprofile")

library(mice)
library(survival)
library(survAUC)
library(pROC)
library(Metrics)
library(purrr)

source("functions/egfr_ckdepi_2021.R")

mi_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dsphyc301_mi_dfs.RDS"))

clean_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>% 
  select(joint_id,age,dpp_intervention,smoking,med_chol_use,med_bp_use,med_dep_use) %>% 
  mutate(dpp_intervention = case_when(
    dpp_intervention == 1 ~ 1,
    TRUE ~ 0
  ))


coxph_dfs <- list()
overall_coxph <- list()
mard_coxph <- list()
mod_coxph <- list()
sidd_coxph <- list()
sird_coxph <- list()


###################### BASELINE ############################

for(i in 1:mi_dfs$m) {
  df <- complete(mi_dfs, action = i) %>% 
    rename(joint_id = original_joint_id) %>%
    mutate(egfr_ckdepi_2021 = egfr_ckdepi_2021(scr = serumcreatinine,female = female,age = age),
           time_to_event = censored_age - age) %>% 
    left_join(clean_df,
              by = c("joint_id","age")) %>% 
    mutate(smoking = case_when(is.na(smoking) ~ "Never",
                               TRUE ~ smoking),
           med_bp_use = case_when(is.na(med_bp_use) ~ 0,
                                  TRUE ~ med_bp_use),
           med_chol_use = case_when(is.na(med_chol_use) ~ 0,
                                    TRUE ~ med_chol_use),
           med_dep_use = case_when(is.na(med_dep_use) ~ 0,
                                   TRUE ~ med_dep_use))
  
  
  analytic_df <- df %>% 
    arrange(study,study_id,joint_id,age) %>% 
    group_by(study,study_id,joint_id) %>%
    ungroup() 
  

  cluster_df <- analytic_df %>% 
    mutate(mard = case_when(cluster == "MARD" ~ 1,
                            TRUE ~ 0),
           mod = case_when(cluster ==  "MOD" ~ 1,
                           TRUE ~ 0),
           sidd = case_when(cluster == "SIDD" ~ 1,
                            TRUE ~ 0),
           sird = case_when(cluster == "SIRD" ~ 1,
                            TRUE ~ 0)) 
  
  coxph_df <- cluster_df %>%
    dplyr::filter(age == earliest_age) %>% 
    # For baseline analysis, use the original newdm_event as the outcome
    mutate(event = newdm_event) %>%
    mutate(across(c(bmi, hba1c, homa2b, homa2ir, ldlc, sbp, egfr_ckdepi_2021), ~replace(., is.infinite(.), NA))) %>% 
    # error due to 0 ppl in NH Other (sidd == 1), ignore this category
    mutate(race = case_when(race == "NH Other" ~ "Other", 
                            TRUE ~ race),
           race = relevel(factor(race), ref = "NH White"),
           race3 = fct_collapse(race,
                                "NH Black" = "NH Black",
                                "Hispanic" = "Hispanic",
                                "White/Other" = c("NH White","Other")),
           race3 = relevel(race3, ref = "White/Other")) %>% 
    # scaling
    mutate(sbp_scaled = sbp/10,
           ldlc_scaled = ldlc/10,
           homa2b_scaled = homa2b/10,
           egfr_ckdepi_2021_scaled = egfr_ckdepi_2021/10) %>% 
    mutate(subtype = case_when(is.na(cluster) ~ "NOT2D",
                               TRUE ~ cluster)) %>%
    mutate(subtype = factor(subtype, levels=c("NOT2D","MARD","MOD","SIDD","SIRD")),
           # 1 if in DPPOS intervention arm, else 0 (also 0 for all non-DPPOS rows)
           dppos_interv = as.numeric(study == "dppos") * as.numeric(dpp_intervention))
  
  coxph_dfs[[i]] <- coxph_df
  
  
  # Cox PH models
  
  overall_coxph[[i]] <- coxph(Surv(time_to_event, newdm_event) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                              + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention
                              + smoking + med_chol_use + med_bp_use + med_dep_use, 
                              data = coxph_df)
  
  mard_coxph[[i]] <- coxph(Surv(time_to_event, mard) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention
                           + smoking + med_chol_use + med_bp_use + med_dep_use, 
                           data = coxph_df)
  
  mod_coxph[[i]] <- coxph(Surv(time_to_event, mod) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                          + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention
                          + smoking + med_chol_use + med_bp_use + med_dep_use, 
                          data = coxph_df)
  
  sidd_coxph[[i]] <- coxph(Surv(time_to_event, sidd) ~ study + female + race3 + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + ridge(dpp_intervention, theta = 100)
                           + smoking + med_chol_use + med_bp_use + med_dep_use,
                           data = coxph_df, ties = "efron", control = coxph.control(iter.max = 100))
  
  # CARDIA has 0 people from SIRD
  sird_coxph[[i]] <- coxph(Surv(time_to_event, sird) ~ strata(study) + female + race3 + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + ridge(dppos_interv, theta = 50)
                           + smoking + med_chol_use + med_bp_use + med_dep_use, 
                           data = coxph_df, ties = "efron", control = coxph.control(iter.max = 200, eps = 1e-09))
}



source("functions/evaluate_coxph_quantile.R")


# Prepare list to store evaluation metrics for each model
model_names <- c("overall_coxph", "mard_coxph", "mod_coxph", "sidd_coxph", "sird_coxph")
pooled_results <- setNames(vector("list", length(model_names)), model_names)

time_horizon <- 8  # Use 8-year risk for baseline analysis (median follow-up time)

for (model_name in model_names) {
  pooled_metrics <- map2(1:10, coxph_dfs, function(i, df) {
    model_list <- get(model_name)
    model <- model_list[[i]]
    # For baseline analysis, the event is whether diabetes occurred within the time horizon
    # Keep the original event definition but ensure we're using the right logic
    evaluate_coxph_quantile(model, df, time_horizon, quantiles = c(0.10, 0.25, 0.50, 0.75, 0.90))
  })
  
  # Stack all 10 imputations
  pooled_df <- bind_rows(pooled_metrics, .id = "imputation")
  pooled_summary <- pooled_df %>%
    group_by(threshold) %>%
    summarise(across(c(threshold_value, sensitivity, specificity, F1, c_index, cal_slope), 
                     list(mean = mean, sd = sd), na.rm = TRUE), .groups = "drop")
  
  pooled_results[[model_name]] <- pooled_summary
}


pooled_long <- purrr::imap_dfr(pooled_results, ~mutate(.x, model = .y)) %>% 
  mutate(
    threshold_value = sprintf("%.3f (%.3f)", threshold_value_mean, threshold_value_sd),
    sensitivity = sprintf("%.3f (%.3f)", sensitivity_mean, sensitivity_sd),
    specificity = sprintf("%.3f (%.3f)", specificity_mean, specificity_sd),
    f1          = sprintf("%.3f (%.3f)", F1_mean, F1_sd),
    c_index     = sprintf("%.3f (%.3f)", c_index_mean, c_index_sd),
    cal_slope   = sprintf("%.3f (%.3f)", cal_slope_mean, cal_slope_sd)
  ) %>%
  select(model, threshold, threshold_value, sensitivity, specificity, f1, c_index, cal_slope)

write.csv(pooled_long, "analysis/dspan06_pooled coxph model metrics.csv", row.names = FALSE)


###################### RANDOM VISIT ############################

set.seed(42) 

for(i in 1:mi_dfs$m) {
  df <- complete(mi_dfs, action = i) %>% 
    rename(joint_id = original_joint_id) %>%
    mutate(egfr_ckdepi_2021 = egfr_ckdepi_2021(scr = serumcreatinine,female = female,age = age),
           time_to_event = censored_age - age) %>% 
    left_join(clean_df,
              by = c("joint_id","age")) %>% 
    mutate(smoking = case_when(is.na(smoking) ~ "Never",
                               TRUE ~ smoking),
           med_bp_use = case_when(is.na(med_bp_use) ~ 0,
                                  TRUE ~ med_bp_use),
           med_chol_use = case_when(is.na(med_chol_use) ~ 0,
                                    TRUE ~ med_chol_use),
           med_dep_use = case_when(is.na(med_dep_use) ~ 0,
                                   TRUE ~ med_dep_use))
  
  
  analytic_df <- df %>% 
    arrange(study,study_id,joint_id,age) %>% 
    group_by(study,study_id,joint_id) %>%
    ungroup() 
  
  
  cluster_df <- analytic_df %>% 
    mutate(mard = case_when(cluster == "MARD" ~ 1,
                            TRUE ~ 0),
           mod = case_when(cluster ==  "MOD" ~ 1,
                           TRUE ~ 0),
           sidd = case_when(cluster == "SIDD" ~ 1,
                            TRUE ~ 0),
           sird = case_when(cluster == "SIRD" ~ 1,
                            TRUE ~ 0)) 
  
  coxph_df <- cluster_df %>%
    # one random visit per person
    group_by(joint_id) %>%
    dplyr::slice_sample(n = 1) %>%
    ungroup() %>% 
    # For random visit analysis, use the original newdm_event as the outcome (consistent with baseline)
    mutate(event = newdm_event) %>%
    mutate(across(c(bmi, hba1c, homa2b, homa2ir, ldlc, sbp, egfr_ckdepi_2021), ~replace(., is.infinite(.), NA))) %>% 
    # error due to 0 ppl in NH Other (sidd == 1), ignore this category
    mutate(race = case_when(race == "NH Other" ~ "Other", 
                            TRUE ~ race),
           race = relevel(factor(race), ref = "NH White"),
           race3 = fct_collapse(race,
                                "NH Black" = "NH Black",
                                "Hispanic" = "Hispanic",
                                "White/Other" = c("NH White","Other")),
           race3 = relevel(race3, ref = "White/Other")) %>% 
    # scaling
    mutate(sbp_scaled = sbp/10,
           ldlc_scaled = ldlc/10,
           homa2b_scaled = homa2b/10,
           egfr_ckdepi_2021_scaled = egfr_ckdepi_2021/10) %>% 
    mutate(subtype = case_when(is.na(cluster) ~ "NOT2D",
                               TRUE ~ cluster)) %>%
    mutate(subtype = factor(subtype, levels=c("NOT2D","MARD","MOD","SIDD","SIRD")),
           # 1 if in DPPOS intervention arm, else 0 (also 0 for all non-DPPOS rows)
           dppos_interv = as.numeric(study == "dppos") * as.numeric(dpp_intervention))
  
  
  coxph_dfs[[i]] = coxph_df
  
  # Cox PH models
  
  overall_coxph[[i]] <- coxph(Surv(time_to_event, newdm_event) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                              + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention
                              + smoking + med_chol_use + med_bp_use + med_dep_use, 
                              data = coxph_df)
  
  mard_coxph[[i]] <- coxph(Surv(time_to_event, mard) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention
                           + smoking + med_chol_use + med_bp_use + med_dep_use, 
                           data = coxph_df)
  
  mod_coxph[[i]] <- coxph(Surv(time_to_event, mod) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                          + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention
                          + smoking + med_chol_use + med_bp_use + med_dep_use, 
                          data = coxph_df)
  
  sidd_coxph[[i]] <- coxph(Surv(time_to_event, sidd) ~ study + female + race3 + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + ridge(dpp_intervention, theta = 100)
                           + smoking + med_chol_use + med_bp_use + med_dep_use,
                           data = coxph_df, ties = "efron", control = coxph.control(iter.max = 100))
  
  # CARDIA has 0 people from SIRD
  sird_coxph[[i]] <- coxph(Surv(time_to_event, sird) ~ strata(study) + female + race3 + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + ridge(dppos_interv, theta = 50)
                           + smoking + med_chol_use + med_bp_use + med_dep_use, 
                           data = coxph_df, ties = "efron", control = coxph.control(iter.max = 200, eps = 1e-09))
}



source("functions/evaluate_coxph_quantile.R")

model_names <- c("overall_coxph", "mard_coxph", "mod_coxph", "sidd_coxph", "sird_coxph")
pooled_results <- setNames(vector("list", length(model_names)), model_names)
names(pooled_results) <- model_names

time_horizon <- 8  # Use 8-year risk for random visit analysis (consistent with baseline)

M <- length(coxph_dfs)

for (model_name in model_names) {
  model_list <- get(model_name, inherits = TRUE)
  if (!is.list(model_list) || length(model_list) != M) {
    stop(sprintf("`%s` must be a list of %d models (one per imputation).", model_name, M))
  }
  
  pooled_metrics <- map2(seq_len(M), coxph_dfs, function(i, df) {
    model <- model_list[[i]]
    # Use the same quantile approach as baseline for consistency
    evaluate_coxph_quantile(model, df, time_horizon, quantiles = c(0.10, 0.25, 0.50, 0.75, 0.90))
  })
  
  pooled_df <- bind_rows(pooled_metrics, .id = "imputation")
  pooled_summary <- pooled_df %>%
    group_by(threshold) %>%
    summarise(across(c(threshold_value, sensitivity, specificity, F1, c_index, cal_slope), 
                     list(mean = mean, sd = sd), na.rm = TRUE), .groups = "drop")
  
  pooled_results[[model_name]] <- pooled_summary
}


pooled_long <- purrr::imap_dfr(pooled_results, ~mutate(.x, model = .y)) %>% 
  mutate(
    threshold_value = sprintf("%.3f (%.3f)", threshold_value_mean, threshold_value_sd),
    sensitivity = sprintf("%.3f (%.3f)", sensitivity_mean, sensitivity_sd),
    specificity = sprintf("%.3f (%.3f)", specificity_mean, specificity_sd),
    f1          = sprintf("%.3f (%.3f)", F1_mean, F1_sd),
    c_index     = sprintf("%.3f (%.3f)", c_index_mean, c_index_sd),
    cal_slope   = sprintf("%.3f (%.3f)", cal_slope_mean, cal_slope_sd)
  ) %>%
  select(model, threshold, threshold_value, sensitivity, specificity, f1, c_index, cal_slope)

write.csv(pooled_long, "analysis/dspan06_pooled coxph model metrics random visit.csv", row.names = FALSE)
