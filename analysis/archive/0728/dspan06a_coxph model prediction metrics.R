rm(list = ls());gc();source(".Rprofile")

library(mice)
library(survival)
library(survAUC)
library(pROC)
library(Metrics)
library(purrr)

source("functions/egfr_ckdepi_2021.R")

mi_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/mi_dfs_new.RDS"))

clean_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>% 
  dplyr::filter(dpp_intervention == 1) %>% 
  distinct(study,study_id,dpp_intervention) # n = 60


coxph_dfs <- list()
overall_coxph <- list()
mard_coxph <- list()
mod_coxph <- list()
sidd_coxph <- list()
sird_coxph <- list()


for(i in 1:mi_dfs$m) {
  df <- complete(mi_dfs, action = i) %>% 
    rename(joint_id = original_joint_id) %>%
    mutate(egfr_ckdepi_2021 = egfr_ckdepi_2021(scr = serumcreatinine,female = female,age = age),
           time_to_event = censored_age - age) %>% 
    left_join(clean_df,
              by = c("study","study_id")) %>% 
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
    )) %>%
    ungroup() 
  

  cluster_df <- analytic_df %>% 
    mutate(mard = case_when(cluster == "MARD" ~ 1,
                            TRUE ~ 0),
           mod = case_when(cluster == "MOD" ~ 1,
                           TRUE ~ 0),
           sidd = case_when(cluster == "SIDD" ~ 1,
                            TRUE ~ 0),
           sird = case_when(cluster == "SIRD" ~ 1,
                            TRUE ~ 0)) 
  
  coxph_df <- cluster_df %>%
    dplyr::filter(age == earliest_age) %>% 
    mutate(across(c(bmi, hba1c, homa2b, homa2ir, ldlc, sbp, egfr_ckdepi_2021), ~replace(., is.infinite(.), NA))) %>% 
    # error due to 0 ppl in NH Other (sidd == 1), ignore this category
    mutate(race = case_when(race == "NH Other" ~ "Other", 
                            TRUE ~ race),
           race = relevel(factor(race), ref = "NH White")) %>% 
    # scaling
    mutate(sbp_scaled = sbp/10,
           ldlc_scaled = ldlc/10,
           homa2b_scaled = homa2b/10,
           egfr_ckdepi_2021_scaled = egfr_ckdepi_2021/10) %>% 
    mutate(subtype = case_when(is.na(cluster) ~ "NOT2D",
                               TRUE ~ cluster)) %>%
    mutate(subtype = factor(subtype, levels=c("NOT2D","MARD","MOD","SIDD","SIRD")),
           study_cardia = case_when (study == "cardia" ~ 1,
                                     TRUE ~ 0),
           study_dppos = case_when(study == "dppos" ~ 1,
                                   TRUE ~ 0),
           study_jhs = case_when(study == "jhs" ~ 1,
                                 TRUE ~ 0),
           study_mesa= case_when(study == "mesa" ~ 1,
                                 TRUE ~ 0))
  
  coxph_dfs[[i]] <- coxph_df
  
  
  # Cox PH models
  
  overall_coxph[[i]] <- coxph(Surv(time_to_event, newdm_event) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                              + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention, 
                              data = coxph_df)
  
  mard_coxph[[i]] <- coxph(Surv(time_to_event, mard) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention, 
                           data = coxph_df)
  
  mod_coxph[[i]] <- coxph(Surv(time_to_event, mod) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                          + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention, 
                          data = coxph_df)
  
  sidd_coxph[[i]] <- coxph(Surv(time_to_event, sidd) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention, 
                           data = coxph_df)
  
  # CARDIA has 0 people from SIRD
  # sird_coxph uses study-specific dummies
  sird_coxph[[i]] <- coxph(Surv(time_to_event, sird) ~ study_dppos + study_jhs + study_cardia + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled, 
                           data = coxph_df)
}



source("functions/evaluate_coxph.R")


# Prepare list to store evaluation metrics for each model
model_names <- c("overall_coxph", "mard_coxph", "mod_coxph", "sidd_coxph", "sird_coxph")
pooled_results <- setNames(vector("list", length(model_names)), model_names)

time_horizon <- 5  # Example: 5-year risk

for (model_name in model_names) {
  pooled_metrics <- map2(1:10, coxph_dfs, function(i, df) {
    model_list <- get(model_name)
    model <- model_list[[i]]
    eval_df <- df %>% mutate(event = ifelse(newdm_event == 1 & time_to_event <= time_horizon, 1, 0))
    evaluate_model(model, eval_df, time_horizon)
  })
  
  # Stack all 10 imputations
  pooled_df <- bind_rows(pooled_metrics, .id = "imputation")
  pooled_summary <- pooled_df %>%
    group_by(threshold) %>%
    summarise(across(c(sensitivity, specificity, F1, c_index, cal_slope), 
                     list(mean = mean, sd = sd), na.rm = TRUE), .groups = "drop")
  
  pooled_results[[model_name]] <- pooled_summary
}


pooled_long <- purrr::imap_dfr(pooled_results, ~mutate(.x, model = .y)) %>% 
  mutate(
    sensitivity = sprintf("%.3f (%.3f)", sensitivity_mean, sensitivity_sd),
    specificity = sprintf("%.3f (%.3f)", specificity_mean, specificity_sd),
    f1          = sprintf("%.3f (%.3f)", F1_mean, F1_sd),
    c_index     = sprintf("%.3f (%.3f)", c_index_mean, c_index_sd),
    cal_slope   = sprintf("%.3f (%.3f)", cal_slope_mean, cal_slope_sd)
  ) %>%
  select(model, threshold, sensitivity, specificity, f1, c_index, cal_slope)

write.csv(pooled_long, "analysis/dspan06a_pooled coxph model metrics.csv", row.names = FALSE)


