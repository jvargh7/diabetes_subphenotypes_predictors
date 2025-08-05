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
  distinct(study,study_id,joint_id,dpp_intervention) # n = 60


tdcm_dfs <- list()
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
    )) %>%
    ungroup() 
  
  
  tdcm_df <- analytic_df %>%
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
  
  tdcm_dfs[[i]] <- tdcm_df
  
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



# ================= Landmark evaluation ==================

landmark_age <- 50
prediction_horizon <- 10
cutpoints <- c(0.10, 0.25, 0.50, 0.75, 0.90)
landmark_metrics_list <- list()

for (i in 1:10) {
  model_list <- list(
    overall = overall_tdcm[[i]],
    mard = mard_tdcm[[i]],
    mod = mod_tdcm[[i]],
    sidd = sidd_tdcm[[i]],
    sird = sird_tdcm[[i]]
  )
  
  df <- tdcm_dfs[[i]]
  
  # Select observation closest to landmark age
  landmark_df <- df %>%
    mutate(age_diff = abs(tstart - landmark_age)) %>%
    group_by(joint_id) %>%
    slice_min(order_by = age_diff, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      fu_time = censored_age - tstart,
      event_within_5yr = ifelse(event_true == 1 & fu_time <= prediction_horizon, 1, 0),
      followup_enough = fu_time >= prediction_horizon | event_within_5yr == 1
    ) %>%
    dplyr::filter(followup_enough)
  
  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    lp <- predict(model, newdata = landmark_df, type = "lp")
    bh <- basehaz(model, centered = FALSE)
    S0 <- approx(bh$time, exp(-bh$hazard), xout = prediction_horizon, rule = 2)$y
    risk <- 1 - S0^exp(lp)
    
    truth <- landmark_df$event_within_5yr
    cal_model <- glm(truth ~ lp, family = binomial)
    cal_slope <- coef(cal_model)[2]
    c_index <- concordance(model)$concordance
    
    metrics_df <- lapply(cutpoints, function(thresh) {
      pred_class <- ifelse(risk >= thresh, 1, 0)
      TP <- sum(pred_class == 1 & truth == 1)
      TN <- sum(pred_class == 0 & truth == 0)
      FP <- sum(pred_class == 1 & truth == 0)
      FN <- sum(pred_class == 0 & truth == 1)
      
      sensitivity <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
      specificity <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
      f1 <- ifelse((2 * TP + FP + FN) == 0, NA, 2 * TP / (2 * TP + FP + FN))
      
      data.frame(
        threshold = thresh,
        sensitivity = sensitivity,
        specificity = specificity,
        f1 = f1,
        c_index = c_index,
        cal_slope = cal_slope,
        model = model_name,
        imputation = i
      )
    }) %>% bind_rows()
    
    landmark_metrics_list[[paste0(model_name, "_", i)]] <- metrics_df
  }
}

landmark_metrics <- bind_rows(landmark_metrics_list)

pooled_metrics <- landmark_metrics %>%
  group_by(model, threshold) %>%
  summarise(across(c(sensitivity, specificity, f1, c_index, cal_slope),
                   list(mean = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE))),
            .groups = "drop")

write.csv(landmark_metrics, "analysis/dspan06b_landmark metrics all.csv", row.names = FALSE)
write.csv(pooled_metrics, "analysis/dspan06b_landmark metrics pooled.csv", row.names = FALSE)


