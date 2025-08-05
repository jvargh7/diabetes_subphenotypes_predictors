rm(list = ls());gc();source(".Rprofile")

library(mice)
library(survival)
library(survminer)
library(ggsurvfit)
library(broom)
library(boot)

source("functions/egfr_ckdepi_2021.R")
mi_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/mi_dfs_new.RDS"))

clean_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>% 
  dplyr::filter(dpp_intervention == 1) %>% 
  distinct(study,study_id,joint_id,dpp_intervention) # n = 60


tdcm_dfs <- list()


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
}


################# SENSITIVITY, SPECIFICITY, F1 SCORE ###########################

# get predicted probabilities, concordance ---- takes more than 12h
source("functions/fit_predict_boot_model.R")

library(boot)
library(arrow)

event_list <- c("event_true", "mod", "mard", "sidd", "sird")

for (event_name in event_list) {
  all_results <- list()
  
  for (m in 1:length(tdcm_dfs)) {
    df <- tdcm_dfs[[m]]
    unique_ids <- unique(df$joint_id)
    
    for (b in 1:200) {
      set.seed(1000 + m * 200 + b)
      sampled_ids <- sample(unique_ids, length(unique_ids), replace = TRUE)
      
      res <- fit_predict_boot_model(df, which(unique_ids %in% sampled_ids), event_name, m, b)
      if (!is.null(res)) {
        all_results[[length(all_results) + 1]] <- res
      }
    }
  }
  
  combined_df <- do.call(rbind, all_results)
  write_parquet(combined_df, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan06_", event_name, "_df.parquet"))
}

# compute metrics ----------------------------------------------

t2d_df <- read_parquet(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan06_event_true_df.parquet"))
mard_df <- read_parquet(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan06_mard_df.parquet"))
mod_df <- read_parquet(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan06_mod_df.parquet"))
sidd_df <- read_parquet(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan06_sidd_df.parquet"))
sird_df <- read_parquet(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan06_sird_df.parquet"))


# Function to compute metrics at multiple cutoffs
source("functions/evaluate_metrics_with_ci.R")

metrics_t2d <- evaluate_metrics_with_ci(t2d_df, "T2D")
metrics_mard <- evaluate_metrics_with_ci(mard_df, "MARD")
metrics_mod  <- evaluate_metrics_with_ci(mod_df, "MOD")
metrics_sidd <- evaluate_metrics_with_ci(sidd_df, "SIDD")
metrics_sird <- evaluate_metrics_with_ci(sird_df, "SIRD")

all_cutoff_metrics <- bind_rows(metrics_t2d, metrics_mard, metrics_mod, metrics_sidd, metrics_sird) %>% 
  mutate(
    sensitivity = glue("{round(sensitivity, 3)} ({round(sens_lower, 3)}, {round(sens_upper, 3)})"),
    specificity = glue("{round(specificity, 3)} ({round(spec_lower, 3)}, {round(spec_upper, 3)})"),
    f1_score    = glue("{round(f1, 3)} ({round(f1_lower, 3)}, {round(f1_upper, 3)})"),
    c_index     = glue("{round(c_index, 3)} ({round(c_index_lower, 3)}, {round(c_index_upper, 3)})")
  ) %>%
  select(event, cutoff, sensitivity, specificity, f1_score, c_index) %>% 
  write.csv(., "analysis/dspan06_metrics at multiple cutoffs.csv", row.names = FALSE)



################# CALIBRATION SLOPE ###########################

source("functions/get_calibration_slope.R")

cal_slope_results <- data.frame()

set.seed(123)
for (m in 1:length(tdcm_dfs)) {
  df <- tdcm_dfs[[m]]
  unique_ids <- unique(df$joint_id)
  
  for (event_name in c("event_true", "mod", "mard", "sidd", "sird")) {
    for (b in 1:200) {
      sampled_ids <- sample(unique_ids, length(unique_ids), replace = TRUE)
      d_boot <- df[df$joint_id %in% sampled_ids, ]
      
      slope <- get_calibration_slope(d_boot, event_name)
      
      cal_slope_results <- rbind(cal_slope_results, data.frame(
        m = m,
        b = b,
        event = event_name,
        calibration_slope = slope
      ))
    }
  }
}

pooled_slope <- cal_slope_results %>%
  group_by(event) %>%
  summarize(mean_slope = round(mean(calibration_slope, na.rm = TRUE), 3),
            lower = round(quantile(calibration_slope, 0.025, na.rm = TRUE), 3),
            upper = round(quantile(calibration_slope, 0.975, na.rm = TRUE), 3)) %>% 
  mutate(slope = glue("{round(mean_slope, 3)} ({round(lower, 3)}, {round(upper, 3)})")) %>% 
  select(event, slope) %>% 
  write.csv(., "analysis/dspan06_pooled calibration slopes.csv", row.names = FALSE)

