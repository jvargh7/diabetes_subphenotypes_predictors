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



############  C-INDEX ############

source("functions/model_evaluate_index.R")

bootstrap_cindex_results <- list()

event_list <- c("sird", "sidd", "mod", "mard", "event_true")

# Step 1: Collect all C-index values from each imputation for each event
all_bootstrap_values <- list()

for (event_var in event_list) {
  all_cindex <- c()
  
  for (i in 1:length(tdcm_dfs)) {
    df <- tdcm_dfs[[i]]
    unique_ids <- unique(df$joint_id)
    
    boot_result <- boot(
      data = unique_ids,
      statistic = function(ids, i) {
        sampled_ids <- ids[i]
        subset_data <- df[df$joint_id %in% sampled_ids, ]
        c_index_function(subset_data, 1:length(sampled_ids), event_var)
      },
      R = 200
    )
    
    all_cindex <- c(all_cindex, na.omit(boot_result$t))
  }
  
  all_bootstrap_values[[event_var]] <- all_cindex
}

# Step 2: Summarize pooled C-index results
final_cindex_summary <- do.call(rbind, lapply(names(all_bootstrap_values), function(event_name) {
  values <- all_bootstrap_values[[event_name]]
  data.frame(
    event = event_name,
    mean_cindex = round(mean(values, na.rm = TRUE), 4),
    CI_lower = round(quantile(values, 0.025, na.rm = TRUE), 4),
    CI_upper = round(quantile(values, 0.975, na.rm = TRUE), 4)
  )
}))


########### SENSITIVITY & SPECIFICITY ########### 

source("functions/model_evaluate_index.R")

# Define events
event_list <- c("sird", "sidd", "mod", "mard", "event_true")

# Initialize storage
pooled_sens_spec_results <- data.frame()

# Loop over events
set.seed(123)
for (event_var in event_list) {
  pooled_sens <- c()
  pooled_spec <- c()
  pooled_f1 <- c()
  
  # Loop over imputations
  for (i in 1:length(tdcm_dfs)) {
    df <- tdcm_dfs[[i]]
    unique_ids <- unique(df$joint_id)
    
    boot_result <- boot(
      data = unique_ids,
      statistic = function(ids, i) {
        sampled_ids <- ids[i]
        subset_data <- df[df$joint_id %in% sampled_ids, ]
        sens_spec_function(subset_data, 1:length(sampled_ids), event_var)
      },
      R = 200
    )
    
    pooled_sens <- c(pooled_sens, na.omit(boot_result$t[, 1]))
    pooled_spec <- c(pooled_spec, na.omit(boot_result$t[, 2]))
    pooled_f1   <- c(pooled_f1,   na.omit(boot_result$t[, 3]))
  }
  
  # Summarize across all imputations
  pooled_sens_spec_results <- rbind(pooled_sens_spec_results, data.frame(
    event = event_var,
    mean_sensitivity = round(mean(pooled_sens, na.rm = TRUE), 4),
    CI_sens_lower = round(quantile(pooled_sens, 0.025, na.rm = TRUE), 4),
    CI_sens_upper = round(quantile(pooled_sens, 0.975, na.rm = TRUE), 4),
    mean_specificity = round(mean(pooled_spec, na.rm = TRUE), 4),
    CI_spec_lower = round(quantile(pooled_spec, 0.025, na.rm = TRUE), 4),
    CI_spec_upper = round(quantile(pooled_spec, 0.975, na.rm = TRUE), 4),
    mean_f1 = round(mean(pooled_f1, na.rm = TRUE), 4),
    CI_f1_lower = round(quantile(pooled_f1, 0.025, na.rm = TRUE), 4),
    CI_f1_upper = round(quantile(pooled_f1, 0.975, na.rm = TRUE), 4)
  ))
}


library(glue)

final_table <- merge(final_cindex_summary, pooled_sens_spec_results, by = "event") %>% 
  mutate(
    c_index = glue("{round(mean_cindex, 3)} ({round(CI_lower, 3)}, {round(CI_upper, 3)})"),
    sensitivity = glue("{round(mean_sensitivity, 3)} ({round(CI_sens_lower, 3)}, {round(CI_sens_upper, 3)})"),
    specificity = glue("{round(mean_specificity, 3)} ({round(CI_spec_lower, 3)}, {round(CI_spec_upper, 3)})"),
    f1_score = glue("{round(mean_f1, 3)} ({round(CI_f1_lower, 3)}, {round(CI_f1_upper, 3)})")
  ) %>% 
  select(event,c_index,sensitivity,specificity,f1_score)

write.csv(final_table, "analysis/dspan06_pooled model evaluation index.csv", row.names = FALSE)










############## calibration curve (need to be fixed later) ##############

bootstrap_cal_slope <- function(data, indices, event_name) {
  sampled_ids <- unique(data$joint_id)[indices]
  d <- data[data$joint_id %in% sampled_ids, ]
  
  formula_str <- as.formula(paste0("Surv(tstart, tstop, ", event_name, ") ~ study + female + race + earliest_age + bmi + hba1c + ",
                                   "homa2b_scaled + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention"))
  fit <- try(coxph(formula_str, data = d), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA)
  
  d$lp <- predict(fit, type = "lp")
  cal_model <- try(coxph(Surv(tstart, tstop, d[[event_name]]) ~ lp, data = d), silent = TRUE)
  if (inherits(cal_model, "try-error")) return(NA)
  
  return(coef(cal_model)[1])
}


source("functions/model_evaluate_index.R")

library(rms)

event_list <- c("sird", "sidd", "mod", "mard", "event_true")
pooled_calibration_slopes <- data.frame()


set.seed(123)
results <- data.frame()

for (event_name in event_list) {
  slopes_all <- c()
  
  for (df in tdcm_dfs) {
    unique_ids <- unique(df$joint_id)
    
    boot_obj <- boot(
      data = unique_ids,
      statistic = function(ids, i) bootstrap_cal_slope(df, i, event_name),
      R = 200
    )
    
    slopes <- na.omit(boot_obj$t)
    slopes_all <- c(slopes_all, slopes)
  }
  
  results <- rbind(results, data.frame(
    event = event_name,
    mean_slope = round(mean(slopes_all), 4),
    CI_lower = round(quantile(slopes_all, 0.025), 4),
    CI_upper = round(quantile(slopes_all, 0.975), 4)
  ))
}


# Save or print result
write.csv(pooled_calibration_slopes, "analysis/pooled calibration slopes.csv", row.names = FALSE)
print(pooled_calibration_slopes)
