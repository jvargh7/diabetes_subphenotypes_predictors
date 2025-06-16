rm(list = ls());gc();source(".Rprofile")

source("functions/egfr_ckdepi_2021.R")

# included participants (n = 10,309)
selected_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>% 
  distinct(joint_id)

# new dm with non-missing age at diagnosis, N = 7,623
final_dataset_temp = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/final_dataset_temp.RDS")) %>%
  mutate(joint_id = paste(study, original_study_id, sep = "_"))

# new T2D (n = 7,623) + No T2D (n = 22,373)
analytic_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dsppre01_analytic df.RDS")) %>% 
  mutate(egfr_ckdepi_2021 = egfr_ckdepi_2021(scr = serumcreatinine,female = female,age = age),
         joint_id = paste(study, study_id, sep = "_"),
         # newly diagnosed T2D or no T2D
         newdm_event = case_when(joint_id %in% final_dataset_temp$joint_id ~ 1,
                                 is.na(dmagediag) ~ 0,
                                 TRUE ~ NA)) %>%
  mutate(study = case_when(study == "dpp" ~ "dppos",
                           TRUE ~ study)) %>% 
  group_by(joint_id) %>%
  mutate(min_age = min(age, na.rm = TRUE)) %>%
  ungroup() %>% 
  dplyr::filter(!is.na(newdm_event)) %>% 
  mutate(included = case_when(joint_id %in% selected_df$joint_id ~ 1,
                              TRUE ~ 0),
         race = case_when(race == "NH Other" ~ "Other", 
                          TRUE ~ race),
         race = factor(race, levels=c("Hispanic","NH Black","NH White","Other")),
         dpp_intervention = case_when(is.na(dpp_intervention) ~ 0,
                                      TRUE ~ 1))

baseline_df <- analytic_df %>% 
  group_by(joint_id) %>% 
  slice_min(order_by = age, with_ties = FALSE) %>% 
  ungroup()

# Fit logistic regression predicting inclusion
ipw_model <- glm(included ~ min_age + female + race + dpp_intervention,
                 data = baseline_df,
                 family = binomial())

p_included <- mean(baseline_df$included == 1)


ipw_df <- baseline_df %>%
  mutate(
    prob_included = predict(ipw_model, newdata = baseline_df, type = "response"),
    ipw_weight = if_else(included == 1, 1 / prob_included, NA_real_),
    # reduce the influence of extreme weights
    ipw_stabilized = if_else(included == 1, p_included / prob_included, NA_real_)
  ) %>% 
  dplyr::filter(!is.na(ipw_stabilized)) %>% 
  distinct(joint_id,ipw_stabilized)


saveRDS(ipw_df, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspse01a_ipw df.RDS"))

