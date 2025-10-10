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
  dplyr::filter(dpp_intervention == 1) %>% 
  distinct(study,study_id,dpp_intervention) # n = 60


coxph_dfs <- list()
train_sets <- list()
test_sets <- list()
overall_coxph <- list()
mard_coxph <- list()
mod_coxph <- list()
sidd_coxph <- list()
sird_coxph <- list()

###################### BASELINE ############################

set.seed(123)

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
  
  
  # Stratified train/test split based on subtype
  train_df <- coxph_df %>%
    group_by(subtype) %>%
    sample_frac(0.75) %>%
    ungroup()
  
  test_df <- anti_join(coxph_df, train_df, by = "joint_id")
  
  train_sets[[i]] <- train_df
  test_sets[[i]] <- test_df
  
  
  
  # Cox PH models
  
  overall_coxph[[i]] <- coxph(Surv(time_to_event, newdm_event) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                              + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention, 
                              data = train_df)
  
  mard_coxph[[i]] <- coxph(Surv(time_to_event, mard) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention, 
                           data = train_df)
  
  mod_coxph[[i]] <- coxph(Surv(time_to_event, mod) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                          + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention, 
                          data = train_df)
  
  sidd_coxph[[i]] <- coxph(Surv(time_to_event, sidd) ~ study + female + race3 + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + ridge(dpp_intervention, theta = 100),
                           data = train_df, ties = "efron", control = coxph.control(iter.max = 100))
  
  # CARDIA has 0 people from SIRD
  sird_coxph[[i]] <- coxph(Surv(time_to_event, sird) ~ strata(study) + female + race3 + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + ridge(dppos_interv, theta = 50), 
                           data = train_df, ties = "efron", control = coxph.control(iter.max = 200, eps = 1e-09))
  
}


# Validate on 25% test set
validation_results <- map2_dfr(1:mi_dfs$m, test_sets, function(i, test_df) {
  tibble(
    imputation = i,
    overall_c_index = concordance(overall_coxph[[i]], newdata = test_df)$concordance,
    mard_c_index    = concordance(mard_coxph[[i]],    newdata = test_df)$concordance,
    mod_c_index     = concordance(mod_coxph[[i]],     newdata = test_df)$concordance,
    sidd_c_index    = concordance(sidd_coxph[[i]],    newdata = test_df)$concordance,
    sird_c_index    = concordance(sird_coxph[[i]],    newdata = test_df)$concordance
  )
}) %>% 
  rename_with(~ c("imputation", "Overall", "MARD", "MOD", "SIDD", "SIRD")) %>%         # Rename columns
  mutate(across(Overall:SIRD, ~ round(.x, 3))) 

# Save as CSV
write.csv(validation_results, "analysis/dspan07_coxph validation cindex.csv", row.names = FALSE)


###################### RANDOM VISIT ############################

set.seed(42)

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
    # one random visit per person
    group_by(joint_id) %>%
    dplyr::slice_sample(n = 1) %>%
    ungroup() %>%
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
  
  
  # Stratified train/test split based on subtype
  train_df <- coxph_df %>%
    group_by(subtype) %>%
    sample_frac(0.75) %>%
    ungroup()
  
  test_df <- anti_join(coxph_df, train_df, by = "joint_id")
  
  train_sets[[i]] <- train_df
  test_sets[[i]] <- test_df
  
  
  
  # Cox PH models
  
  overall_coxph[[i]] <- coxph(Surv(time_to_event, newdm_event) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                              + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention, 
                              data = train_df)
  
  mard_coxph[[i]] <- coxph(Surv(time_to_event, mard) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention, 
                           data = train_df)
  
  mod_coxph[[i]] <- coxph(Surv(time_to_event, mod) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                          + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention, 
                          data = train_df)
  
  sidd_coxph[[i]] <- coxph(Surv(time_to_event, sidd) ~ study + female + race3 + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + ridge(dpp_intervention, theta = 100),
                           data = train_df, ties = "efron", control = coxph.control(iter.max = 100))
  
  # CARDIA has 0 people from SIRD
  sird_coxph[[i]] <- coxph(Surv(time_to_event, sird) ~ strata(study) + female + race3 + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + ridge(dppos_interv, theta = 50), 
                           data = train_df, ties = "efron", control = coxph.control(iter.max = 200, eps = 1e-09))
  
}


# Validate on 25% test set
validation_results <- map2_dfr(1:mi_dfs$m, test_sets, function(i, test_df) {
  tibble(
    imputation = i,
    overall_c_index = concordance(overall_coxph[[i]], newdata = test_df)$concordance,
    mard_c_index    = concordance(mard_coxph[[i]],    newdata = test_df)$concordance,
    mod_c_index     = concordance(mod_coxph[[i]],     newdata = test_df)$concordance,
    sidd_c_index    = concordance(sidd_coxph[[i]],    newdata = test_df)$concordance,
    sird_c_index    = concordance(sird_coxph[[i]],    newdata = test_df)$concordance
  )
}) %>% 
  rename_with(~ c("imputation", "Overall", "MARD", "MOD", "SIDD", "SIRD")) %>%         # Rename columns
  mutate(across(Overall:SIRD, ~ round(.x, 3))) 

# Save as CSV
write.csv(validation_results, "analysis/dspan07_coxph validation cindex radom visit.csv", row.names = FALSE)


  