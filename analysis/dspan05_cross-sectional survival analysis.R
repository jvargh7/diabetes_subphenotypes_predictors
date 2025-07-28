rm(list = ls());gc();source(".Rprofile")

library(mice)
library(survival)
library(survminer)
library(ggsurvfit)
library(broom)
library(cmprsk)
library(nnet)
source("functions/egfr_ckdepi_2021.R")

mi_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/mi_dfs_new.RDS"))

clean_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>% 
  dplyr::filter(dpp_intervention == 1) %>% 
  distinct(study,study_id,dpp_intervention) # n = 60

analytic_dfs <- list()
overall_coxph <- list()
mard_coxph <- list()
mod_coxph <- list()
sidd_coxph <- list()
sird_coxph <- list()

mod_sdh <- mard_sdh <- sidd_sdh <- sird_sdh <- list()

multinom_all <- list()
multinom_ref <- list()

ph_tests_overall <- vector("list", length(overall_coxph))
ph_tests_mard    <- vector("list", length(mard_coxph))
ph_tests_mod     <- vector("list", length(mod_coxph))
ph_tests_sidd    <- vector("list", length(sidd_coxph))
ph_tests_sird    <- vector("list", length(sird_coxph))


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
    )
    ) %>%
    ungroup() 
  
  analytic_dfs[[i]] <- analytic_df
}


for (i in 1:length(analytic_dfs)) {
  df <- analytic_dfs[[i]]  
  
  cluster_df <- analytic_df %>% 
    mutate(mard = case_when(cluster == "MARD" ~ 1,
                            TRUE ~ 0),
           mod = case_when(cluster == "MOD" ~ 1,
                           TRUE ~ 0),
           sidd = case_when(cluster == "SIDD" ~ 1,
                            TRUE ~ 0),
           sird = case_when(cluster == "SIRD" ~ 1,
                            TRUE ~ 0)) %>% 
    mutate(cluster_numeric = case_when(cluster == "MARD" ~ 1,
                                       cluster == "MOD" ~ 2,
                                       cluster == "SIDD" ~ 3,
                                       cluster == "SIRD" ~ 4,
                                       TRUE ~ 5),
           race_white = case_when(race == "NH White" ~ 1,
                                  TRUE ~ 0),
           race_black = case_when(race == "NH Black" ~ 1,
                                     TRUE ~ 0),
           race_hispanic = case_when(race == "Hispanic" ~ 1,
                                     TRUE ~ 0),
           race_other = case_when(race %in% c("NH Other", "Other") ~ 1,
                                  TRUE ~ 0),
           study_cardia = case_when (study == "cardia" ~ 1,
                                     TRUE ~ 0),
           study_dppos = case_when(study == "dppos" ~ 1,
                                   TRUE ~ 0),
           study_jhs = case_when(study == "jhs" ~ 1,
                                 TRUE ~ 0),
           study_mesa= case_when(study == "mesa" ~ 1,
                                 TRUE ~ 0)
           )
  
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
    mutate(subtype = factor(subtype, levels=c("NOT2D","MARD","MOD","SIDD","SIRD")))
  
  # Multinomial models
  multinom_all[[i]] <- nnet::multinom(subtype ~ log(time_to_event) + study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                                      + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention,
                                      data = coxph_df)
  
  coxph_ref <- coxph_df %>% 
    mutate(subtype <- relevel(subtype, ref = "MARD"))
  # reference: MARD
  multinom_ref[[i]] <- nnet::multinom(subtype ~ log(time_to_event) + study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                                      + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention,
                                      data = coxph_ref)
  
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
  
  # SDH competing risk models (crr)
  mard_sdh[[i]] <- crr(coxph_df$time_to_event, coxph_df$cluster_numeric, cengroup = 0, failcode = 1,
                       cov1 = coxph_df[, c("study_cardia", "study_dppos", "study_jhs", "female", "race_black", "race_hispanic", "race_other", "earliest_age",
                                           "bmi", "hba1c", "homa2b_scaled", "homa2ir", "ldlc_scaled",
                                           "sbp_scaled", "egfr_ckdepi_2021_scaled", "dpp_intervention")])
  
  mod_sdh[[i]] <- crr(coxph_df$time_to_event, coxph_df$cluster_numeric, cengroup = 0, failcode = 2,
                      cov1 = coxph_df[, c("study_cardia", "study_dppos", "study_jhs", "female", "race_black", "race_hispanic", "race_other", "earliest_age",
                                          "bmi", "hba1c", "homa2b_scaled", "homa2ir", "ldlc_scaled",
                                          "sbp_scaled", "egfr_ckdepi_2021_scaled", "dpp_intervention")])
  
  sidd_sdh[[i]] <- crr(coxph_df$time_to_event, coxph_df$cluster_numeric, cengroup = 0, failcode = 3,
                       cov1 = coxph_df[, c("study_cardia", "study_dppos", "study_jhs", "female", "race_black", "race_hispanic", "race_other", "earliest_age",
                                           "bmi", "hba1c", "homa2b_scaled", "homa2ir", "ldlc_scaled",
                                           "sbp_scaled", "egfr_ckdepi_2021_scaled", "dpp_intervention")])
  
  sird_sdh[[i]] <- crr(coxph_df$time_to_event, coxph_df$cluster_numeric, cengroup = 0, failcode = 4,
                       cov1 = coxph_df[, c("study_dppos", "study_jhs", "female", "race_black", "race_hispanic", "race_other", "earliest_age",
                                           "bmi", "hba1c", "homa2b_scaled", "homa2ir", "ldlc_scaled",
                                           "sbp_scaled", "egfr_ckdepi_2021_scaled")])
  
  # PH assumption tests
  ph_tests_overall[[i]] <- cox.zph(overall_coxph[[i]])
  ph_tests_mard[[i]]    <- cox.zph(mard_coxph[[i]])
  ph_tests_mod[[i]]     <- cox.zph(mod_coxph[[i]])
  ph_tests_sidd[[i]]    <- cox.zph(sidd_coxph[[i]])
  ph_tests_sird[[i]]    <- cox.zph(sird_coxph[[i]])
  
  
}


# COXPH ------------
unpooled_resuts_coxph <- bind_rows(
  imap_dfr(overall_coxph,function(m,index){
    broom::tidy(m) %>% mutate(model = "Overall",iteration = index) %>% 
      return(.)
  }),
  imap_dfr(mard_coxph,function(m,index){
    broom::tidy(m) %>% mutate(model = "MARD",iteration = index) %>% 
      return(.)
  }),
  imap_dfr(mod_coxph,function(m,index){
    broom::tidy(m) %>% mutate(model = "MOD",iteration = index) %>% 
      return(.)
  }),
  imap_dfr(sidd_coxph,function(m,index){
    broom::tidy(m) %>% mutate(model = "SIDD",iteration = index) %>% 
      return(.)
  }),
  imap_dfr(sird_coxph,function(m,index){
    broom::tidy(m) %>% mutate(model = "SIRD",iteration = index) %>% 
      return(.)
  })
  
  
  
  
)

write_csv(unpooled_resuts_coxph, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan05_unpooled coxph results.csv"))



source("functions/pool_results.R")

coxph_results <- bind_rows(
  pool_results(overall_coxph) %>% mutate(model = "Overall"),
  pool_results(mard_coxph) %>% mutate(model = "MARD"),
  pool_results(mod_coxph) %>% mutate(model = "MOD"),
  pool_results(sidd_coxph) %>% mutate(model = "SIDD"),
  pool_results(sird_coxph) %>% mutate(model = "SIRD")
) %>% 
  write_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan05_pooled cox ph results.csv"))



# create Schoenfeld residual plots
# save_cox_zph_plots <- function(ph_test, file_prefix = "model", width = 1200, height = 1200, res = 150) {
#   # Number of covariates (exclude "GLOBAL")
#   n_covariates <- nrow(ph_test$table) - 1
#   n_row <- ceiling(sqrt(n_covariates))
#   n_col <- ceiling(n_covariates / n_row)
#   
#   png(paste0(file_prefix, ".png"), width = width, height = height, res = res)
#   par(mfrow = c(n_row, n_col))
#   plot(ph_test)
#   par(mfrow = c(1, 1))
#   dev.off()
# }
# 
# i <- 1  # choose the imputation index
# 
# save_cox_zph_plots(ph_tests_overall[[i]], "cox_zph_overall")
# save_cox_zph_plots(ph_tests_mard[[i]],    "cox_zph_mard")
# save_cox_zph_plots(ph_tests_mod[[i]],     "cox_zph_mod")
# save_cox_zph_plots(ph_tests_sidd[[i]],    "cox_zph_sidd")
# save_cox_zph_plots(ph_tests_sird[[i]],    "cox_zph_sird")
# 
# # for each model
# n_imputations <- length(ph_tests_overall)
# 
# for (i in seq_len(n_imputations)) {
#   save_cox_zph_plots(ph_tests_overall[[i]], paste0("cox_zph_overall_imp", i))
#   save_cox_zph_plots(ph_tests_mard[[i]],    paste0("cox_zph_mard_imp", i))
#   save_cox_zph_plots(ph_tests_mod[[i]],     paste0("cox_zph_mod_imp", i))
#   save_cox_zph_plots(ph_tests_sidd[[i]],    paste0("cox_zph_sidd_imp", i))
#   save_cox_zph_plots(ph_tests_sird[[i]],    paste0("cox_zph_sird_imp", i))
# }


# Fine-Gray SDH ------------
unpooled_resuts_sdh <- bind_rows(
  imap_dfr(mard_sdh,function(m,index){
    broom::tidy(m) %>% mutate(model = "MARD",iteration = index) %>% 
      return(.)
  }),
  imap_dfr(mod_sdh,function(m,index){
    broom::tidy(m) %>% mutate(model = "MOD",iteration = index) %>% 
      return(.)
  }),
  imap_dfr(sidd_sdh,function(m,index){
    broom::tidy(m) %>% mutate(model = "SIDD",iteration = index) %>% 
      return(.)
  }),
  imap_dfr(sird_sdh,function(m,index){
    broom::tidy(m) %>% mutate(model = "SIRD",iteration = index) %>% 
      return(.)
  })
  
  
  
  
)

write_csv(unpooled_resuts_sdh, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan05_unpooled sdh results.csv"))


# Fine-Gray Models similar to Cox PH models 
# Aalen-Johnson Estimator similar to Kaplan-Meier Estimator (unadjusted)
# https://statisticalhorizons.com/for-causal-analysis-of-competing-risks/
# https://www.ahajournals.org/doi/10.1161/CIRCOUTCOMES.121.008368
source("functions/pool_results.R")

sdh_results <- bind_rows(
  pool_results(mard_sdh) %>% mutate(model = "MARD"),
  pool_results(mod_sdh) %>% mutate(model = "MOD"),
  pool_results(sidd_sdh) %>% mutate(model = "SIDD"),
  pool_results(sird_sdh) %>% mutate(model = "SIRD")
) %>% 
  write_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan05_pooled sdh results.csv"))


# Multinomial models -----------

multinom_results <- pool_results(multinom_all) 
multinom_ref_results <- pool_results(multinom_ref) 

bind_rows(multinom_results %>% mutate(reference = "NOT2D"),
          multinom_ref_results %>% mutate(reference = "MARD")) %>% 
  write_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan05_pooled multinomial results.csv"))
