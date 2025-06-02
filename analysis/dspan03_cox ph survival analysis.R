rm(list = ls());gc();source(".Rprofile")

library(mice)
library(survival)
library(survminer)
library(ggsurvfit)
library(broom)

source("functions/egfr_ckdepi_2021.R")
mi_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/mi_dfs.RDS"))

analytic_dfs <- list()
overall_coxph <- list()
mard_coxph <- list()
mod_coxph <- list()
sidd_coxph <- list()
sird_coxph <- list()

ph_tests_overall <- vector("list", length(overall_coxph))
ph_tests_mard    <- vector("list", length(mard_coxph))
ph_tests_mod     <- vector("list", length(mod_coxph))
ph_tests_sidd    <- vector("list", length(sidd_coxph))
ph_tests_sird    <- vector("list", length(sird_coxph))


for(i in 1:mi_dfs$m) {
  df <- complete(mi_dfs, action = i) %>% 
    mutate(egfr_ckdepi_2021 = egfr_ckdepi_2021(scr = serumcreatinine,female = female,age = age),
           time_to_event = censored_age - age) 
  
  
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
                            TRUE ~ 0))
  
  coxph_df <- cluster_df %>%
    dplyr::filter(age == earliest_age) %>% 
    mutate(across(c(bmi, hba1c, homa2b, homa2ir, ldlc, sbp, egfr_ckdepi_2021), ~replace(., is.infinite(.), NA))) %>% 
    # error due to 0 ppl in NH Other (sidd == 1), ignore this category
    mutate(race = case_when(race == "NH Other" ~ "Other", 
                            TRUE ~ race)) %>% 
    # scaling
    mutate(sbp_scaled = sbp/10,
           ldlc_scaled = ldlc/10,
           homa2b_scaled = homa2b/10,
           egfr_ckdepi_2021_scaled = egfr_ckdepi_2021/10)
  
  
  # Cox PH
  
  overall_coxph[[i]] <- coxph(Surv(time_to_event, newdm_event) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                              + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled, 
                              data = coxph_df)
  
  mard_coxph[[i]] <- coxph(Surv(time_to_event, mard) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled, 
                           data = coxph_df)
  
  mod_coxph[[i]] <- coxph(Surv(time_to_event, mod) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                          + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled, 
                          data = coxph_df)
  
  sidd_coxph[[i]] <- coxph(Surv(time_to_event, sidd) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled, 
                           data = coxph_df)
  
  # CARDIA has 0 people from SIRD 
  coxph_sird <- coxph_df %>% dplyr::filter(!(study == "cardia"))
  sird_coxph[[i]] <- coxph(Surv(time_to_event, sird) ~ study + female + race + earliest_age + bmi + hba1c + homa2b_scaled 
                           + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled, 
                           data = coxph_sird)
  
  
  ph_tests_overall[[i]] <- cox.zph(overall_coxph[[i]])

  ph_tests_mard[[i]]    <- cox.zph(mard_coxph[[i]])

  ph_tests_mod[[i]]     <- cox.zph(mod_coxph[[i]])

  ph_tests_sidd[[i]]    <- cox.zph(sidd_coxph[[i]])

  ph_tests_sird[[i]]    <- cox.zph(sird_coxph[[i]])
  
}


source("functions/pool_results.R")

coxph_results <- bind_rows(
  pool_results(overall_coxph) %>% mutate(model = "Overall"),
  pool_results(mard_coxph) %>% mutate(model = "MARD"),
  pool_results(mod_coxph) %>% mutate(model = "MOD"),
  pool_results(sidd_coxph) %>% mutate(model = "SIDD"),
  pool_results(sird_coxph) %>% mutate(model = "SIRD")
) %>% 
  write_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan03_pooled cox ph results.csv"))



# create Schoenfeld residual plots
save_cox_zph_plots <- function(ph_test, file_prefix = "model", width = 1200, height = 1200, res = 150) {
  # Number of covariates (exclude "GLOBAL")
  n_covariates <- nrow(ph_test$table) - 1
  n_row <- ceiling(sqrt(n_covariates))
  n_col <- ceiling(n_covariates / n_row)
  
  png(paste0(file_prefix, ".png"), width = width, height = height, res = res)
  par(mfrow = c(n_row, n_col))
  plot(ph_test)
  par(mfrow = c(1, 1))
  dev.off()
}

i <- 1  # choose the imputation index

save_cox_zph_plots(ph_tests_overall[[i]], "cox_zph_overall")
save_cox_zph_plots(ph_tests_mard[[i]],    "cox_zph_mard")
save_cox_zph_plots(ph_tests_mod[[i]],     "cox_zph_mod")
save_cox_zph_plots(ph_tests_sidd[[i]],    "cox_zph_sidd")
save_cox_zph_plots(ph_tests_sird[[i]],    "cox_zph_sird")

# for each model
n_imputations <- length(ph_tests_overall)

for (i in seq_len(n_imputations)) {
  save_cox_zph_plots(ph_tests_overall[[i]], paste0("cox_zph_overall_imp", i))
  save_cox_zph_plots(ph_tests_mard[[i]],    paste0("cox_zph_mard_imp", i))
  save_cox_zph_plots(ph_tests_mod[[i]],     paste0("cox_zph_mod_imp", i))
  save_cox_zph_plots(ph_tests_sidd[[i]],    paste0("cox_zph_sidd_imp", i))
  save_cox_zph_plots(ph_tests_sird[[i]],    paste0("cox_zph_sird_imp", i))
}
