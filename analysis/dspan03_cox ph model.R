rm(list = ls());gc();source(".Rprofile")

library(survival)
library(survminer)
library(ggsurvfit)
library(broom)
library(tidyr)

analytic_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01_analytic df.RDS")) %>% 
  group_by(study_id) %>% 
  mutate(min_age = min(age),
         max_age = max(age)) %>%
  dplyr::filter(age <= min_age + 15) %>%
  select(-max_age) %>% 
  ungroup() %>% 
  mutate(dm = case_when(!is.na(dmagediag) | hba1c >= 6.5 | glucosef >= 126 | glucose2h >= 200 ~ 1,
                        TRUE ~ 0))
# dm: 8103, non-dm: 17446
# MARD: 1386, MOD: 1040, SIDD: 99, SIRD: 370

time_df <- analytic_df %>% 
  group_by(study_id) %>% 
  mutate(lab_diagage = case_when(hba1c >= 6.5 | glucosef >= 126 | glucose2h >= 200 ~ age, 
                                 TRUE ~ NA_real_),
         min_lab_diagage = min(lab_diagage),
         min_dmagediag = min(dmagediag)) %>% 
  
  mutate(time_to_event = case_when(
    !is.na(min_dmagediag) ~ abs(min_dmagediag - min_age),
    is.na(min_dmagediag) & !is.na(min_lab_diagage) ~ abs(min_lab_diagage - min_age),
    TRUE ~ NA_real_
  )) %>% 
  ungroup() %>% 
  select(-c(min_age,min_lab_diagage,min_dmagediag))

#------------------------------------------------------------------------------------------------------------------------
cox_mod <- coxph(Surv(time_to_event, dm) ~ study + age + female + race_eth + bmi + hba1c + sbp + dbp + hdlc + ldlc + homa2b + homa2ir, data = time_df)
summary(cox_mod)

ggsurvplot(survfit(cox_mod, data = time_df), color = "#2E9FDF",
           ggtheme = theme_minimal())


cox_model <- coxph(Surv(time_to_event, dm) ~ age + female + race_eth + bmi + hba1c + sbp + dbp + hdlc + ldlc + homa2b + homa2ir + strata(study), 
                   data = time_df)

summary(cox_model)

#------------------------------------------------------------------------------------------------------------------------
# cause-specific HR
analytic_sample <- time_df %>% 
  mutate(mard = case_when(cluster == "MARD" ~ 1,
                          TRUE ~ 0),
         mod = case_when(cluster == "MOD" ~ 1,
                          TRUE ~ 0),
         sidd = case_when(cluster == "SIDD" ~ 1,
                          TRUE ~ 0),
         sird = case_when(cluster == "SIRD" ~ 1,
                          TRUE ~ 0))

diseases <- c("mard", "mod", "sidd", "sird")

# overall rates
analytic_sample %>% 
  summarize(across(diseases, ~mean(.)))

# Define the regression function with survey weights
regression_dm <- function(outcome_var, df) {
  
  # Cox models with survey weights
  m1 <- coxph(as.formula(paste0("Surv(time_to_event, ", outcome_var, ") ~ study + female + race_eth")),
              data = df)
  
  m2 <- coxph(as.formula(paste0("Surv(time_to_event, ", outcome_var, ") ~ study + female + race_eth + bmi + hba1c + sbp + dbp + hdlc + ldlc + homa2b + homa2ir")),
              data = df)
  m0 <- coxph(as.formula(paste0("Surv(time_to_event, ", outcome_var, ") ~ study")),
              data = df)
  
  # Combine model results
  bind_rows(
            broom::tidy(m0) %>% mutate(model = "m0"),
            broom::tidy(m1) %>% mutate(model = "m1"),
            broom::tidy(m2) %>% mutate(model = "m2")) %>% 
    mutate(outcome = outcome_var) %>% 
    return(.)
}

# Run regression for each cluster and save results
regression_results <- list()
for (i in seq_along(diseases)) {
  print(i)
  regression_results[[diseases[i]]] <- regression_dm(outcome_var = diseases[i],
                                                                df = analytic_sample)
}

# Save regression results to CSV
regression_results %>% 
  bind_rows() %>% 
  write_csv("analysis/dspan03_survival analysis.csv")


# Format regression results for viewing (Age- and sex-adjusted model)
regression_results %>% 
  bind_rows() %>% 
  mutate(HR = exp(estimate),
         lci = exp(estimate - 1.96 * std.error),
         uci = exp(estimate + 1.96 * std.error),
         model = case_when(model == "m1" ~ "Sex- and race-adjusted",
                           model == "m2" ~ "Sex-, race- and lab-adjusted",
                           TRUE ~ NA_character_)) %>% 
  # dplyr::filter(model == "Age- and sex-adjusted",
  #               str_detect(term, "cluster")) %>% 
  # mutate(cluster = str_replace(term, "cluster", "")) %>% 
  mutate(coef_ci = paste0(round(HR, 2), " (", round(lci, 2), ", ", round(uci, 2), ")"),
         outcome = factor(outcome, levels = diseases, labels = diseases)) %>% 
  dplyr::select(term, outcome, coef_ci) %>% 
  # pivot_wider(names_from = cluster, values_from = coef_ci) %>% 
  # mutate(MARD = "Ref") %>% 
  # dplyr::select(outcome, MARD, MOD, SIDD, SIRD) %>% 
  # View() %>% 
  write_csv("analysis/dspan03_survival analysis results.csv")

























