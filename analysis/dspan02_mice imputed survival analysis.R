rm(list = ls());gc();source(".Rprofile")

library(survival)
library(survminer)
library(ggsurvfit)
library(broom)
library(tidyr)

source("functions/egfr_ckdepi_2021.R")

# follow-up time >15y & dm == 0 --- event == 0
# mice_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/mi_dfs.RDS")) %>% 
imputed_df <- read_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dsppre02_knn imputation.csv")) %>%
  mutate(egfr_ckdepi_2021 = egfr_ckdepi_2021(scr = serumcreatinine, female = female, age = age)) %>% 
  group_by(study, study_id) %>% 
  mutate(min_age = min(age)) %>%
  # Restrict to observations within 15 years of earliest age
  dplyr::filter(age <= (min_age + 15)) %>%
  # There are individuals whose dmagediag < minimum age
  dplyr::filter(is.na(dmagediag) | (age <= dmagediag)) %>% 
  
  mutate(event = case_when(# Individuals who are never diagnosed
                           is.na(dmagediag) ~ 0,
                           
                           # Individuals who are diagnosed within 15 years of earliest wave
                           dmagediag <= (min_age + 15) ~ 1,
                           
                           # Individuals who are diagnosed after 15 years of earliest wave
                           TRUE ~ 0
                           
                           
                           ),
         censored_age = case_when(is.na(dmagediag) ~ max(age),
                                  
                                  !is.na(dmagediag) & dmagediag <= (min_age + 15) ~ dmagediag,
                                  
                                  TRUE ~ max(age)
                                  
                                  )) %>%
  ungroup() %>% 
  # 5 categories in study: aric, cardia, dppos, jhs, mesa
  mutate(study_aric = case_when(study == "aric" ~ 1,
                                TRUE ~ 0),
         study_cardia = case_when(study == "cardia" ~ 1,
                                TRUE ~ 0),
         study_dppos = case_when(study == "dppos" ~ 1,
                                TRUE ~ 0),
         study_jhs = case_when(study == "jhs" ~ 1,
                                TRUE ~ 0),
    # 5 categories in race: NH White, NH Black, NH Other, Hispanic, Other  
         race_nhwhi = case_when(race == "NH White" ~ 1,
                                TRUE ~ 0),
         race_nhbla = case_when(race == "NH Black" ~ 1,
                                TRUE ~ 0),
         race_nhoth = case_when(race == "NH Other" ~ 1,
                                TRUE ~ 0),
         race_hisp = case_when(race == "Hispanic" ~ 1,
                                TRUE ~ 0)) %>%
  # define time to event
  mutate(time_to_event = censored_age - age) %>% 
  dplyr::filter(time_to_event > 0) 

#------------------------------------------------------------------------------------------------------------------------
# cross-sectional datatset, cox PH model
# overall, include all variables, 1 obs for each person
cross_df <- imputed_df %>% 
  group_by(study_id,study) %>% 
  dplyr::filter(age == min(age)) %>% 
  ungroup()

table(cross_df$event)

# Visualize distribution of time to event or censoring
ggplot(data=cross_df,aes(x=time_to_event,group=event,fill=factor(event))) +
  geom_histogram(alpha=0.5,position=position_dodge(width=0.9),bins=10)


cox_mod <- coxph(Surv(time_to_event, event) ~ study_aric + study_cardia + study_dppos + study_jhs + age + female 
                                           + race_nhwhi + race_nhbla + race_nhoth + race_hisp + bmi + hba1c + sbp + dbp + hdlc + ldlc 
                                           + homa2b + homa2ir + height + wc + insulinf + glucosef + glucose2h + tgl
                                           + serumcreatinine + urinecreatinine + egfr + apo_a + apo_b + uric_acid 
                                           + vldlc + hc + triceps + iliac + abdominal + medial + ast + alt + urinealbumin
                                           + weight, data = cross_df)

cox_mod_var8 <- coxph(Surv(time_to_event, event) ~ age + bmi + hba1c + homa2b + homa2ir + ldlc + sbp + egfr_ckdepi_2021, 
                      data = cross_df)

summary(cox_mod)
summary(cox_mod_var8)
# p>.05 covariates: insulinf, abdominal

#------------------------------------------------------------------------------------------------------------------------
# cause-specific HR
analytic_sample <- imputed_df %>% 
  ### exclude diabetes patients with no cluster data
  dplyr::filter(!(event == 1 & is.na(cluster))) %>% 
  mutate(mard = case_when(cluster == "MARD" ~ 1,
                          TRUE ~ 0),
         mod = case_when(cluster == "MOD" ~ 1,
                         TRUE ~ 0),
         sidd = case_when(cluster == "SIDD" ~ 1,
                          TRUE ~ 0),
         sird = case_when(cluster == "SIRD" ~ 1,
                          TRUE ~ 0))

saveRDS(analytic_sample, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan02_cluster df.RDS"))

diseases <- c("mard", "mod", "sidd", "sird")


# overall rates
analytic_sample %>% 
  summarize(across(diseases, ~mean(.)))

#--------------------------------------------------------------------------------------------------------------------------
# variable selection - LASSO, SCAD
# chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cran.r-project.org/web/packages/ncvreg/ncvreg.pdf

# exclude missing values
variables_to_check <- c("study_aric", "study_cardia", "study_dppos", "study_jhs", 
                        "race_nhwhi", "race_nhbla", "race_nhoth", "race_hisp",
                        "age", "female", "bmi", "hba1c", "sbp", "dbp", "hdlc", "ldlc",
                        "homa2b", "homa2ir", "height", "wc", "glucosef", "glucose2h", "tgl",
                        "serumcreatinine", "urinecreatinine", "egfr_ckdepi_2021", "apo_a", "apo_b", "uric_acid",
                        "vldlc", "hc", "triceps", "iliac", "medial", "ast", "alt", "urinealbumin",
                        "weight", "insulinf", "abdominal")


df_clean <- analytic_sample %>%
  dplyr::filter(complete.cases(.[variables_to_check])) %>%
  # require all variables to be numeric
  mutate(female = as.numeric(as.character(female)))

y <- Surv(df_clean$time_to_event, df_clean$event)  

# Select predictors
X <- as.matrix(df_clean[, c("study_aric", "study_cardia", "study_dppos", "study_jhs", 
                      "race_nhwhi", "race_nhbla", "race_nhoth", "race_hisp",
                      "age", "female", "bmi", "hba1c", "sbp", "dbp", "hdlc", "ldlc",
                      "homa2b", "homa2ir", "height", "wc", "glucosef", "glucose2h", "tgl",
                      "serumcreatinine", "urinecreatinine", "egfr_ckdepi_2021", "apo_a", "apo_b", "uric_acid",
                      "vldlc", "hc", "triceps", "iliac", "medial", "ast", "alt", "urinealbumin",
                      "weight", "insulinf", "abdominal")])

### lasso ###
library(ncvreg) 

# Variables in your model
n <- nrow(X)   # Number of observations
p <- ncol(X)   # Number of predictors

cv_lasso <- cv.ncvsurv(X, y, penalty = "lasso",  # Choose among "MCP", "SCAD", "lasso"
                   gamma = 3.7,  # Relevant for SCAD and MCP
                   alpha = 1, 
                   lambda.min = ifelse(n > p, 0.001, 0.05), 
                   nlambda = 100,
                   eps = 1e-04, 
                   max.iter = 10000, 
                   convex = TRUE, 
                   dfmax = p, 
                   penalty.factor = rep(1, p),
                   warn = TRUE)

plot(cv_lasso)
best_lambda <- cv_lasso$lambda.min

final_lasso_model <- ncvsurv(X, y, penalty = "lasso", lambda = best_lambda)
print(coef(final_lasso_model))


### SCAD ###
# SCAD: more accurately estimate significant coefficients, thus potentially keeping more variables 
# in the model and often providing better estimation accuracy for non-zero coefficients.
cv_scad <- cv.ncvsurv(X, y, penalty = "SCAD",  # Choose among "MCP", "SCAD", "lasso"
                   gamma = 3.7,  # Relevant for SCAD and MCP
                   alpha = 1, 
                   lambda.min = ifelse(n > p, 0.001, 0.05), 
                   nlambda = 100,
                   eps = 1e-04, 
                   max.iter = 10000, 
                   convex = TRUE, 
                   dfmax = p, 
                   penalty.factor = rep(1, p),
                   warn = TRUE)


plot(cv_scad)
best_lambda <- cv_scad$lambda.min

final_scad_model <- ncvsurv(X, y, penalty = "SCAD", lambda = best_lambda)
print(coef(final_scad_model))

# non-zero: study, homa2ir, race, hba1c, uric_acid, wc, glucosef, age, glucose2h, weight, homa2b, insulinf, alt, female
# abdominal, egfr_ckdepi_2021, tgl, sbp, hdlc, apo_a
#---------------------------------------------------------------------------------------------------------------------------------
# cox PH model
cross_df <- analytic_sample %>% 
  group_by(study_id,study) %>% 
  dplyr::filter(age == min(age)) %>% 
  ungroup()

table(cross_df$event)

# check collinearity
# >0.7: homa2ir ~ homa2b, homa2ir ~ insulinf, wc ~ weight, insulinf ~ homa2b
# keep: insulinf, wc
cor(cross_df[, c("homa2ir", "hba1c", "uric_acid", "wc", "glucosef", "age", "glucose2h", "weight", 
                 "homa2b", "insulinf", "alt", "female")]) 

mard_cp <- coxph(Surv(time_to_event, mard) ~ study_aric + study_cardia + study_dppos + study_jhs + race
                 + hba1c + uric_acid + wc + glucosef + age + glucose2h + insulinf + alt + female, 
                 data = cross_df)

mod_cp <- coxph(Surv(time_to_event, mod) ~ study_aric + study_cardia + study_dppos + study_jhs + race
                + hba1c + uric_acid + wc + glucosef + age + glucose2h + insulinf + alt + female, 
                data = cross_df)

# get convergence warning from adding race as a covariate
# sidd == 1: Hisp, 9; NH Bla, 51; NH Oth, 0; NH Whi, 30; Other, 2
# merge "NH Other" category to "Other"
sidd_df <- cross_df %>%
  mutate(race = ifelse(race == "NH Other", "Other", race))

sidd_cp <- coxph(Surv(time_to_event, sidd) ~ study_aric + study_cardia + study_dppos + study_jhs + strata(race)
                 + hba1c + uric_acid + wc + glucosef + age + glucose2h + insulinf + alt + female, 
                 data = sidd_df)

sird_cp <- coxph(Surv(time_to_event, sird) ~ study_aric + study_cardia + study_dppos + study_jhs + race
                 + hba1c + uric_acid + wc + glucosef + age + glucose2h + insulinf + alt + female, 
                 data = cross_df)


coxph_results <- bind_rows(
  broom::tidy(mard_cp) %>% mutate(model = "MARD"),
  broom::tidy(mod_cp) %>% mutate(model = "MOD"),
  broom::tidy(sidd_cp) %>% mutate(model = "SIDD"),
  broom::tidy(sird_cp) %>% mutate(model = "SIRD")) 

# convert to HR
coxph_output <- coxph_results %>% 
  mutate(HR = exp(estimate),
         lci = exp(estimate - 1.96 * std.error),
         uci = exp(estimate + 1.96 * std.error)) %>% 
  mutate(coef_ci = paste0(round(HR, 2), " (", round(lci, 2), ", ", round(uci, 2), ")")) %>% 
  pivot_wider(names_from = model, values_from = coef_ci) %>% 
  dplyr::select(term, MARD, MOD, SIDD, SIRD)

df1 <- coxph_output %>% 
  dplyr::select(term, MARD) 
df2 <- coxph_output %>% 
  dplyr::select(term, MOD)
df3 <- coxph_output %>% 
  dplyr::select(term, SIDD)
df4 <- coxph_output %>% 
  dplyr::select(term, SIRD)

output <- na.omit(df1) %>% 
  left_join(na.omit(df2), by = "term") %>% 
  left_join(na.omit(df3), by = "term") %>% 
  left_join(na.omit(df4), by = "term")  %>% 
  write_csv("analysis/dspan02_imputed cox ph results.csv")

#--------------------------------------------------------------------
overall_cp_var8 <- coxph(Surv(time_to_event, event) ~ age + bmi + hba1c + homa2b + homa2ir + ldlc + sbp + egfr_ckdepi_2021, 
                      data = cross_df)

mard_cp_var8 <- coxph(Surv(time_to_event, mard) ~ age + bmi + hba1c + homa2b + homa2ir + ldlc + sbp + egfr_ckdepi_2021, 
                 data = cross_df)

mod_cp_var8 <- coxph(Surv(time_to_event, mod) ~ age + bmi + hba1c + homa2b + homa2ir + ldlc + sbp + egfr_ckdepi_2021, 
                data = cross_df)

sidd_cp_var8 <- coxph(Surv(time_to_event, sidd) ~ age + bmi + hba1c + homa2b + homa2ir + ldlc + sbp + egfr_ckdepi_2021, 
                 data = cross_df)

sird_cp_var8 <- coxph(Surv(time_to_event, sird) ~ age + bmi + hba1c + homa2b + homa2ir + ldlc + sbp + egfr_ckdepi_2021, 
                 data = cross_df)


coxph_results <- bind_rows(
  broom::tidy(overall_cp_var8) %>% mutate(model = "Overall"),
  broom::tidy(mard_cp_var8) %>% mutate(model = "MARD"),
  broom::tidy(mod_cp_var8) %>% mutate(model = "MOD"),
  broom::tidy(sidd_cp_var8) %>% mutate(model = "SIDD"),
  broom::tidy(sird_cp_var8) %>% mutate(model = "SIRD")) 

# convert to HR
coxph_output <- coxph_results %>% 
  mutate(HR = exp(estimate),
         lci = exp(estimate - 1.96 * std.error),
         uci = exp(estimate + 1.96 * std.error)) %>% 
  mutate(coef_ci = paste0(round(HR, 2), " (", round(lci, 2), ", ", round(uci, 2), ")")) %>% 
  pivot_wider(names_from = model, values_from = coef_ci) %>% 
  dplyr::select(term, Overall, MARD, MOD, SIDD, SIRD)

df0 <- coxph_output %>% 
  dplyr::select(term, Overall)
df1 <- coxph_output %>% 
  dplyr::select(term, MARD) 
df2 <- coxph_output %>% 
  dplyr::select(term, MOD)
df3 <- coxph_output %>% 
  dplyr::select(term, SIDD)
df4 <- coxph_output %>% 
  dplyr::select(term, SIRD)

output <- na.omit(df0) %>% 
  left_join(na.omit(df1), by = "term") %>% 
  left_join(na.omit(df2), by = "term") %>% 
  left_join(na.omit(df3), by = "term") %>% 
  left_join(na.omit(df4), by = "term")  %>% 
  write_csv("analysis/dspan02_imputed cox ph results with 8 variables.csv")
#------------------------------------------------------------------------------------------------------------------------
# TDCM analysis
# Start time: previous age; End time: current age at visit

tdcm_df <- analytic_sample %>%
  arrange(study, study_id, age) %>%
  group_by(study, study_id) %>%
  mutate(
    # Creating baseline variables by using the `lag()` function to get the previous age's value
    # Baseline values
    baseline_hba1c = if_else(row_number() == 1, hba1c, dplyr::lag(hba1c, 1, default = NA)),
    baseline_uric_acid = if_else(row_number() == 1, uric_acid, dplyr::lag(uric_acid, 1, default = NA)),
    baseline_wc = if_else(row_number() == 1, wc, dplyr::lag(wc, 1, default = NA)),
    baseline_glucosef = if_else(row_number() == 1, glucosef, dplyr::lag(glucosef, 1, default = NA)),
    baseline_glucose2h = if_else(row_number() == 1, glucose2h, dplyr::lag(glucose2h, 1, default = NA)),
    baseline_insulinf = if_else(row_number() == 1, insulinf, dplyr::lag(insulinf, 1, default = NA)),
    baseline_alt = if_else(row_number() == 1, alt, dplyr::lag(alt, 1, default = NA)),
    
    # Change calculations
    change_hba1c = hba1c - baseline_hba1c,
    change_uric_acid = uric_acid - baseline_uric_acid,
    change_wc = wc - baseline_wc,
    change_glucosef = glucosef - baseline_glucosef,
    change_glucose2h = glucose2h - baseline_glucose2h,
    change_insulinf = insulinf - baseline_insulinf,
    change_alt = alt - baseline_alt,
    
    # for 8 variables
    baseline_bmi = if_else(row_number() == 1, bmi, dplyr::lag(bmi, 1, default = NA)),
    baseline_egfr = if_else(row_number() == 1, egfr_ckdepi_2021, dplyr::lag(egfr_ckdepi_2021, 1, default = NA)),
    baseline_ldlc = if_else(row_number() == 1, ldlc, dplyr::lag(ldlc, 1, default = NA)),
    baseline_homa2b = if_else(row_number() == 1, homa2b, dplyr::lag(homa2b, 1, default = NA)),
    baseline_homa2ir = if_else(row_number() == 1, homa2ir, dplyr::lag(homa2ir, 1, default = NA)),
    baseline_sbp = if_else(row_number() == 1, sbp, dplyr::lag(sbp, 1, default = NA)),
    
    change_bmi = bmi - baseline_bmi,
    change_egfr = egfr - baseline_egfr,
    change_ldlc = ldlc - baseline_ldlc,
    change_homa2b = homa2b - baseline_homa2b,
    change_homa2ir = homa2ir - baseline_homa2ir,
    change_sbp = sbp - baseline_sbp,
    
  ) %>%
  mutate(
    tstart = case_when(row_number() == 1 ~ age, 
                       TRUE ~ dplyr::lag(age, n = 1)), 
    tstop = age
  ) %>%
  ungroup() %>% 
  # dplyr::filter(tstart < tstop)
  dplyr::filter((tstart < tstop) & (tstop <= censored_age))

# check
test <- tdcm_df %>% select(study, study_id, age, alt,baseline_alt,change_alt, tstart,tstop)

# "age" leads to convergence issues, exclude it
# exclude study_aric and study_cardia because obs = 0 in these studies (when event == 1)
tdcm_mod <- coxph(Surv(tstart, tstop, event) ~ study_jhs + study_dppos + race
                  + hba1c + uric_acid + wc + glucosef + insulinf + alt + female, data = tdcm_df)

mard_tdcm <- coxph(Surv(tstart, tstop, mard) ~ study_jhs + study_dppos + race
                   + hba1c + uric_acid + wc + glucosef + insulinf + alt + female, data = tdcm_df)
mod_tdcm <- coxph(Surv(tstart, tstop, mod) ~ study_jhs + study_dppos + race
                  + hba1c + uric_acid + wc + glucosef + insulinf + alt + female, data = tdcm_df)

# sidd == 1: Hisp, 9; NH Bla, 51; NH Oth, 0; NH Whi, 30; Other, 2
# merge "NH Other" category to "Other"
sidd_df <- tdcm_df %>%
  mutate(race = ifelse(race == "NH Other", "Other", race))

sidd_tdcm <- coxph(Surv(tstart, tstop, sidd) ~ study_jhs + study_dppos + race
                   + hba1c + uric_acid + wc + glucosef + insulinf + alt + female, data = sidd_df)
sird_tdcm <- coxph(Surv(tstart, tstop, sird) ~ study_jhs + study_dppos + race
                   + hba1c + uric_acid + wc + glucosef + insulinf + alt + female, data = tdcm_df)

tdcm_results <- bind_rows(
  broom::tidy(mard_tdcm) %>% mutate(model = "MARD"),
  broom::tidy(mod_tdcm) %>% mutate(model = "MOD"),
  broom::tidy(sidd_tdcm) %>% mutate(model = "SIDD"),
  broom::tidy(sird_tdcm) %>% mutate(model = "SIRD")) 

# covert to Hazard Ratio
tdcm_output <- tdcm_results %>% 
  mutate(HR = exp(estimate),
         lci = exp(estimate - 1.96 * std.error),
         uci = exp(estimate + 1.96 * std.error)) %>% 
  mutate(coef_ci = paste0(round(HR, 2), " (", round(lci, 2), ", ", round(uci, 2), ")")) %>% 
  pivot_wider(names_from = model, values_from = coef_ci) %>% 
  dplyr::select(term, MARD, MOD, SIDD, SIRD) 

df1 <- tdcm_output %>% 
  dplyr::select(term, MARD) 
df2 <- tdcm_output %>% 
  dplyr::select(term, MOD)
df3 <- tdcm_output %>% 
  dplyr::select(term, SIDD)
df4 <- tdcm_output %>% 
  dplyr::select(term, SIRD)

output <- na.omit(df1) %>% 
  left_join(na.omit(df2), by = "term") %>% 
  left_join(na.omit(df3), by = "term") %>% 
  left_join(na.omit(df4), by = "term") %>% 
  write_csv("analysis/dspan02_imputed tdcm results.csv")

#----------------------------------------------------------------
overall_tdcm_var7 <- coxph(Surv(tstart, tstop, event) ~ bmi + hba1c + homa2b + homa2ir + ldlc + sbp + egfr_ckdepi_2021, data = tdcm_df)

mard_tdcm_var7 <- coxph(Surv(tstart, tstop, mard) ~ bmi + hba1c + homa2b + homa2ir + ldlc + sbp + egfr_ckdepi_2021, data = tdcm_df)
mod_tdcm_var7 <- coxph(Surv(tstart, tstop, mod) ~ bmi + hba1c + homa2b + homa2ir + ldlc + sbp + egfr_ckdepi_2021, data = tdcm_df)
sidd_tdcm_var7 <- coxph(Surv(tstart, tstop, sidd) ~ bmi + hba1c + homa2b + homa2ir + ldlc + sbp + egfr_ckdepi_2021, data = tdcm_df)
sird_tdcm_var7 <- coxph(Surv(tstart, tstop, sird) ~ bmi + hba1c + homa2b + homa2ir + ldlc + sbp + egfr_ckdepi_2021, data = tdcm_df)

tdcm_results <- bind_rows(
  broom::tidy(overall_tdcm_var7) %>% mutate(model = "Overall"),
  broom::tidy(mard_tdcm_var7) %>% mutate(model = "MARD"),
  broom::tidy(mod_tdcm_var7) %>% mutate(model = "MOD"),
  broom::tidy(sidd_tdcm_var7) %>% mutate(model = "SIDD"),
  broom::tidy(sird_tdcm_var7) %>% mutate(model = "SIRD")) 

# covert to Hazard Ratio
tdcm_output <- tdcm_results %>% 
  mutate(HR = exp(estimate),
         lci = exp(estimate - 1.96 * std.error),
         uci = exp(estimate + 1.96 * std.error)) %>% 
  mutate(coef_ci = paste0(round(HR, 2), " (", round(lci, 2), ", ", round(uci, 2), ")")) %>% 
  pivot_wider(names_from = model, values_from = coef_ci) %>% 
  dplyr::select(term, Overall, MARD, MOD, SIDD, SIRD) 

df0 <- tdcm_output %>% 
  dplyr::select(term, Overall) 
df1 <- tdcm_output %>% 
  dplyr::select(term, MARD) 
df2 <- tdcm_output %>% 
  dplyr::select(term, MOD)
df3 <- tdcm_output %>% 
  dplyr::select(term, SIDD)
df4 <- tdcm_output %>% 
  dplyr::select(term, SIRD)

output <- na.omit(df0) %>% 
  left_join(na.omit(df1), by = "term") %>%
  left_join(na.omit(df2), by = "term") %>% 
  left_join(na.omit(df3), by = "term") %>% 
  left_join(na.omit(df4), by = "term") %>% 
  write_csv("analysis/dspan02_imputed tdcm results with 7 variables.csv")
#------------------------------------------------------------------------------------------------------------------------
# TDCM change analysis - baseline + change
 
tdcm_change <- coxph(Surv(tstart, tstop, event) ~ study_jhs + study_dppos + race
                     + baseline_hba1c + baseline_uric_acid + baseline_wc + baseline_glucosef 
                     + baseline_insulinf + baseline_alt + female
                     + change_hba1c + change_uric_acid + change_wc + change_glucosef 
                     + change_insulinf + change_alt, data = tdcm_df)

mard_change <- coxph(Surv(tstart, tstop, mard) ~ study_jhs + study_dppos + race
                     + baseline_hba1c + baseline_uric_acid + baseline_wc + baseline_glucosef 
                     + baseline_insulinf + baseline_alt + female
                     + change_hba1c + change_uric_acid + change_wc + change_glucosef 
                     + change_insulinf + change_alt, data = tdcm_df)

mod_change <- coxph(Surv(tstart, tstop, mod) ~ study_jhs + study_dppos + race
                    + baseline_hba1c + baseline_uric_acid + baseline_wc + baseline_glucosef 
                    + baseline_insulinf + baseline_alt + female
                    + change_hba1c + change_uric_acid + change_wc + change_glucosef 
                    + change_insulinf + change_alt, data = tdcm_df)

sidd_df <- tdcm_df %>%
  mutate(race = ifelse(race == "NH Other", "Other", race))

sidd_change <- coxph(Surv(tstart, tstop, sidd) ~ study_jhs + study_dppos + race
                     + baseline_hba1c + baseline_uric_acid + baseline_wc + baseline_glucosef 
                     + baseline_insulinf + baseline_alt + female
                     + change_hba1c + change_uric_acid + change_wc + change_glucosef 
                     + change_insulinf + change_alt, data = sidd_df)

sird_change <- coxph(Surv(tstart, tstop, sird) ~ study_jhs + study_dppos + race
                     + baseline_hba1c + baseline_uric_acid + baseline_wc + baseline_glucosef 
                     + baseline_insulinf + baseline_alt + female
                     + change_hba1c + change_uric_acid + change_wc + change_glucosef 
                     + change_insulinf + change_alt, data = tdcm_df)

tdcm_change_results <- bind_rows(
  broom::tidy(mard_change) %>% mutate(model = "MARD"),
  broom::tidy(mod_change) %>% mutate(model = "MOD"),
  broom::tidy(sidd_change) %>% mutate(model = "SIDD"),
  broom::tidy(sird_change) %>% mutate(model = "SIRD")) 

tdcm_change_output <- tdcm_change_results %>% 
  mutate(HR = exp(estimate),
         lci = exp(estimate - 1.96 * std.error),
         uci = exp(estimate + 1.96 * std.error)) %>% 
  mutate(coef_ci = paste0(round(HR, 2), " (", round(lci, 2), ", ", round(uci, 2), ")")) %>% 
  pivot_wider(names_from = model, values_from = coef_ci) %>% 
  dplyr::select(term, MARD, MOD, SIDD, SIRD) 

df1 <- tdcm_change_output %>% 
  dplyr::select(term, MARD) 
df2 <- tdcm_change_output %>% 
  dplyr::select(term, MOD)
df3 <- tdcm_change_output %>% 
  dplyr::select(term, SIDD)
df4 <- tdcm_change_output %>% 
  dplyr::select(term, SIRD)

output <- na.omit(df1) %>% 
  left_join(na.omit(df2), by = "term") %>% 
  left_join(na.omit(df3), by = "term") %>% 
  left_join(na.omit(df4), by = "term") %>% 
  write_csv("analysis/dspan02_imputed tdcm change results.csv")

#-------------------------------------------------------------------
overall_change_var7 <- coxph(Surv(tstart, tstop, event) ~ baseline_hba1c + baseline_homa2b + baseline_homa2ir + baseline_ldlc 
                     + baseline_bmi + baseline_sbp + baseline_egfr
                     + change_hba1c + change_homa2b + change_homa2ir + change_ldlc 
                     + change_bmi + change_sbp + change_egfr, data = tdcm_df)

mard_change_var7 <- coxph(Surv(tstart, tstop, mard) ~ baseline_hba1c + baseline_homa2b + baseline_homa2ir + baseline_ldlc 
                     + baseline_bmi + baseline_sbp + baseline_egfr
                     + change_hba1c + change_homa2b + change_homa2ir + change_ldlc 
                     + change_bmi + change_sbp + change_egfr, data = tdcm_df)

mod_change_var7 <- coxph(Surv(tstart, tstop, mod) ~ baseline_hba1c + baseline_homa2b + baseline_homa2ir + baseline_ldlc 
                    + baseline_bmi + baseline_sbp + baseline_egfr
                    + change_hba1c + change_homa2b + change_homa2ir + change_ldlc 
                    + change_bmi + change_sbp + change_egfr, data = tdcm_df)

sidd_change_var7 <- coxph(Surv(tstart, tstop, sidd) ~ baseline_hba1c + baseline_homa2b + baseline_homa2ir + baseline_ldlc 
                     + baseline_bmi + baseline_sbp + baseline_egfr
                     + change_hba1c + change_homa2b + change_homa2ir + change_ldlc 
                     + change_bmi + change_sbp + change_egfr, data = tdcm_df)

sird_change_var7 <- coxph(Surv(tstart, tstop, sird) ~ baseline_hba1c + baseline_homa2b + baseline_homa2ir + baseline_ldlc 
                     + baseline_bmi + baseline_sbp + baseline_egfr
                     + change_hba1c + change_homa2b + change_homa2ir + change_ldlc 
                     + change_bmi + change_sbp + change_egfr, data = tdcm_df)

tdcm_change_results <- bind_rows(
  broom::tidy(overall_change_var7) %>% mutate(model = "Overall"),
  broom::tidy(mard_change_var7) %>% mutate(model = "MARD"),
  broom::tidy(mod_change_var7) %>% mutate(model = "MOD"),
  broom::tidy(sidd_change_var7) %>% mutate(model = "SIDD"),
  broom::tidy(sird_change_var7) %>% mutate(model = "SIRD")) 

tdcm_change_output <- tdcm_change_results %>% 
  mutate(HR = exp(estimate),
         lci = exp(estimate - 1.96 * std.error),
         uci = exp(estimate + 1.96 * std.error)) %>% 
  mutate(coef_ci = paste0(round(HR, 2), " (", round(lci, 2), ", ", round(uci, 2), ")")) %>% 
  pivot_wider(names_from = model, values_from = coef_ci) %>% 
  dplyr::select(term, Overall, MARD, MOD, SIDD, SIRD) 

df0 <- tdcm_change_output %>% 
  dplyr::select(term, Overall) 
df1 <- tdcm_change_output %>% 
  dplyr::select(term, MARD) 
df2 <- tdcm_change_output %>% 
  dplyr::select(term, MOD)
df3 <- tdcm_change_output %>% 
  dplyr::select(term, SIDD)
df4 <- tdcm_change_output %>% 
  dplyr::select(term, SIRD)

output <- na.omit(df0) %>% 
  left_join(na.omit(df1), by = "term") %>% 
  left_join(na.omit(df2), by = "term") %>% 
  left_join(na.omit(df3), by = "term") %>% 
  left_join(na.omit(df4), by = "term") %>% 
  write_csv("analysis/dspan02_imputed tdcm change results with 7 variables.csv")
#------------------------------------------------------------------------------------------------------------------------
# causal survival forest
# chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://arxiv.org/pdf/2312.02482

cluster_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan02_cluster df.RDS"))

# outcome
Y = cluster_df$time_to_event 
# subgroup indicator
W_mard = cluster_df$mard
W_mod = cluster_df$mod
W_sidd = cluster_df$sidd
W_sird = cluster_df$sird
# non-censoring indicator
D = cluster_df$event 
# covariates
X = cbind(
  study_aric = cluster_df$study_aric,
  study_cardia = cluster_df$study_cardia,
  study_dppos = cluster_df$study_dppos,
  study_jhs = cluster_df$study_jhs,
  race_nhwhi = cluster_df$race_nhwhi,
  race_nhbla = cluster_df$race_nhbla,
  race_nhoth = cluster_df$race_nhoth,
  race_hisp = cluster_df$race_hisp,
  hba1c = cluster_df$hba1c,
  uric_acid = cluster_df$uric_acid,
  wc = cluster_df$wc,
  glucosef = cluster_df$glucosef,
  age = cluster_df$age,
  glucose2h = cluster_df$glucose2h,
  insulinf = cluster_df$insulinf,
  alt = cluster_df$alt,
  female = cluster_df$female
  
)

failure_times <- seq(0, 15, by=0.5)

# set the truncation time to 20 years, at which point most subsequent observations are censored.
hist(Y[D==1],main="HistogramofY", xlab="")
hist(Y[D==0], col=adjustcolor("red", 0.5), add=TRUE);
legend("topright", c("Event", "Censored"), fill=c("gray", adjustcolor("red", 0.5)), lwd=4)

abline(v=720,lty=2)

library(grf)
csf_mard <- causal_survival_forest(X, Y, W_mard, D, W.hat = 0.5, target = "RMST", horizon = 15, failure.times = failure_times)
csf_mod <- causal_survival_forest(X, Y, W_mod, D, W.hat = 0.5, target = "RMST", horizon = 15, failure.times = failure_times)
csf_sidd <- causal_survival_forest(X, Y, W_sidd, D, W.hat = 0.5, target = "RMST", horizon = 15, failure.times = failure_times)
csf_sird <- causal_survival_forest(X, Y, W_sird, D, W.hat = 0.5, target = "RMST", horizon = 15, failure.times = failure_times)

blp_mard <- best_linear_projection(csf_mard, X)
blp_mod <- best_linear_projection(csf_mod, X)
blp_sidd <- best_linear_projection(csf_sidd, X)
blp_sird <- best_linear_projection(csf_sird, X)

csf_results <- bind_rows(as.data.frame.matrix(blp_mard) %>% 
                          mutate(model = "MARD"),
                        as.data.frame.matrix(blp_mod) %>% 
                          mutate(model = "MOD"),
                        as.data.frame.matrix(blp_sidd) %>% 
                          mutate(model = "SIDD"),
                        as.data.frame.matrix(blp_sird) %>% 
                          mutate(model = "SIRD")) %>% 
  rename(estimate = Estimate,
         std_error = `Std. Error`,
         t_value = `t value`,
         p_value = `Pr(>|t|)`) %>% 
  rownames_to_column(var = "term") %>% 
  mutate(term = gsub("\\.{3}\\d+$", "", term))

csf_output <- csf_results %>% 
  mutate(HR = exp(estimate),
         lci = exp(estimate - 1.96 * std_error),
         uci = exp(estimate + 1.96 * std_error)) %>% 
  mutate(coef_ci = paste0(round(HR, 2), " (", round(lci, 2), ", ", round(uci, 2), ")")) %>% 
  pivot_wider(names_from = model, values_from = coef_ci) %>% 
  dplyr::select(term, MARD, MOD, SIDD, SIRD) 

df1 <- csf_output %>% 
  dplyr::select(term, MARD) 
df2 <- csf_output %>% 
  dplyr::select(term, MOD)
df3 <- csf_output %>% 
  dplyr::select(term, SIDD)
df4 <- csf_output %>% 
  dplyr::select(term, SIRD)

output <- na.omit(df1) %>% 
  left_join(na.omit(df2), by = "term") %>% 
  left_join(na.omit(df3), by = "term") %>% 
  left_join(na.omit(df4), by = "term") %>% 
  write_csv("analysis/dspan02_imputed causal survival forest results.csv")

average_treatment_effect(csf)
# Retrieve out-of-bag CATE estimates 
tau.hat= predict(csf)$predictions 
summary(tau.hat)

#-------------------------------------------------------------------------------
# covariates 
X = cbind(
  hba1c = cluster_df$hba1c,
  bmi = cluster_df$bmi,
  ldlc = cluster_df$ldlc,
  homa2b = cluster_df$homa2b,
  homa2ir = cluster_df$homa2ir,
  age = cluster_df$age,
  sbp = cluster_df$sbp,
  egfr_ckdepi_2021 = cluster_df$egfr_ckdepi_2021
  
)

failure_times <- seq(0, 15, by=0.5)

library(grf)
csf_mard_var8 <- causal_survival_forest(X, Y, W_mard, D, W.hat = 0.5, target = "RMST", horizon = 15, failure.times = failure_times)
csf_mod_var8 <- causal_survival_forest(X, Y, W_mod, D, W.hat = 0.5, target = "RMST", horizon = 15, failure.times = failure_times)
csf_sidd_var8 <- causal_survival_forest(X, Y, W_sidd, D, W.hat = 0.5, target = "RMST", horizon = 15, failure.times = failure_times)
csf_sird_var8 <- causal_survival_forest(X, Y, W_sird, D, W.hat = 0.5, target = "RMST", horizon = 15, failure.times = failure_times)

blp_mard_var8 <- best_linear_projection(csf_mard_var8, X)
blp_mod_var8 <- best_linear_projection(csf_mod_var8, X)
blp_sidd_var8 <- best_linear_projection(csf_sidd_var8, X)
blp_sird_var8 <- best_linear_projection(csf_sird_var8, X)

csf_results <- bind_rows(as.data.frame.matrix(blp_mard_var8) %>% 
                           mutate(model = "MARD"),
                         as.data.frame.matrix(blp_mod_var8) %>% 
                           mutate(model = "MOD"),
                         as.data.frame.matrix(blp_sidd_var8) %>% 
                           mutate(model = "SIDD"),
                         as.data.frame.matrix(blp_sird_var8) %>% 
                           mutate(model = "SIRD")) %>% 
  rename(estimate = Estimate,
         std_error = `Std. Error`,
         t_value = `t value`,
         p_value = `Pr(>|t|)`) %>% 
  rownames_to_column(var = "term") %>% 
  mutate(term = gsub("\\.{3}\\d+$", "", term))

csf_output <- csf_results %>% 
  mutate(HR = exp(estimate),
         lci = exp(estimate - 1.96 * std_error),
         uci = exp(estimate + 1.96 * std_error)) %>% 
  mutate(coef_ci = paste0(round(HR, 2), " (", round(lci, 2), ", ", round(uci, 2), ")")) %>% 
  pivot_wider(names_from = model, values_from = coef_ci) %>% 
  dplyr::select(term, MARD, MOD, SIDD, SIRD) 

df1 <- csf_output %>% 
  dplyr::select(term, MARD) 
df2 <- csf_output %>% 
  dplyr::select(term, MOD)
df3 <- csf_output %>% 
  dplyr::select(term, SIDD)
df4 <- csf_output %>% 
  dplyr::select(term, SIRD)

output <- na.omit(df1) %>% 
  left_join(na.omit(df2), by = "term") %>% 
  left_join(na.omit(df3), by = "term") %>% 
  left_join(na.omit(df4), by = "term") %>% 
  write_csv("dspan02_imputed causal survival forest results with 8 variables.csv")
