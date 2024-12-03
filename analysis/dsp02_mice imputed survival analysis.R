rm(list = ls());gc();source(".Rprofile")

library(survival)
library(survminer)
library(ggsurvfit)
library(broom)
library(tidyr)

# follow-up time >15y & dm == 0 --- event == 0
mice_df <- readRDS("analysis/mi_dfs.RDS") %>% 
  group_by(study, study_id) %>% 
  mutate(min_age = min(age),
         max_age = max(age)) %>%
  # dplyr::filter(age <= min_age + 15) %>%
  select(-max_age) %>% 
  ungroup() %>% 
  mutate(dm = case_when(!is.na(dmagediag) | hba1c >= 6.5 | glucosef >= 126 | glucose2h >= 200 ~ 1,
                        TRUE ~ 0),
         event = case_when(age <= min_age + 15 ~ dm,
                           TRUE ~ 0)) %>% 
  mutate(study_aric = case_when(study == "aric" ~ 1,
                                TRUE ~ 0),
         study_cardia = case_when(study == "cardia" ~ 1,
                                TRUE ~ 0),
         study_dppos = case_when(study == "dppos" ~ 1,
                                TRUE ~ 0),
         study_jhs = case_when(study == "jhs" ~ 1,
                                TRUE ~ 0),
         study_mesa = case_when(study == "mesa" ~ 1,
                                TRUE ~ 0),
         race_nhwhi = case_when(race_eth == "NH White" ~ 1,
                                TRUE ~ 0),
         race_nhbla = case_when(race_eth == "NH Black" ~ 1,
                                TRUE ~ 0),
         race_nhoth = case_when(race_eth == "NH Other" ~ 1,
                                TRUE ~ 0),
         race_hisp = case_when(race_eth == "Hispanic" ~ 1,
                                TRUE ~ 0))

# define time to event
time_df <- mice_df %>% 
  group_by(study, study_id) %>% 
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
# cross-section datatse, cox PH model
# overall, include all variables, 1 obs for each person
cross_df <- time_df %>% 
  dplyr::filter(age == min(age))

cox_mod <- coxph(Surv(time_to_event, event) ~ study_aric + study_cardia + study_dppos + study_jhs + study_mesa + age + female 
                                           + race_nhwhi + race_nhbla + race_nhoth + race_hisp + bmi + hba1c + sbp + dbp + hdlc + ldlc 
                                           + homa2b + homa2ir + height + wc + insulinf + glucosef + glucose2h + tgl
                                           + serumcreatinine + urinecreatinine + egfr + apo_a + apo_b + uric_acid 
                                           + vldlc + hc + triceps + iliac + abdominal + medial + ast + alt + urinealbumin
                                           + weight, data = cross_df)
summary(cox_mod)

# p>.05 covariates: insulinf, abdominal

#------------------------------------------------------------------------------------------------------------------------
final_dataset_temp = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/final_dataset_temp.RDS"))
clusters = read_csv(paste0(path_diabetes_subphenotypes_adults_folder,"/working/processed/dec_an02_clean_kmeans_5var_mi_knn_cluster.csv")) %>% 
  dplyr::select(-one_of("...1")) %>% 
  left_join(final_dataset_temp %>% 
              dplyr::select(study_id,original_study_id),
            by=c("study_id")) %>% 
  rename(cluster_study_id = study_id)

# assign cluster
cluster_df <- time_df %>% 
  left_join(clusters %>% 
              dplyr::select(cluster_study_id,original_study_id,cluster,study),
            by=c("study"="study","study_id" = "original_study_id", "cluster_study_id")) %>% 
  ### exclude dm == 1 with no cluster data
  dplyr::filter(!(dm == 1 & is.na(cluster)))


# cause-specific HR
analytic_sample <- cluster_df %>% 
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

#--------------------------------------------------------------------------------------------------------------------------
# variable selection - LASSO, SCAD
# chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cran.r-project.org/web/packages/ncvreg/ncvreg.pdf
# require all variables to be numeric

df <- analytic_sample %>%  
  dplyr::filter(time_to_event > 0)
  
y <- Surv(df$time_to_event, df$mard)  

# Select predictors
X <- as.matrix(df[, c("study_aric", "study_cardia", "study_dppos", "study_jhs", "study_mesa",
                      "race_nhwhi", "race_nhbla", "race_nhoth", "race_hisp",
                      "age", "female", "bmi", "hba1c", "sbp", "dbp", "hdlc", "ldlc",
                      "homa2b", "homa2ir", "height", "wc", "glucosef", "glucose2h", "tgl",
                      "serumcreatinine", "urinecreatinine", "egfr", "apo_a", "apo_b", "uric_acid",
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


### lasso + Ridge regression ###

# glmnet is designed to implement LASSO (Least Absolute Shrinkage and Selection Operator) and Ridge regression methods
library(glmnet)
# Assuming y is your response vector and X is the predictor matrix
cv_model <- cv.glmnet(X, y, family = "cox", alpha = 1)  # LASSO regression

# Plot the cross-validation results
plot(cv_model)

# Best lambda (smallest error)
best_lambda <- cv_model$lambda.min
print(best_lambda)

# Fit model on selected lambda
final_model <- glmnet(X, y, family = "cox", alpha = 1, lambda = best_lambda)
print(coef(final_model))


# influential vars: 
# >0.05: female, study, age, serumcreatinine
# >0.005: race_eth, hba1c, hc, glucosef, weight, dbp, ast, sbp
# >0.001: glucose2h,medial,ldlc,wc,iliac,alt,homa2b,homa2ir

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

# non-zero: study, serumcreatinine, race_eth, female, uric_acid, hba1c, bmi, hc, age, homa2ir, medial
# weight, dbp, wc, egfr, iliac, ldlc, triceps, ast, insulinf, alt
# hdlc, abdominal, sbp, glucose2h, tgl, apo_b, apo_a, urinealbumin
#---------------------------------------------------------------------------------------------------------------------------------
# cox PH model
mard_cp <- coxph(Surv(time_to_event, mard) ~ study_aric + study_cardia + study_dppos + study_jhs + study_mesa 
                 + serumcreatinine + race_nhwhi + race_nhbla + race_nhoth + race_hisp 
                 + female + uric_acid + hba1c + bmi + hc + age
                 + homa2ir + medial, data = analytic_sample)

mod_cp <- coxph(Surv(time_to_event, mod) ~ study_aric + study_cardia + study_dppos + study_jhs + study_mesa
                + serumcreatinine + race_nhwhi + race_nhbla + race_nhoth + race_hisp
                + female + uric_acid + hba1c + bmi + hc + age
                + homa2ir + medial, data = analytic_sample)

sidd_cp <- coxph(Surv(time_to_event, sidd) ~ study_aric + study_cardia + study_dppos + study_jhs + study_mesa
                 + serumcreatinine + race_nhwhi + race_nhbla + race_nhoth + race_hisp
                 + female + uric_acid + hba1c + bmi + hc + age
                 + homa2ir + medial, data = analytic_sample)

sird_cp <- coxph(Surv(time_to_event, sird) ~ study_aric + study_cardia + study_dppos + study_jhs + study_mesa
                 + serumcreatinine + race_nhwhi + race_nhbla + race_nhoth + race_hisp
                 + female + uric_acid + hba1c + bmi + hc + age
                 + homa2ir + medial, data = analytic_sample)


coxph_results <- bind_rows(
  broom::tidy(mard_cp) %>% mutate(model = "mard_cp"),
  broom::tidy(mod_cp) %>% mutate(model = "mod_cp"),
  broom::tidy(sidd_cp) %>% mutate(model = "sidd_cp"),
  broom::tidy(sird_cp) %>% mutate(model = "sird_cp")) 

# convert to HR
coxph_output <- coxph_results %>% 
  mutate(HR = exp(estimate),
         lci = exp(estimate - 1.96 * std.error),
         uci = exp(estimate + 1.96 * std.error),
         model = case_when(model == "mard_cp" ~ "MARD",
                           model == "mod_cp" ~ "MOD",
                           model == "sidd_cp" ~ "SIDD",
                           model == "sird_cp" ~ "SIRD",
                           TRUE ~ NA_character_)) %>% 
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
  write_csv("analysis/dsp02_imputed cox ph results.csv")
#------------------------------------------------------------------------------------------------------------------------
# TDCM analysis
# Start time: previous age; End time: current age at visit

tdcm_df <- analytic_sample %>%
  arrange(study, study_id, age) %>%
  group_by(study, study_id) %>%
  mutate(
    # Creating baseline variables by using the `lag()` function to get the previous age's value
    baseline_serumcreatinine = case_when(
      row_number() == 1 ~ serumcreatinine,
      TRUE ~ dplyr::lag(serumcreatinine, 1, default = NA)
    ),
    baseline_uric_acid = case_when(
      row_number() == 1 ~ uric_acid,
      TRUE ~ dplyr::lag(uric_acid, 1, default = NA)
    ),
    baseline_hba1c = case_when(
      row_number() == 1 ~ hba1c,
      TRUE ~ dplyr::lag(hba1c, 1, default = NA)
    ),
    baseline_bmi = case_when(
      row_number() == 1 ~ bmi,
      TRUE ~ dplyr::lag(bmi, 1, default = NA)
    ),
    baseline_homa2ir = case_when(
      row_number() == 1 ~ homa2ir,
      TRUE ~ dplyr::lag(homa2ir, 1, default = NA)
    ),
    baseline_medial = case_when(
      row_number() == 1 ~ medial,
      TRUE ~ dplyr::lag(medial, 1, default = NA)
    ),
    
    # Calculating change by subtracting the baseline value from the current value
    change_serumcreatinine = serumcreatinine - baseline_serumcreatinine,
    change_uric_acid = uric_acid - baseline_uric_acid,
    change_hba1c = hba1c - baseline_hba1c,
    change_bmi = bmi - baseline_bmi,
    change_homa2ir = homa2ir - baseline_homa2ir,
    change_medial = medial - baseline_medial
  ) %>%
  mutate(
    tstart = case_when(row_number() == 1 ~ age, 
                       TRUE ~ dplyr::lag(age, n = 1)), 
    tstop = age 
  ) %>%
  ungroup() %>% 
  dplyr::filter(tstart < tstop)

# check
test <- tdcm_df %>% 
  select(study, study_id, age, tstart,tstop)

test <- tdcm_df %>% select(study, study_id, age, serumcreatinine,baseline_serumcreatinine,change_serumcreatinine, tstart,tstop)

# collinearity
cor(tdcm_df[, c("serumcreatinine", "female", "uric_acid", "hba1c", "bmi", "hc", "age", "homa2ir", "medial")]) 


# exclude correlation coefficient > 0.7 -- highly correlated
# hba1c ~ glucosef, glucose 2h ~ glucosef, sbp ~ dbp, hc ~ bmi

tdcm_mod <- coxph(Surv(tstart, tstop, event) ~ study + serumcreatinine + race_eth + female + uric_acid + hba1c + bmi + age
                  + homa2ir + medial, data = tdcm_df)

mard_tdcm <- coxph(Surv(tstart, tstop, mard) ~ study + serumcreatinine + race_eth + female + uric_acid + hba1c + bmi + age
            + homa2ir + medial, data = tdcm_df)
mod_tdcm <- coxph(Surv(tstart, tstop, mod) ~ study + serumcreatinine + race_eth + female + uric_acid + hba1c + bmi + age
            + homa2ir + medial, data = tdcm_df)
sidd_tdcm <- coxph(Surv(tstart, tstop, sidd) ~ study + serumcreatinine + race_eth + female + uric_acid + hba1c + bmi + age
            + homa2ir + medial, data = tdcm_df)
sird_tdcm <- coxph(Surv(tstart, tstop, sird) ~ study + serumcreatinine + race_eth + female + uric_acid + hba1c + bmi + age
            + homa2ir + medial, data = tdcm_df)

tdcm_results <- bind_rows(
  broom::tidy(mard_tdcm) %>% mutate(model = "mard_tdcm"),
  broom::tidy(mod_tdcm) %>% mutate(model = "mod_tdcm"),
  broom::tidy(sidd_tdcm) %>% mutate(model = "sidd_tdcm"),
  broom::tidy(sird_tdcm) %>% mutate(model = "sird_tdcm")) 

# covert to Hazard Ratio
tdcm_output <- tdcm_results %>% 
  mutate(HR = exp(estimate),
         lci = exp(estimate - 1.96 * std.error),
         uci = exp(estimate + 1.96 * std.error),
         model = case_when(model == "mard_tdcm" ~ "MARD",
                           model == "mod_tdcm" ~ "MOD",
                           model == "sidd_tdcm" ~ "SIDD",
                           model == "sird_tdcm" ~ "SIRD",
                           TRUE ~ NA_character_)) %>% 
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
  write_csv("analysis/dsp02_imputed tdcm results.csv")
#------------------------------------------------------------------------------------------------------------------------
# TDCM change analysis - baseline + change
 
tdcm_change <- coxph(Surv(tstart, tstop, event) ~ study + age + race_eth + female 
                  + baseline_serumcreatinine + change_serumcreatinine + baseline_uric_acid + change_uric_acid + baseline_hba1c + change_hba1c 
                  + baseline_bmi + change_bmi + baseline_homa2ir + change_homa2ir + baseline_medial + change_medial, data = tdcm_df)

mard_change <- coxph(Surv(tstart, tstop, mard) ~ study + age + race_eth + female 
            + baseline_serumcreatinine + change_serumcreatinine + baseline_uric_acid + change_uric_acid + baseline_hba1c + change_hba1c 
            + baseline_bmi + change_bmi + baseline_homa2ir + change_homa2ir + baseline_medial + change_medial, data = tdcm_df)
mod_change <- coxph(Surv(tstart, tstop, mod) ~ study + age + race_eth + female 
            + baseline_serumcreatinine + change_serumcreatinine + baseline_uric_acid + change_uric_acid + baseline_hba1c + change_hba1c 
            + baseline_bmi + change_bmi + baseline_homa2ir + change_homa2ir + baseline_medial + change_medial, data = tdcm_df)
sidd_change <- coxph(Surv(tstart, tstop, sidd) ~ study + age + race_eth + female 
            + baseline_serumcreatinine + change_serumcreatinine + baseline_uric_acid + change_uric_acid + baseline_hba1c + change_hba1c 
            + baseline_bmi + change_bmi + baseline_homa2ir + change_homa2ir + baseline_medial + change_medial, data = tdcm_df)
sird_change <- coxph(Surv(tstart, tstop, sird) ~ study + age + race_eth + female 
            + baseline_serumcreatinine + change_serumcreatinine + baseline_uric_acid + change_uric_acid + baseline_hba1c + change_hba1c 
            + baseline_bmi + change_bmi + baseline_homa2ir + change_homa2ir + baseline_medial + change_medial, data = tdcm_df)

tdcm_change_results <- bind_rows(
  broom::tidy(mard_change) %>% mutate(model = "mard_change"),
  broom::tidy(mod_change) %>% mutate(model = "mod_change"),
  broom::tidy(sidd_change) %>% mutate(model = "sidd_change"),
  broom::tidy(sird_change) %>% mutate(model = "sird_change")) 

tdcm_change_output <- tdcm_change_results %>% 
  mutate(HR = exp(estimate),
         lci = exp(estimate - 1.96 * std.error),
         uci = exp(estimate + 1.96 * std.error),
         model = case_when(model == "mard_change" ~ "MARD",
                           model == "mod_change" ~ "MOD",
                           model == "sidd_change" ~ "SIDD",
                           model == "sird_change" ~ "SIRD",
                           TRUE ~ NA_character_)) %>% 
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
  write_csv("analysis/dsp02_imputed tdcm change results.csv")
#------------------------------------------------------------------------------------------------------------------------
# causal survival forest
# chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://arxiv.org/pdf/2312.02482

# outcome
Y = tdcm_df$time_to_event 
# subgroup indicator
W_mard = tdcm_df$mard
W_mod = tdcm_df$mod
W_sidd = tdcm_df$sidd
W_sird = tdcm_df$sird
# non-censoring indicator
D = tdcm_df$event 
# covariates
X = cbind(
  study_aric = tdcm_df$study_aric,
  study_cardia = tdcm_df$study_cardia,
  study_dppos = tdcm_df$study_dppos,
  study_jhs = tdcm_df$study_jhs,
  study_mesa = tdcm_df$study_mesa,
  race_nhwhi = tdcm_df$race_nhwhi,
  race_nhbla = tdcm_df$race_nhbla,
  race_nhoth = tdcm_df$race_nhoth,
  race_hisp = tdcm_df$race_hisp,
  female = tdcm_df$female,
  serumcreatinine = tdcm_df$serumcreatinine,
  uric_acid = tdcm_df$uric_acid,
  age = tdcm_df$age,
  bmi = tdcm_df$bmi,
  hba1c = tdcm_df$hba1c,
  homa2ir = tdcm_df$homa2ir,
  medial = tdcm_df$medial
)

# set the truncation time to 20 years, at which point most subsequent observations are censored.
hist(Y[D==1],main="HistogramofY", xlab="")
# hist(Y[D ==0],col= adjustcolor("red",0.5),add=TRUE) legend("topright", c("Event","Censored"), col=c("gray",adjustcolor("red",0.5)), lwd=4)
# abline(v=720,lty=2)

library(grf)
csf_mard <- causal_survival_forest(X, Y, W_mard, D, W.hat = 0.5, target = "RMST", horizon = 20)
csf_mod <- causal_survival_forest(X, Y, W_mod, D, W.hat = 0.5, target = "RMST", horizon = 20)
csf_sidd <- causal_survival_forest(X, Y, W_sidd, D, W.hat = 0.5, target = "RMST", horizon = 20)
csf_sird <- causal_survival_forest(X, Y, W_sird, D, W.hat = 0.5, target = "RMST", horizon = 20)

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
  write_csv("analysis/dsp02_imputed csf results.csv")

average_treatment_effect(csf)
# Retrieve out-of-bag CATE estimates 
tau.hat= predict(csf)$predictions 
summary(tau.hat)



## presentation comparing random forest, csf, cox PH ##


