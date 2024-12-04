rm(list = ls());gc();source(".Rprofile")

library(survival)
library(survminer)
library(ggsurvfit)
library(broom)
library(tidyr)

# follow-up time >15y & dm == 0 --- event == 0
mice_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/mi_dfs.RDS")) %>% 
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
                                TRUE ~ 0)) %>%
  # define time to event
  mutate(time_to_event = censored_age - age) %>% 
  dplyr::filter(time_to_event > 0)


#------------------------------------------------------------------------------------------------------------------------
# cross-sectional datatset, cox PH model
# overall, include all variables, 1 obs for each person
cross_df <- mice_df %>% 
  group_by(study_id,study) %>% 
  dplyr::filter(age == min(age)) %>% 
  ungroup()

table(cross_df$event)

# Visualize distribution of time to event or censoring
ggplot(data=cross_df,aes(x=time_to_event,group=event,fill=factor(event))) +
  geom_histogram(alpha=0.5,position=position_dodge(width=0.9),bins=10)


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
cluster_df <- mice_df %>% 
  left_join(clusters %>% 
              dplyr::select(cluster_study_id,original_study_id,cluster,study),
            by=c("study"="study","study_id" = "original_study_id", "cluster_study_id")) %>% 
  ### exclude diabetes patients with no cluster data
  dplyr::filter(!(event == 1 & is.na(cluster)))


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

saveRDS(analytic_sample, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/analytic_df.RDS"))

diseases <- c("mard", "mod", "sidd", "sird")


# overall rates
analytic_sample %>% 
  summarize(across(diseases, ~mean(.)))

#--------------------------------------------------------------------------------------------------------------------------
# variable selection - LASSO, SCAD
# chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cran.r-project.org/web/packages/ncvreg/ncvreg.pdf

# exclude missing values
variables_to_check <- c("study_aric", "study_cardia", "study_dppos", "study_jhs", "study_mesa",
                        "race_nhwhi", "race_nhbla", "race_nhoth", "race_hisp",
                        "age", "female", "bmi", "hba1c", "sbp", "dbp", "hdlc", "ldlc",
                        "homa2b", "homa2ir", "height", "wc", "glucosef", "glucose2h", "tgl",
                        "serumcreatinine", "urinecreatinine", "egfr", "apo_a", "apo_b", "uric_acid",
                        "vldlc", "hc", "triceps", "iliac", "medial", "ast", "alt", "urinealbumin",
                        "weight", "insulinf", "abdominal")


df_clean <- analytic_sample %>%
  dplyr::filter(complete.cases(.[variables_to_check])) %>%
  # require all variables to be numeric
  mutate(female = as.numeric(as.character(female)))

y <- Surv(df_clean$time_to_event, df_clean$event)  

# Select predictors
X <- as.matrix(df_clean[, c("study_aric", "study_cardia", "study_dppos", "study_jhs", "study_mesa",
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

# non-zero: study, homa2ir, hba1c, race_eth, female, serumcreatinine, uric_acid, bmi, insulinf, age, hc, glucosef
# glucose2h, homa2b, wc, height, hdlc, vldlc, dbp, triceps, alt, sbp, egfr, urinealbumin, urinecreatinine
#---------------------------------------------------------------------------------------------------------------------------------
# cox PH model
cross_df <- analytic_sample %>% 
  group_by(study_id,study) %>% 
  dplyr::filter(age == min(age)) %>% 
  ungroup()

table(cross_df$event)

mard_cp <- coxph(Surv(time_to_event, mard) ~ study_aric + study_cardia + study_dppos + study_jhs + study_mesa 
                 + race_nhwhi + race_nhbla + race_nhoth + race_hisp 
                 + homa2ir + hba1c + female + serumcreatinine + uric_acid + bmi + insulinf + age + hc + glucosef, 
                 data = cross_df)

mod_cp <- coxph(Surv(time_to_event, mod) ~ study_aric + study_cardia + study_dppos + study_jhs + study_mesa 
                + race_nhwhi + race_nhbla + race_nhoth + race_hisp 
                + homa2ir + hba1c + female + serumcreatinine + uric_acid + bmi + insulinf + age + hc + glucosef, 
                data = cross_df)

sidd_cp <- coxph(Surv(time_to_event, sidd) ~ study_aric + study_cardia + study_dppos + study_jhs + study_mesa 
                 + race_nhwhi + race_nhbla + race_nhoth + race_hisp 
                 + homa2ir + hba1c + female + serumcreatinine + uric_acid + bmi + insulinf + age + hc + glucosef, 
                 data = cross_df)

sird_cp <- coxph(Surv(time_to_event, sird) ~ study_aric + study_cardia + study_dppos + study_jhs + study_mesa 
                 + race_nhwhi + race_nhbla + race_nhoth + race_hisp 
                 + homa2ir + hba1c + female + serumcreatinine + uric_acid + bmi + insulinf + age + hc + glucosef, 
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
  write_csv("analysis/dsp02_imputed cox ph results.csv")
#------------------------------------------------------------------------------------------------------------------------
# TDCM analysis
# Start time: previous age; End time: current age at visit

tdcm_df <- analytic_sample %>%
  arrange(study, study_id, age) %>%
  group_by(study, study_id) %>%
  mutate(
    # Creating baseline variables by using the `lag()` function to get the previous age's value
    # Baseline values
    baseline_homa2ir = case_when(
      row_number() == 1 ~ homa2ir,
      TRUE ~ dplyr::lag(homa2ir, 1, default = NA)
    ),
    baseline_hba1c = case_when(
      row_number() == 1 ~ hba1c,
      TRUE ~ dplyr::lag(hba1c, 1, default = NA)
    ),
    baseline_serumcreatinine = case_when(
      row_number() == 1 ~ serumcreatinine,
      TRUE ~ dplyr::lag(serumcreatinine, 1, default = NA)
    ),
    baseline_uric_acid = case_when(
      row_number() == 1 ~ uric_acid,
      TRUE ~ dplyr::lag(uric_acid, 1, default = NA)
    ),
    baseline_bmi = case_when(
      row_number() == 1 ~ bmi,
      TRUE ~ dplyr::lag(bmi, 1, default = NA)
    ),
    baseline_insulinf = case_when(
      row_number() == 1 ~ insulinf,
      TRUE ~ dplyr::lag(insulinf, 1, default = NA)
    ),
    baseline_hc = case_when(
      row_number() == 1 ~ hc,
      TRUE ~ dplyr::lag(hc, 1, default = NA)
    ),
    baseline_glucosef = case_when(
      row_number() == 1 ~ glucosef,
      TRUE ~ dplyr::lag(glucosef, 1, default = NA)
    ),
    
    # Change calculations
    change_homa2ir = homa2ir - baseline_homa2ir,
    change_hba1c = hba1c - baseline_hba1c,
    change_serumcreatinine = serumcreatinine - baseline_serumcreatinine,
    change_uric_acid = uric_acid - baseline_uric_acid,
    change_bmi = bmi - baseline_bmi,
    change_insulinf = insulinf - baseline_insulinf,
    change_hc = hc - baseline_hc,
    change_glucosef = glucosef - baseline_glucosef
  ) %>%
  mutate(
    tstart = case_when(row_number() == 1 ~ age, 
                       TRUE ~ dplyr::lag(age, n = 1)), 
    tstop = age
  ) %>%
  ungroup() %>% 
  dplyr::filter(tstart < tstop)
  # dplyr::filter((tstart < tstop) & (tstop <= censored_age))

# check
test <- tdcm_df %>% select(study, study_id, age, serumcreatinine,baseline_serumcreatinine,change_serumcreatinine, tstart,tstop)

# collinearity
cor(tdcm_df[, c("homa2ir", "hba1c", "serumcreatinine", "uric_acid", "bmi", "insulinf", "hc", "glucosef")]) 


# exclude correlation coefficient > 0.7 -- highly correlated
# hba1c ~ glucosef, glucose 2h ~ glucosef, sbp ~ dbp, hc ~ bmi

tdcm_mod <- coxph(Surv(tstart, tstop, event) ~ study + race_eth + homa2ir + hba1c + female + serumcreatinine + uric_acid + bmi 
                  + insulinf + age + hc + glucosef, data = tdcm_df)

mard_tdcm <- coxph(Surv(tstart, tstop, mard) ~ study + race_eth + homa2ir + hba1c + female + serumcreatinine + uric_acid + bmi 
                   + insulinf + age + hc + glucosef, data = tdcm_df)
mod_tdcm <- coxph(Surv(tstart, tstop, mod) ~ study + race_eth + homa2ir + hba1c + female + serumcreatinine + uric_acid + bmi 
                  + insulinf + age + hc + glucosef, data = tdcm_df)
sidd_tdcm <- coxph(Surv(tstart, tstop, sidd) ~ study + race_eth + homa2ir + hba1c + female + serumcreatinine + uric_acid + bmi 
                   + insulinf + age + hc + glucosef, data = tdcm_df)
sird_tdcm <- coxph(Surv(tstart, tstop, sird) ~ study + race_eth + homa2ir + hba1c + female + serumcreatinine + uric_acid + bmi 
                   + insulinf + age + hc + glucosef, data = tdcm_df)

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
  write_csv("analysis/dsp02_imputed tdcm results.csv")
#------------------------------------------------------------------------------------------------------------------------
# TDCM change analysis - baseline + change
 
tdcm_change <- coxph(Surv(tstart, tstop, event) ~ study + age + race_eth + female 
                     + baseline_homa2ir + baseline_hba1c + baseline_serumcreatinine + 
                       baseline_uric_acid + baseline_bmi + baseline_insulinf + 
                       baseline_hc + baseline_glucosef
                     + change_homa2ir + change_hba1c + change_serumcreatinine + 
                       change_uric_acid + change_bmi + change_insulinf + 
                       change_hc + change_glucosef, data = tdcm_df)

mard_change <- coxph(Surv(tstart, tstop, mard) ~ study + age + race_eth + female 
                     + baseline_homa2ir + baseline_hba1c + baseline_serumcreatinine + 
                       baseline_uric_acid + baseline_bmi + baseline_insulinf + 
                       baseline_hc + baseline_glucosef
                     + change_homa2ir + change_hba1c + change_serumcreatinine + 
                       change_uric_acid + change_bmi + change_insulinf + 
                       change_hc + change_glucosef, data = tdcm_df)
mod_change <- coxph(Surv(tstart, tstop, mod) ~ study + age + race_eth + female 
                    + baseline_homa2ir + baseline_hba1c + baseline_serumcreatinine + 
                      baseline_uric_acid + baseline_bmi + baseline_insulinf + 
                      baseline_hc + baseline_glucosef
                    + change_homa2ir + change_hba1c + change_serumcreatinine + 
                      change_uric_acid + change_bmi + change_insulinf + 
                      change_hc + change_glucosef, data = tdcm_df)
sidd_change <- coxph(Surv(tstart, tstop, sidd) ~ study + age + race_eth + female 
                     + baseline_homa2ir + baseline_hba1c + baseline_serumcreatinine + 
                       baseline_uric_acid + baseline_bmi + baseline_insulinf + 
                       baseline_hc + baseline_glucosef
                     + change_homa2ir + change_hba1c + change_serumcreatinine + 
                       change_uric_acid + change_bmi + change_insulinf + 
                       change_hc + change_glucosef, data = tdcm_df)
sird_change <- coxph(Surv(tstart, tstop, sird) ~ study + age + race_eth + female 
                     + baseline_homa2ir + baseline_hba1c + baseline_serumcreatinine + 
                       baseline_uric_acid + baseline_bmi + baseline_insulinf + 
                       baseline_hc + baseline_glucosef
                     + change_homa2ir + change_hba1c + change_serumcreatinine + 
                       change_uric_acid + change_bmi + change_insulinf + 
                       change_hc + change_glucosef, data = tdcm_df)

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
  write_csv("analysis/dsp02_imputed tdcm change results.csv")
#------------------------------------------------------------------------------------------------------------------------
# causal survival forest
# chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://arxiv.org/pdf/2312.02482

analytic_sample <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/analytic_df.RDS"))

# outcome
Y = analytic_sample$time_to_event 
# subgroup indicator
W_mard = analytic_sample$mard
W_mod = analytic_sample$mod
W_sidd = analytic_sample$sidd
W_sird = analytic_sample$sird
# non-censoring indicator
D = analytic_sample$event 
# covariates
X = cbind(
  study_aric = analytic_sample$study_aric,
  study_cardia = analytic_sample$study_cardia,
  study_dppos = analytic_sample$study_dppos,
  study_jhs = analytic_sample$study_jhs,
  study_mesa = analytic_sample$study_mesa,
  race_nhwhi = analytic_sample$race_nhwhi,
  race_nhbla = analytic_sample$race_nhbla,
  race_nhoth = analytic_sample$race_nhoth,
  race_hisp = analytic_sample$race_hisp,
  female = analytic_sample$female,
  age = analytic_sample$age,
  homa2ir = analytic_sample$homa2ir,
  hba1c = analytic_sample$hba1c,
  serumcreatinine = analytic_sample$serumcreatinine,
  uric_acid = analytic_sample$uric_acid,
  bmi = analytic_sample$bmi
  # insulinf = analytic_sample$insulinf,
  # hc = analytic_sample$hc,
  # glucosef = analytic_sample$glucosef
)

failure_times <- seq(0, 15, by=0.5)

# set the truncation time to 20 years, at which point most subsequent observations are censored.
hist(Y[D==1],main="HistogramofY", xlab="")
hist(Y[D==0], col=adjustcolor("red", 0.5), add=TRUE);
legend("topright", c("Event", "Censored"), fill=c("gray", adjustcolor("red", 0.5)), lwd=4)

abline(v=720,lty=2)

library(grf)
csf_mard <- causal_survival_forest(X, Y, W_mard, D, W.hat = 0.5, target = "RMST", horizon = 15, failure.times = failure_times)
csf_mod <- causal_survival_forest(X, Y, W_mod, D, W.hat = 0.5, target = "RMST", horizon = 15)
csf_sidd <- causal_survival_forest(X, Y, W_sidd, D, W.hat = 0.5, target = "RMST", horizon = 15)
csf_sird <- causal_survival_forest(X, Y, W_sird, D, W.hat = 0.5, target = "RMST", horizon = 15)

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

