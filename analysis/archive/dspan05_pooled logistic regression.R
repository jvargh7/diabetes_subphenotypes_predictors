rm(list = ls());gc();source(".Rprofile")

library(survival)

options(scipen = 999)

# participants who develop diabetes after baseline, or never develop diabetes
analytic_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01_analytic df.RDS")) %>% 
  mutate(dm = case_when(!is.na(dmagediag) | hba1c >= 6.5 | glucosef >= 126 | glucose2h >= 200 ~ 1,
                        TRUE ~ 0))
# dm: 8594, non-dm: 17820


anthro_vars <- c("sbp", "dbp", "height", "wc", "bmi")
lab_vars <- c("hba1c", "insulinf", "glucosef", "glucose2h", "tgl", "hdlc", "ldlc", 
              "serumcreatinine", "urinecreatinine", "egfr", "apo_a", "apo_b", "uric_acid", "homa2b", "homa2ir")

# Create a pooled logistic regression model
model <- glm(dm ~ factor(visit) + sbp + dbp + bmi + hba1c + hdlc + ldlc + homa2b + homa2ir
             + serumcreatinine + urinecreatinine,
             family = binomial(link = "logit"), data = analytic_df)


summary(model)

model_summary <- broom::tidy(model)

# height + wc + glucose2h + tgl + hdlc + ldlc +
# egfr + apo_a + apo_b + uric_acid







