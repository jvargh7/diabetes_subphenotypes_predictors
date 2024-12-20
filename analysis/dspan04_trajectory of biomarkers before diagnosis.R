rm(list = ls());gc();source(".Rprofile")

source("functions/egfr_ckdepi_2021.R")

# 1. Restrict to individuals with newly diagnosed T2D or no T2D. Set max_age as dmagediag if new T2D or max(age) if no T2D

# 2. Take their data retrospectively up to 15 years (i.e. max_age - age <= 15)
# Create a 't' variable: max_age - age
# Retain the wave of max_age (t = 0)

# 3. Include individuals who have 
# - at least one wave before T2D diagnosis (if newly diagnosed) AND have cluster info  
# - OR have at least two waves (if no T2D)
# Create an exposure variable called 'subtype' with 5 categories: SIDD, SIRD, MOD, MARD and NoT2D


# 4. Fit a model that accounts for clustering
# study_id: Variable that is a unique ID for each individual

# 5. Estimate predicted trajectories for each subtype over time along with confidence intervals
# https://stats.stackexchange.com/questions/109369/how-can-i-estimate-model-predicted-means-a-k-a-marginal-means-lsmeans-or-em

# new dm
final_dataset_temp = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/final_dataset_temp.RDS")) %>% 
  mutate(joint_id = paste(study, study_id, sep = "_"))

analytic_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dsppre01_analytic df.RDS")) %>% 
  mutate(egfr_ckdepi_2021 = egfr_ckdepi_2021(scr = serumcreatinine,female = female,age = age),
         joint_id = paste(study, study_id, sep = "_"),
         # newly diagnosed T2D or no T2D
         newdm_event = case_when(joint_id %in% final_dataset_temp$joint_id ~ 1,
                           is.na(dmagediag) ~ 0,
                           TRUE ~ NA)) %>%
  dplyr::filter(!is.na(newdm_event)) %>% 
  group_by(study,study_id,joint_id) %>% 
  mutate(max_age = case_when(newdm_event == 1 ~ dmagediag,
                             TRUE ~ max(age)),
         t = max_age - age) %>% 
  dplyr::filter(t >= 0 & t <= 15) %>% 
  ungroup()

wave_df <- analytic_df %>% 
  group_by(study,study_id,joint_id) %>%
  mutate(has_age_before_max = any(age < max_age)) %>% 
  dplyr::filter((newdm_event == 1 & has_age_before_max & !is.na(cluster)) 
                | (newdm_event == 0 & has_age_before_max)) %>% 
  mutate(subtype = case_when(is.na(dmagediag) ~ "NOT2D",
                             !is.na(cluster) ~ cluster,
                             TRUE ~ NA_character_)) %>% 
  ungroup()
  

library(splines)
library(geepack)
m1 = geeglm(homa2b ~ subtype*ns(t) + max_age + female + study + race, family = gaussian(),data = wave_df,id = joint_id,corstr = "exchangeable")
m2 = geeglm(bmi ~ subtype*ns(t) + max_age + female + study + race, family = gaussian(),data = wave_df,id = joint_id,corstr = "exchangeable")
m3 = geeglm(hba1c ~ subtype*ns(t) + max_age + female + study + race, family = gaussian(),data = wave_df,id = joint_id,corstr = "exchangeable")

































































