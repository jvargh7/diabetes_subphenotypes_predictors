rm(list = ls());gc();source(".Rprofile")

source("functions/table1_summary.R")

# predictors analysis ---------------------------

# Baseline
baseline_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>%
  # arrange(study, study_id, joint_id, age) %>%
  # group_by(study, study_id, joint_id) %>% 
  # mutate(
  #   tstart = age, 
  #   tstop = dplyr::lead(age)
  # ) %>% 
  # ungroup() %>% 
  # dplyr::filter(age < censored_age) %>% 
  # dplyr::filter((tstart < tstop) & (tstop <= censored_age)) %>%  
  mutate(tg_hdl = tgl/hdlc) %>% 
  group_by(joint_id) %>% 
  mutate(fu_time = censored_age - earliest_age) %>% 
  ungroup() %>% 
  add_count(joint_id, name = "visit_number") %>%  # count total observations per joint_id
  # only baseline
  dplyr::filter(age == earliest_age)


c_vars = c("age","dmagediag","bmi","hba1c","homa2b","homa2ir","sbp","dbp","ldlc","hdlc",
           "insulinf","glucosef","egfr_ckdepi_2021","tgl","tg_hdl","fu_time","visit_number")
p_vars = c("female")
g_vars = c("study","race")


table_df = baseline_df %>% 
  bind_rows(.,
            {.} %>% 
              mutate(subtype="Total")) %>% 
  table1_summary(.,c_vars = c_vars,p_vars = p_vars,g_vars = g_vars,id_vars = "subtype") %>% 
  write_csv(.,"analysis/dspan04_predictors descriptive characteristics by subtype at baseline.csv")


# last FU
lastfu_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>%
  mutate(tg_hdl = tgl/hdlc) %>% 
  # only last follow-up
  dplyr::filter(age == censored_age)


c_vars = c("age","dmagediag","bmi","hba1c","homa2b","homa2ir","sbp","dbp","ldlc","hdlc",
           "insulinf","glucosef","egfr_ckdepi_2021","tgl","tg_hdl")
p_vars = c("female")
g_vars = c("study","race")


table_df = lastfu_df %>% 
  bind_rows(.,
            {.} %>% 
              mutate(subtype="Total")) %>% 
  table1_summary(.,c_vars = c_vars,p_vars = p_vars,g_vars = g_vars,id_vars = "subtype") %>% 
  write_csv(.,"analysis/dspan04_predictors descriptive characteristics by subtype at last follow-up.csv")
