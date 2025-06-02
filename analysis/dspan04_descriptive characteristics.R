rm(list = ls());gc();source(".Rprofile")

source("functions/table1_summary.R")

# trajectory analysis ---------------------------

analytic_df = readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>% 
  dplyr::filter(t <= 0 & t >= -15) %>% 
  arrange(joint_id,t) %>% 
  distinct(joint_id,t,.keep_all=TRUE) %>% 
  mutate(tg_hdl = tgl/hdlc)


c_vars = c("age","dmagediag","bmi","hba1c","homa2b","homa2ir","sbp","dbp","ldlc","hdlc",
           "insulinf","glucosef","egfr","tgl","tg_hdl")
p_vars = c("female")
g_vars = c("study","race")


table_df = analytic_df %>% 
  bind_rows(.,
            {.} %>% 
              mutate(subtype="Total")) %>% 
  table1_summary(.,c_vars = c_vars,p_vars = p_vars,g_vars = g_vars,id_vars = "subtype")

write_csv(table_df,"analysis/dspan04_trajectory descriptive characteristics by subtype.csv")


# predictors analysis ---------------------------

analytic_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>%
  arrange(study, study_id, joint_id, age) %>%
  group_by(study, study_id, joint_id) %>% 
  mutate(
    tstart = case_when(row_number() == 1 ~ age, 
                       TRUE ~ dplyr::lag(age, n = 1)), 
    tstop = age
  ) %>% 
  ungroup() %>% 
  dplyr::filter((tstart < tstop) & (tstop <= censored_age))


c_vars = c("age","dmagediag","bmi","hba1c","homa2b","homa2ir","sbp","dbp","ldlc","hdlc",
           "insulinf","glucosef","egfr","tgl","tg_hdl")
p_vars = c("female")
g_vars = c("study","race")


table_df = analytic_df %>% 
  bind_rows(.,
            {.} %>% 
              mutate(subtype="Total")) %>% 
  table1_summary(.,c_vars = c_vars,p_vars = p_vars,g_vars = g_vars,id_vars = "subtype")

write_csv(table_df,"analysis/dspan04_predictors descriptive characteristics by subtype.csv")





