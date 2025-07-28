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
p_vars = c("female","med_bp_use","med_dep_use","med_chol_use")
g_vars = c("study","race","smoking")


table_df = baseline_df %>% 
  bind_rows(.,
            {.} %>% 
              mutate(subtype="Total")) %>% 
  table1_summary(.,c_vars = c_vars,p_vars = p_vars,g_vars = g_vars,id_vars = "subtype") %>% 
  write_csv(.,"analysis/dspan04_predictors descriptive characteristics by subtype at baseline.csv")

# check N (%) -------
lastfu_df %>%
  group_by(med_chol_use,subtype) %>%
  summarise(
    n_joint_id = n_distinct(joint_id),
    .groups = "drop"
  ) 



# last FU
lastfu_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>%
  group_by(joint_id) %>% 
  mutate(fu_time = censored_age - earliest_age) %>% 
  ungroup() %>% 
  mutate(tg_hdl = tgl/hdlc) %>% 
  # only last follow-up
  dplyr::filter(age == censored_age)


c_vars = c("age","dmagediag","bmi","hba1c","homa2b","homa2ir","sbp","dbp","ldlc","hdlc",
           "insulinf","glucosef","egfr_ckdepi_2021","tgl","tg_hdl","fu_time")
p_vars = c("female","med_bp_use","med_dep_use","med_chol_use")
g_vars = c("study","race","smoking")


table_df = lastfu_df %>% 
  bind_rows(.,
            {.} %>% 
              mutate(subtype="Total")) %>% 
  table1_summary(.,c_vars = c_vars,p_vars = p_vars,g_vars = g_vars,id_vars = "subtype") %>% 
  write_csv(.,"analysis/dspan04_predictors descriptive characteristics by subtype at last follow-up.csv")




# Percentage of missing values by subtype ---------------------------------------

clean_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS"))

vars <- c("hba1c", "bmi", "homa2b", "homa2ir", "sbp", "sbp", "ldlc", "hdlc", 
          "insulinf","glucosef","egfr_ckdepi_2021","tgl")


na_summary <- clean_df %>%
  select(subtype, all_of(vars)) %>%
  mutate(subtype = ifelse(is.na(subtype), "Missing", as.character(subtype))) %>%
  pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "value") %>%
  group_by(variable, subtype) %>%
  summarise(
    na_percent = mean(is.na(value)) * 100,
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = subtype, values_from = na_percent) %>%
  # Add overall NA % column
  left_join(
    clean_df %>%
      select(all_of(vars)) %>%
      pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
      group_by(variable) %>%
      summarise(Overall = mean(is.na(value)) * 100, .groups = "drop"),
    by = "variable"
  ) %>%
  # Reorder columns
  select(variable, Overall, NOT2D, MARD, MOD, SIDD, SIRD) %>%
  # Format % and relabel variables
  mutate(across(-variable, ~sprintf("%.1f%%", .))) %>% 
  write.csv(.,"analysis/dspan04_percentage of missing values by subtype.csv")













