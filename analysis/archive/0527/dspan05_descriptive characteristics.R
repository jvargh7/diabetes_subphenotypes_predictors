rm(list = ls());gc();source(".Rprofile")

source("functions/table1_summary.R")

# trajectory analysis ---------------------------

wave_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan04_trajectory analysis analytic df.RDS")) %>% 
  mutate(tg_hdl = tgl/hdlc)


c_vars = c("age","dmagediag","bmi","hba1c","homa2b","homa2ir","sbp","dbp","ldlc","hdlc",
           "insulinf","glucosef","egfr","tgl","tg_hdl")
p_vars = c("female")
g_vars = c("study","race")


table_df = wave_df %>% 
  bind_rows(.,
            {.} %>% 
              mutate(subtype="Total")) %>% 
  table1_summary(.,c_vars = c_vars,p_vars = p_vars,g_vars = g_vars,id_vars = "subtype")

write_csv(table_df,"analysis/dspan05_trajectory descriptive characteristics by subtype.csv")


# predictors analysis ---------------------------

ipcw_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan02_predictors analysis ipcw df.RDS")) %>% 
  mutate(subtype = case_when(is.na(dmagediag) ~ "NOT2D",
                             !is.na(cluster) ~ cluster,
                             TRUE ~ NA_character_),
         tg_hdl = tgl/hdlc) %>% 
  dplyr::filter(!is.na(subtype))

c_vars = c("age","dmagediag","bmi","hba1c","homa2b","homa2ir","sbp","dbp","ldlc","hdlc",
           "insulinf","glucosef","egfr","tgl","tg_hdl")
p_vars = c("female")
g_vars = c("study","race")


table_df = ipcw_df %>% 
  bind_rows(.,
            {.} %>% 
              mutate(subtype="Total")) %>% 
  table1_summary(.,c_vars = c_vars,p_vars = p_vars,g_vars = g_vars,id_vars = "subtype")

write_csv(table_df,"analysis/dspan05_predictors descriptive characteristics by subtype.csv")





