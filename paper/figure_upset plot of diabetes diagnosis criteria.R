rm(list = ls());gc();source(".Rprofile")

library(haven)
library(UpSetR)

final_dataset_temp = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/final_dataset_temp.RDS")) 

final_dataset_6c_clean = read.csv(paste0(path_diabetes_subphenotypes_adults_folder,"/working/processed/final_dataset_6c_clean_mi_imputed.csv")) %>%
  left_join(final_dataset_temp %>% select(study, study_id, original_study_id), by = c("study","study_id")) %>% 
  # study_id is the row number
  select(-study_id) %>% 
  mutate(study = case_when(study == "dpp" ~ "dppos",
                           TRUE ~ study)) %>% 
  rename(study_id = original_study_id) %>% 
  mutate(joint_id = paste(study, study_id, sep = "_"))


t2d_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>% 
  dplyr::filter(newdm_event == 1) %>% 
  distinct(joint_id)
  
# N = 1,569
t2d_criteria <- final_dataset_6c_clean %>% 
  dplyr::filter(joint_id %in% t2d_df$joint_id) %>% 
  select(joint_id, report_dx, a1c_dx, glucf_dx, gluc2h_dx) %>% 
  rename(
    Self_report = report_dx,
    Elevated_HbA1c = a1c_dx,
    Elevated_FPG = glucf_dx,
    Elevated_2hOGTT = gluc2h_dx
  )


# Mutually exclusive Upset Plot ---------------------------------------------

# Plot & save with a graphics device (UpSetR does not return a ggplot object)
outfile <- file.path(path_diabetes_subphenotypes_predictors_folder,
                     "figures", "upset_plot_diabetes_diagnosis_criteria.png")

png(outfile, width = 2400, height = 1800, res = 300)
UpSetR::upset(
  t2d_criteria,
  sets = c("Self_report", "Elevated_HbA1c", "Elevated_FPG", "Elevated_2hOGTT"),
  order.by = "freq",
  nintersects = NA,
  sets.bar.color = "#6baed6",
  main.bar.color = "#2171b5",
  text.scale = c(1.3, 1.3, 1.1, 1.1, 1.3, 1.2),
  mb.ratio = c(0.55, 0.45),
  keep.order = TRUE,
  number.angles = 0
)
dev.off()

