rm(list = ls());gc();source(".Rprofile")

library(mice)
library(survival)
library(survminer)
library(ggsurvfit)
library(broom)
library(boot)

source("functions/egfr_ckdepi_2021.R")
mi_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/mi_dfs_new.RDS"))

clean_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>% 
  dplyr::filter(dpp_intervention == 1) %>% 
  distinct(study,study_id,joint_id,dpp_intervention) # n = 60

mi_dfs_data = mi_dfs$data 

elevated_a1c = mi_dfs_data %>% 
  dplyr::filter(is.na(cluster)) %>% 
  dplyr::filter(hba1c >=6.5)

length(unique(elevated_a1c$study_id))


mi_dfs_data %>% 
  dplyr::filter(is.na(cluster)) %>% 
  # dplyr::filter(!study_id %in% elevated_a1c$study_id) %>% 
  nrow()

mi_dfs_data %>% 
  dplyr::filter(is.na(cluster)) %>% 
  dplyr::filter(!study_id %in% elevated_a1c$study_id) %>% 
  nrow()
