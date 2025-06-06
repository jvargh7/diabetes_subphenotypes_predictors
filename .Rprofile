library(tidyverse)
library(openxlsx)

duration_cutoff <- 1
lab_cutoff <- c(0:365)



if(Sys.info()["user"] == "JVARGH7"){
  path_diabetes_subphenotypes_youth_folder <- "C:/Cloud/OneDrive - Emory University/Papers/Subphenotypes in Youth-onset T2DM"
  path_diabetes_subphenotypes_adults_folder <- "C:/Cloud/Emory University/li, zhongyu - Diabetes Endotypes Project (JV and ZL)"
  cluster_labels = c("MOD","SIRD","MARD","SIDD")
  cluster_colors = c("MOD"="#F8BDA4","SIRD"="#A1C3AC","SIDD"="#ACD9EA","MARD"="#D0ACC9")
  cluster_colors_cosmos = c("MOD"="#F8BDA4","SIRD"="darkgreen","SIDD"="#4682b4","MARD"="#D0ACC9")
  cluster_colors_ada = c("MOD"="#D55E00","SIRD"="darkgreen","SIDD"="#56B4E9","MARD"="#CC79A7")
  
  path_diabetes_subphenotypes_predictors_folder <- "C:/Cloud/OneDrive - Emory University/Papers/Predictors of Subphenotypes"
  
  path_prediabetes_subphenotypes_folder <- "C:/Cloud/OneDrive - Emory University/Papers/Subphenotypes of Prediabetes"
  
}

if(Sys.info()["user"] == "zhongyuli"){
  path_diabetes_subphenotypes_adults_folder <- "/Users/zhongyuli/Library/CloudStorage/OneDrive-EmoryUniversity/Diabetes Endotypes Project (JV and ZL)"
  cluster_labels = c("MOD","SIRD","MARD","SIDD")
  cluster_colors = c("MOD"="#F8BDA4","SIRD"="#A1C3AC","SIDD"="#ACD9EA","MARD"="#D0ACC9")
}


if(Sys.info()["user"] == "JGUO258"){
  cluster_labels = c("MOD","SIRD","MARD","SIDD")
  cluster_colors = c("MOD"="#F8BDA4","SIRD"="#A1C3AC","SIDD"="#ACD9EA","MARD"="#D0ACC9")
  cluster_colors_cosmos = c("MOD"="#F8BDA4","SIRD"="darkgreen","SIDD"="#4682b4","MARD"="#D0ACC9")
  cluster_colors_ada = c("MOD"="#D55E00","SIRD"="darkgreen","SIDD"="#56B4E9","MARD"="#CC79A7")
  
  path_nhanes_ckm_folder <- "C:/Users/JGUO258/OneDrive - Emory/Papers/NHANES CKM Cascade"
  path_g2a_longitudinal_folder <- "C:/Users/JGUO258/OneDrive - Emory/Papers/Crossnational Longitudinal Concordance"
  path_diabetes_subphenotypes_predictors_folder <- "C:/Users/JGUO258/OneDrive - Emory/Papers/Predictors of Subphenotypes"
  path_prediabetes_subphenotypes_folder <- "C:/Users/JGUO258/OneDrive - Emory/Papers/Subphenotypes of Prediabetes"
  path_diabetes_subphenotypes_youth_folder <- "C:/Users/JGUO258/OneDrive - Emory/Papers/Subphenotypes in Youth-onset T2DM"
  path_diabetes_subphenotypes_adults_folder <- "C:/Users/JGUO258/OneDrive - Emory/Papers/Diabetes Endotypes Project (JV and ZL)"
  path_g2a_longitudinal_folder <- "C:/Users/JGUO258/OneDrive - Emory/Papers/Crossnational Longitudinal Concordance"
  path_pasc_cmr_folder <- "C:/Users/JGUO258/OneDrive - Emory/Papers/PASC Cardiometabolic Risk Factors"
  path_pasc_diabetes_folder <- "C:/Users/JGUO258/OneDrive - Emory/Papers/PASC Diabetes Incidence"
  path_pasc_subgroups_folder <- "C:/Users/JGUO258/OneDrive - Emory/Papers/PASC Subgroups"
}
