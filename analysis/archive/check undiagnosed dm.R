
# define anyone in "final_dataset_temp.RDS" as newdm; anyone with missing dmagedia as NoT2D; exclude the rest
# "final_dataset_temp.RDS" consists of all study_newdm.RDS 
# study_newdm is extracted from study_analysis; that means, anyone who's not in study_newdm should be DM free
# but they are not, some of them can be defined as undiagnosed T2D

# JHS -------------------------------------
jhs_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/jhs_newdm.RDS")) # N = 268
jhs_baselinedm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/jhs_baseline_dm.RDS")) 

jhs_analysis <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/jhspre01_jhs_analysis.RDS")) %>% 
  dplyr::filter(aric == 0) %>% 
  dplyr::filter(!study_id %in% jhs_baselinedm$study_id)

jhs_undm <- jhs_analysis %>% 
  dplyr::filter(!study_id %in% jhs_newdm$study_id) %>% 
  dplyr::filter(is.na(dmagediag)) %>% 
  dplyr::filter(hba1c>=6.5 | glucosef>=126)


# ARIC -------------------------------------
aric_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/aric_newdm.RDS")) 
aric_baselinedm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/aric_baseline_dm.RDS")) 

aric_analysis <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/aric_analysis.RDS")) %>%
  dplyr::filter(!study_id %in% aric_baselinedm$study_id) %>%
  dplyr::filter(!is.na(age))

aric_undm <- aric_analysis %>% 
  dplyr::filter(!study_id %in% aric_newdm$study_id) %>% 
  dplyr::filter(is.na(dmagediag)) %>% 
  dplyr::filter(hba1c>=6.5 | glucosef>=126)


# CARDIA -------------------------------------
cardia_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/cardia_newdm.RDS")) 
cardia_baselinedm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/cardia_baseline_dm.RDS")) 

cardia_dat_all <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/cardia_dat_all.RDS")) %>%
  dplyr::filter(!study_id %in% cardia_baselinedm$study_id)

cardia_undm <- cardia_dat_all %>% 
  dplyr::filter(!study_id %in% cardia_newdm$study_id) %>% 
  dplyr::filter(is.na(dmagediag)) %>% 
  dplyr::filter(hba1c>=6.5 | glucosef>=126)


# MESA -------------------------------------
mesa_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/mesa_newdm.RDS")) 
mesa_baselinedm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/mesa_baseline_dm.RDS")) 

mesa_dat_all <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/mesa_dat_all.RDS")) %>% 
  dplyr::filter(!study_id %in% mesa_baselinedm$study_id) %>% 
  dplyr::filter(!is.na(age))

mesa_undm <- mesa_dat_all %>% 
  dplyr::filter(!study_id %in% mesa_newdm$study_id) %>% 
  # dplyr::filter(is.na(dmagediag)) %>% 
  dplyr::filter(hba1c>=6.5 | glucosef>=126)







