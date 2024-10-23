rm(list=ls());gc();source(".Rprofile")

#--------------------------------------------------------------------------------------------------------------
## Youth ##
#--------------------------------------------------------------------------------------------------------------
# N = 307, obs = 2834
search <- readRDS(paste0(path_diabetes_subphenotypes_youth_folder,"/working/search/search_etiologic.RDS"))
# N = 334, obs = 699
today <- readRDS(paste0(path_diabetes_subphenotypes_youth_folder,"/working/today/today_baseline.RDS")) 