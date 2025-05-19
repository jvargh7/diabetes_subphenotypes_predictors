rm(list=ls());gc();source(".Rprofile")

#--------------------------------------------------------------------------------------------------------------
### ACCORD ### all participants were diagnosed with T2D at baseline
#--------------------------------------------------------------------------------------------------------------
# N = 10251
accord_total = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/accord_baseline.RDS")) 
# N = 601
accord_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/accord_newdm.RDS")) 


# N = 235
accord_female <- accord_newdm %>%
  dplyr::filter(female == 1)  

# N = 237
accord_racemin <- accord_newdm %>%
  dplyr::filter(!is.na(race_eth)) %>% 
  dplyr::filter(race_eth != "NH White") 

# longitudinal 
accord_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01g_accord.RDS"))

accord_newdm_long <-accord_longitudinal[accord_longitudinal$dmduration%in% c(0, 1), ] 
accord_newdm_long$study = "accord" #N=601 

saveRDS(accord_newdm_long, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp01_accord new dm.RDS"))

