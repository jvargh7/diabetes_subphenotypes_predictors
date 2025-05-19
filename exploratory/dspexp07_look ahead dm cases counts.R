rm(list=ls());gc();source(".Rprofile")

#--------------------------------------------------------------------------------------------------------------
### look ahead ### all participants were diagnosed with T2D at baseline
#--------------------------------------------------------------------------------------------------------------
# N = 4901
la_total = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/la_baseline.RDS")) 
# N = 877
la_newdm = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/la_newdm.RDS")) 


# N = 565
la_female <- la_newdm %>%
  dplyr::filter(female == 1)  

# N = 299
la_racemin <- la_newdm %>%
  dplyr::filter(!is.na(race_eth)) %>% 
  dplyr::filter(race_eth != "NH White") 



# longitudinal 
la_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01h_la.RDS"))

la_newdm_long <-la_longitudinal[la_longitudinal$dmduration%in% c(0, 1), ] 
la_newdm_long$study = "look ahead" #N=877 

saveRDS(la_newdm_long, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp07_look ahead new dm.RDS"))
