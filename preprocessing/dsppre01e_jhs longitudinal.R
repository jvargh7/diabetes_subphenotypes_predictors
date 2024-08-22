rm(list = ls());gc();source(".Rprofile")



jhs_newdm <-readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/jhs.RDS")) %>% 
  dplyr::mutate(ratio_th=tgl/hdlc,
                glucosef2=glucosef*0.0555,
                insulinf2=insulinf*6,
                race_rev = "AA")%>% 
  rename(race = race_eth)%>% 
  #new_id = row_number())%>% 
  #rename(study_id = new_id)%>% 
  dplyr::select(study_id,visit,aric,bmi,hba1c,ldlc,hdlc,tgl,sbp,dbp,ratio_th,dmagediag,dmduration,glucosef2,insulinf2,
                serumcreatinine, urinealbumin, urinecreatinine, egfr, totalc,female,race,race_rev) %>% 
  dplyr::filter(dmduration %in% c(0,1),!(aric == 1 & visit == 1))


