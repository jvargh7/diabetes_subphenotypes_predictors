rm(list = ls());gc();source(".Rprofile")



aric_newdm <- readRDS(paste0(path_endotypes_folder,"/working/cleaned/aric_newdm.RDS"))%>% 
  dplyr::mutate(ratio_th=tgl/hdlc,
                glucosef2=glucosef*0.0555,
                insulinf2=insulinf*6,
                urinealbumin = urinealbumin/10,
                race_rev = case_when(
                  race == "W" ~ "White",
                  race == "B" ~ "AA",
                  TRUE ~ NA_character_  
                ),
                female = case_when(
                  female == "M" ~ 0,  
                  female == "F" ~ 1,  
                  TRUE ~ NA_integer_
                ),
                race = case_when(
                  race == "W" ~ "NH White",
                  race == "B" ~ "NH Black",
                  TRUE ~ NA_character_  
                ),
  )%>% 
  dplyr::select(bmi,hba1c,ldlc,hdlc,tgl,sbp,dbp,ratio_th,dmagediag,glucosef2,insulinf2,
                serumcreatinine,urinealbumin,totalc,female,race,race_rev)