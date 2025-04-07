rm(list=ls());gc();source(".Rprofile")

biomarkers <- readRDS(paste0(path_g2a_longitudinal_folder,"/working/hrs biomarkers.RDS")) %>% 
  dplyr::select(hhid,pn,wave,year,a1c_adj,tc_adj,hdl_adj)

process_hrs_wave <- function(wave_num, path_folder) {
  male_file <- paste0(path_folder, "/working/hrs/G2A HRS Wave ", wave_num, " male.RDS")
  female_file <- paste0(path_folder, "/working/hrs/G2A HRS Wave ", wave_num, " female.RDS")
  
  bind_rows(
    readRDS(male_file), 
    readRDS(female_file)
  ) %>%
    dplyr::select(hhid,pn,hhidpn,age,gender,race,ethnicity,raceeth,
                  diagnosed_dm,agediagnosed_dm,medication_dm,medication_bp,
                  # smoke,smokeever,smokecurr,alcohol,alcohol_days,
                  height,weight,sbp,dbp,waistcircumference) %>%
    mutate(wave = wave_num)
}


hrs_longitudinal <- bind_rows(
  lapply(8:14, process_hrs_wave, path_folder = path_g2a_longitudinal_folder)
) %>% 
  left_join(biomarkers,
            by = c("hhid","pn","wave")) %>% 
  mutate(bmi = weight/(height^2),
         female = case_when(gender == "Female" ~ 1,
                            TRUE ~ 0),
         study_id = as.numeric(hhidpn)) %>% 
  rename(dmagediag = agediagnosed_dm,
         race_eth = raceeth,
         hba1c = a1c_adj,
         totalc = tc_adj,
         hdlc = hdl_adj,
         wc = waistcircumference,
         med_bp_use = medication_bp,
         med_dm_use = medication_dm)

saveRDS(hrs_longitudinal,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01i_hrs.RDS"))
