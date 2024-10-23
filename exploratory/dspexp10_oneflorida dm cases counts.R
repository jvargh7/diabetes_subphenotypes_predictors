rm(list=ls());gc();source(".Rprofile")

#--------------------------------------------------------------------------------------------------------------
### EHR - OneFlorida ###
#--------------------------------------------------------------------------------------------------------------
# all diagnosed dm, N = 14948
demographic <- readRDS(paste0(path_pasc_cmr_folder,"/working/cleaned/pcrpre101_demographic.RDS"))
supremedm_diagnosis <- readRDS(paste0(path_pasc_subgroups_folder,"/working/pspre01_clinical characteristics for phenotyping.RDS")) %>% 
  left_join(demographic %>% 
              dplyr::select(-matchid),
            by = "ID") %>% 
  mutate(age_diagnosis = as.numeric(difftime(diagnosis_date,index_date,units="days")/365.25) + age)


# Among diagnosed DM, duration <= 1 year, N = 12502
supremedm_newdiag <- supremedm_diagnosis %>% 
  group_by(ID) %>% 
  dplyr::filter(!is.na(age_diagnosis) & !is.na(age) & 
                  (age_diagnosis - age) >= 0 & 
                  (age_diagnosis - age) <= 1) %>%
  ungroup()


### OneFlorida Newly diagnosed dm: 12502 ###

# distinct(ID) %>%
# nrow()

# female, N = 6761
supremedm_female <- supremedm_newdiag %>%
  dplyr::filter(female == 1) 

# race minority, N = 7043
supremedm_racemin <- supremedm_newdiag %>%
  dplyr::filter(nhwhite == 0) 