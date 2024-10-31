rm(list=ls());gc();source(".Rprofile")

library(arrow)
#--------------------------------------------------------------------------------------------------------------
### EHR - OneFlorida ###
#--------------------------------------------------------------------------------------------------------------
demographic <- readRDS(paste0(path_pasc_cmr_folder,"/working/cleaned/pcrpre101_demographic.RDS"))
# index_date <- readRDS(paste0(path_pasc_cmr_folder,"/working/cleaned/pcrpre201_index date.RDS"))
options(arrow.unsafe_metadata = TRUE)
encounter <- read_parquet(paste0(path_pasc_cmr_folder,"/working/raw/encounter_chakkalakal_v1.parquet")) %>% 
  dplyr::select(ID,ENCOUNTERID,ADMIT_DATE) %>% 
  distinct(ENCOUNTERID, .keep_all = TRUE)
nonsupremedm_lastfup <- readRDS(paste0(path_pasc_diabetes_folder,"/working/cleaned/pdpre202_nonsupremedm last followup.RDS"))
supremedm_diagnosis <- readRDS(paste0(path_pasc_subgroups_folder,"/working/pspre01_clinical characteristics for phenotyping.RDS"))

# Non-dm
nonsupremedm <- nonsupremedm_lastfup %>% 
  left_join(demographic %>% 
              dplyr::select(-matchid),
            by = "ID") %>% 
  left_join(encounter,
            by = "ID") %>% 
  group_by(ID) %>% 
  mutate(fuptime = as.numeric(difftime(last_followup_date, min(ADMIT_DATE), units = "days")) / 365.25,
         age_latest = as.numeric(difftime(last_followup_date, index_date, units = "days")) / 365.25 + age) %>% 
  ungroup()

# N obs of non-dm -- number of unique encounter ID, N = 15171513
nonsupremedm %>% 
  distinct(ENCOUNTERID) %>%
  nrow()

# N of non-dm, N = 317043
nonsupremedm %>% 
  distinct(ID) %>%
  nrow()


# DM
supremedm_diagnosis <- readRDS(paste0(path_pasc_subgroups_folder,"/working/pspre01_clinical characteristics for phenotyping.RDS")) %>% 
  left_join(demographic %>% 
              dplyr::select(-matchid),
            by = "ID") %>% 
  mutate(age_diagnosis = as.numeric(difftime(diagnosis_date,index_date,units="days")/365.25) + age) %>% 
  left_join(encounter,
            by = "ID") 


# Among diagnosed DM, duration <= 1 year, N = 12502
supremedm_newdiag <- supremedm_diagnosis %>% 
  group_by(ID) %>% 
  dplyr::filter(!is.na(age_diagnosis) & !is.na(age) & 
                  (age_diagnosis - age) >= 0 & 
                  (age_diagnosis - age) <= 1) %>%
  mutate(fuptime = as.numeric(difftime(diagnosis_date, min(ADMIT_DATE), units = "days")) / 365.25,
         age_latest = age_diagnosis) %>% 
  ungroup()


### OneFlorida Newly diagnosed dm: N = 12502 ###
supremedm_newdiag %>% 
  distinct(ID) %>%
  nrow()

# N obs of dm -- number of unique encounter ID, N = 1218973
supremedm_newdiag %>% 
  distinct(ENCOUNTERID) %>%
  nrow()

## All: dm + non-dm, N = 329545 
supremeall <- bind_rows(supremedm_newdiag %>% 
                          mutate(dm = 1),
                        nonsupremedm %>% 
                          mutate(dm = 0)) 

# distinct(ID) %>%
# nrow()

# female, N = 205618
supremeall_female <- supremeall %>%
  dplyr::filter(female == 1) 

# race minority, N = 164830
supremeall_racemin <- supremeall %>%
  dplyr::filter(nhwhite == 0) 

# mean age_latest
mean(supremeall$age_latest)
sd(supremeall$age_latest)