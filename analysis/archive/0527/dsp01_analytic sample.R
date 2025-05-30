rm(list = ls());gc();source(".Rprofile")

library(purrr)
library(readr)
library(writexl)

final_dataset_temp = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/final_dataset_temp.RDS"))
aric_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01a_aric.RDS")) %>% 
  mutate(study = "aric")
cardia_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01b_cardia.RDS")) %>% 
  mutate(study = "cardia")
jhs_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01e_jhs.RDS")) %>% 
  mutate(study = "jhs",
         race = "NH Black")
dppos_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01c_dppos.RDS")) %>% 
  mutate(study = "dppos",
         race = race_eth)
mesa_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01f_mesa.RDS")) %>% 
  mutate(study = "mesa")

clusters = read_csv(paste0(path_diabetes_subphenotypes_adults_folder,"/working/processed/dec_an02_clean_kmeans_5var_mi_knn_cluster.csv")) %>% 
  dplyr::select(-one_of("...1")) %>% 
  left_join(final_dataset_temp %>% 
              dplyr::select(study_id,original_study_id),
            by=c("study_id")) %>% 
  rename(cluster_study_id = study_id)

longitudinal_df = bind_rows(aric_longitudinal %>% mutate(study_id = as.numeric(str_replace(study_id,"C",""))),
                            cardia_longitudinal,
                            jhs_longitudinal,
                            dppos_longitudinal %>% mutate(female = sex - 1),
                            mesa_longitudinal) %>%
  
  bind_rows(final_dataset_temp %>% dplyr::select(-study_id) %>% 
              rename(study_id = original_study_id) %>% 
              mutate(age = case_when(is.na(age) ~ dmagediag,
                                     TRUE ~ age)) %>% 
              dplyr::select(one_of(unique(c(
                names(aric_longitudinal),
                names(cardia_longitudinal),
                names(jhs_longitudinal),
                names(dppos_longitudinal),
                names(mesa_longitudinal))
                
              )
              ))) %>% 
  
  mutate(glucosef2 = case_when(is.na(glucosef2) ~ glucosef*0.0555,
                               TRUE ~ glucosef2),
         # 1 μIU/mL = 6.00 pmol/L
         insulinf2 = case_when(is.na(insulinf2) ~ insulinf*6,
                               TRUE ~ insulinf2)) %>% 
  mutate(glucosef = case_when(is.na(glucosef) ~ glucosef2/0.0555,
                              TRUE ~ glucosef),
         insulinf = case_when(is.na(insulinf) ~ insulinf2/6,
                              TRUE ~ insulinf),
         weight = case_when(!is.na(height) ~ bmi*(height/100)^2,
                            TRUE ~ NA_real_)
  ) %>% 
  mutate(study = case_when(study == "dpp" ~ "dppos",
                           TRUE ~ study))

saveRDS(longitudinal_df,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01_longitudinal df.RDS"))

#--------------------------------------------------------------------------------------------------------------------
# one cohort one sheet; keep study_id, study, age
homa2_df <- longitudinal_df %>% 
  dplyr::select(study_id,study,age,glucosef2,insulinf2) %>% 
  dplyr::filter(!is.na(glucosef2) & !is.na(insulinf2)) %>% 
  mutate(
    glucosef2 = case_when(
      glucosef2 < 3 ~ 3,
      glucosef2 > 25 ~ 25,
      TRUE ~ glucosef2
    ),
    insulinf2 = case_when(
      insulinf2 < 20 ~ 20,
      insulinf2 > 400 ~ 400,
      TRUE ~ insulinf2
    ))

study_names <- c("aric","cardia","jhs","dppos","mesa")
list_of_dataframes <- lapply(study_names, function(name) {
  homa2_df %>% 
    dplyr::filter(study == name)
})

names(list_of_dataframes) <- study_names


excel_path <- paste0(path_diabetes_subphenotypes_predictors_folder, "/working/cleaned/homa2 calculation/homa2 indices calculation df.xlsx")
write_xlsx(list_of_dataframes, excel_path)

#--------------------------------------------------------------------------------------------------------------------
# add homa2 in the dataset by study_id, study, age
path_to_file <- paste0(path_diabetes_subphenotypes_predictors_folder, "/working/cleaned/homa2 calculation/homa2 indices values.xlsx")
sheets <- c("aric", "cardia", "jhs", "dppos", "mesa")
list_of_data <- lapply(sheets, function(sheet) readxl::read_excel(path_to_file, sheet = sheet))

homa2_combined <- bind_rows(list_of_data, .id = "study")
homa2_combined$study <- rep(sheets, sapply(list_of_data, nrow))

analytic_df <- longitudinal_df %>% 
  left_join(homa2_combined %>% 
              rename(homa2b = `HOMA2 %B`,
                     homa2ir = `HOMA2 IR`) %>% 
              dplyr::select(-c(`HOMA2 %S`,glucosef2,insulinf2)) %>% 
              distinct(study_id,study,age,.keep_all=TRUE),
            by = c("study_id","study","age")) %>% 
  left_join(clusters %>% 
              dplyr::select(cluster_study_id,original_study_id,cluster,study,female),
            by=c("study"="study","study_id" = "original_study_id","female"))

saveRDS(analytic_df,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dsppre01_analytic df.RDS"))

