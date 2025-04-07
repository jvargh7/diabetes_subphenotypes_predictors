rm(list=ls());gc();source(".Rprofile")

# cleaned longitudinal data
pooled_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dspexp13_final longitudinal data.RDS")) 

#--------------------------------------------------------------------------------------------------------------------
# add homa2 in the dataset by study_id, study, age
path_to_file <- paste0(path_diabetes_subphenotypes_predictors_folder, "/working/cleaned/homa2 calculation/homa2 indices values.xlsx")
sheets <- c("aric", "cardia", "jhs", "dppos", "mesa")
list_of_data <- lapply(sheets, function(sheet) readxl::read_excel(path_to_file, sheet = sheet))

homa2_combined <- bind_rows(list_of_data, .id = "study")
homa2_combined$study <- rep(sheets, sapply(list_of_data, nrow))

# cluster data based on 5 cohorts
final_dataset_temp = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/final_dataset_temp.RDS"))

clusters = read_csv(paste0(path_diabetes_subphenotypes_adults_folder,"/working/processed/dec_an02_clean_kmeans_5var_mi_knn_cluster.csv")) %>% 
  dplyr::select(-one_of("...1")) %>% 
  left_join(final_dataset_temp %>% 
              dplyr::select(study_id,original_study_id),
            by=c("study_id")) %>% 
  rename(cluster_study_id = study_id)


analytic_df <- pooled_df %>% 
  left_join(homa2_combined %>% 
              rename(homa2b = `HOMA2 %B`,
                     homa2ir = `HOMA2 IR`) %>% 
              dplyr::select(-c(`HOMA2 %S`,glucosef2,insulinf2)) %>% 
              distinct(study_id,study,age,.keep_all=TRUE),
            by = c("study_id","study","age")) %>% 
  left_join(clusters %>% 
              dplyr::select(cluster_study_id,original_study_id,cluster,study,female) %>% 
              mutate(study = case_when(study == "dpp" ~ "dppos",
                                       TRUE ~ study)),
            by=c("study"="study","study_id" = "original_study_id","female"))

saveRDS(analytic_df,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/8 cohorts/dsppre01_analytic df.RDS"))

write.csv(analytic_df, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/8 cohorts/dsppre01_analytic df.csv"))



### delete soon ---

impute_bmi_df <- analytic_df %>% 
  dplyr::filter(study %in% c("jhs","mesa","dppos"))

write.csv(impute_bmi_df, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/8 cohorts/dsppre01_analytic df for bmi imputation.csv"))
# bmi NA: 229
