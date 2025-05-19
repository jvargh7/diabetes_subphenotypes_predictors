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
            by=c("study"="study","study_id" = "original_study_id","female")) %>% 
  # exclude outliers
  mutate(
    weight = case_when(weight <= 150 ~ weight,
                    TRUE ~ NA_real_),
    bmi = case_when(bmi >= 15 & bmi <= 65 ~ bmi,
                    TRUE ~ NA_real_),
    wc = case_when(wc >= 60 & wc <= 185 ~ bmi,
                    TRUE ~ NA_real_),
    sbp = case_when(sbp >= 90 & sbp <= 200 ~ sbp,
                    TRUE ~ NA_real_),
    dbp = case_when(dbp >= 60 & dbp <= 120 ~ dbp,
                    TRUE ~ NA_real_),
    hba1c = case_when(hba1c <= 14 ~ hba1c,
                       TRUE ~ NA_real_),
    totalc = case_when(totalc >= 100 & totalc <= 300 ~ totalc,
                    TRUE ~ NA_real_),
    ldlc = case_when(ldlc >= 50 & ldlc <= 200 ~ ldlc,
                       TRUE ~ NA_real_),
    hdlc = case_when(hdlc >= 20 & hdlc <= 100 ~ hdlc,
                     TRUE ~ NA_real_),
    vldlc = case_when(vldlc >= 5 & vldlc <= 40 ~ vldlc,
                     TRUE ~ NA_real_),
    glucosef = case_when(glucosef >= 50 & glucosef <= 300 ~ glucosef,
                      TRUE ~ NA_real_),
    insulinf = case_when(insulinf >= 2 & insulinf <= 50 ~ insulinf,
                         TRUE ~ NA_real_),
    glucose2h = case_when(glucose2h >= 70 & glucose2h <= 300 ~ glucose2h,
                         TRUE ~ NA_real_),
    tgl = case_when(tgl >= 50 & tgl <= 600 ~ tgl,
                          TRUE ~ NA_real_),
    serumcreatinine = case_when(serumcreatinine >= 0.4 & serumcreatinine <= 2 ~ serumcreatinine,
                    TRUE ~ NA_real_),
    urinecreatinine = case_when(urinecreatinine >= 20 & urinecreatinine <= 400 ~ urinecreatinine,
                                TRUE ~ NA_real_),
    urinealbumin = case_when(urinealbumin <= 100 ~ urinealbumin,
                                TRUE ~ NA_real_),
    uacr = case_when(uacr <= 16.6 ~ uacr,
                             TRUE ~ NA_real_),
    egfr = case_when(egfr >= 15 & egfr <= 120 ~ egfr,
                     TRUE ~ NA_real_),
    homa2b = case_when(homa2b >= 80 & homa2b <= 120 ~ homa2b,
                     TRUE ~ NA_real_),
    homa2ir = case_when(homa2ir >= 0.5 & homa2ir <= 2.72 ~ homa2ir,
                     TRUE ~ NA_real_)
  )



saveRDS(analytic_df,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/8 cohorts/dsppre01_analytic df.RDS"))

write.csv(analytic_df, paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/8 cohorts/dsppre01_analytic df.csv"))


# check outliers--------------------------------------
continuous_vars <- c("age", "height","weight","bmi","wc","sbp", "dbp","hba1c", 
                     "totalc","ldlc","hdlc","vldlc","glucosef","insulinf","glucose2h",
                     "tgl", "serumcreatinine","urinecreatinine","urinealbumin",
                     "uacr","egfr", "homa2b", "homa2ir")

summary_stats <- analytic_df %>%
  summarise(across(all_of(continuous_vars), list(
    min = ~min(., na.rm = TRUE),
    qua1st = ~quantile(., 0.25, na.rm = TRUE),
    median = ~median(., na.rm = TRUE),
    mean = ~mean(., na.rm = TRUE),
    qua3rd = ~quantile(., 0.75, na.rm = TRUE),
    max = ~max(., na.rm = TRUE),
    NAs = ~sum(is.na(.))
  ), .names = "{.col}_{.fn}")) %>%
  pivot_longer(cols = everything(), names_to = "name", values_to = "value") %>%
  separate(name, into = c("variable", "stat"), sep = "_") %>%
  pivot_wider(names_from = stat, values_from = value) %>% 
  mutate(iqr = qua3rd - qua1st,
         lower_bound = qua1st - 1.5*iqr,
         upper_bound = qua3rd + 1.5*iqr)


write.csv(summary_stats, "analysis/check outliers before mice.csv")

