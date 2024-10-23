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
  mutate(study = "jhs")
dppos_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01c_dppos.RDS")) %>% 
  mutate(study = "dppos")
mesa_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01f_mesa.RDS")) %>% 
  mutate(study = "mesa")


longitudinal_df = bind_rows(aric_longitudinal %>% mutate(study_id = as.numeric(str_replace(study_id,"C",""))),
                            cardia_longitudinal,
                            jhs_longitudinal,
                            dppos_longitudinal,
                            mesa_longitudinal) %>%
  
  bind_rows(final_dataset_temp %>% dplyr::select(-study_id) %>% rename(study_id = original_study_id) %>% 
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
list_of_data <- lapply(sheets, function(sheet) read_excel(path_to_file, sheet = sheet))

homa2_combined <- bind_rows(list_of_data, .id = "study")
homa2_combined$study <- rep(sheets, sapply(list_of_data, nrow))

analytic_df <- longitudinal_df %>% 
  left_join(homa2_combined %>% 
              rename(homa2b = `HOMA2 %B`,
                     homa2ir = `HOMA2 IR`) %>% 
              dplyr::select(-c(`HOMA2 %S`,glucosef2,insulinf2)),
            by = c("study_id","study","age"))

#---------------------------------------------------------------------------------------------------------------------------------------



# DM cases, 5y before diagnosis records, N = 1129
dm_df <- analytic_df %>% 
  dplyr::filter(!is.na(dmagediag)) %>% 
  mutate(
    time_diff = dmagediag - age,  
    time_diff5 = abs(time_diff - 5) 
  ) %>%
  group_by(study_id) %>%
  dplyr::filter(time_diff > 0 & time_diff5 <= 1) %>%  
  slice_min(order_by = time_diff5, with_ties = FALSE) %>%  
  ungroup() %>% 
  dplyr::filter(!is.na(bmi) & !is.na(hba1c)) %>% 
  dplyr::select(-c("time_diff","time_diff5"))


# Non-DM cases
ndm_df <- analytic_df %>% 
  dplyr::filter(is.na(dmagediag)) %>% 
  group_by(study_id) %>%
  mutate(age_range = max(age) - min(age)) %>%  
  dplyr::filter(age_range >= 4) %>% 
  ungroup() %>%
  dplyr::select(-age_range)

# pair with previous year, N = 8244
ndm_df_joined <- ndm_df %>% 
  rename(current_age = age) %>% 
  mutate(current_age_minus5 = current_age - 5) %>% 
  left_join(ndm_df %>% 
              dplyr::select(study_id,age) %>% 
              rename(previous_age = age), by = "study_id") %>%
  dplyr::filter(previous_age < current_age) %>%
  mutate(age_diff = abs(previous_age - current_age_minus5)) %>%
  dplyr::filter(age_diff <= 1) %>% 
  group_by(study_id, current_age) %>%
  slice(1) %>% 
  ungroup() %>% 
  dplyr::filter(!is.na(bmi) & !is.na(hba1c)) %>% 
  dplyr::select(-c("current_age", "current_age_minus5", "age_diff"))

#-------------------------------------------------------------------------------------------------------------------------------
clusters = read_csv(paste0(path_diabetes_subphenotypes_adults_folder,"/working/processed/dec_an02_clean_kmeans_5var_mi_knn_cluster.csv")) %>% 
  dplyr::select(-one_of("...1")) %>% 
  left_join(final_dataset_temp %>% 
              dplyr::select(study_id,original_study_id),
            by=c("study_id")) %>% 
  rename(cluster_study_id = study_id)

# match DM and non-DM cases based on age, bmi, hba1c, study
full_df <- bind_rows(dm_df %>% 
                         mutate(dm = 1), 
                       ndm_df_joined %>% 
                         mutate(dm = 0) %>% 
                       rename(age = previous_age)) %>% 
  left_join(clusters %>% 
              dplyr::select(cluster_study_id,original_study_id,cluster,study,female),
            by=c("study"="study","study_id" = "original_study_id")) %>% 
  
  dplyr::filter(!is.na(cluster_study_id)) %>% 
  mutate(t = age - dmagediag) %>% 
  dplyr::filter(t >= -15, t < 0)

# add cluster varfiable


library(MatchIt)

match_model <- matchit(dm ~ age + bmi + hba1c + study, data = full_df,
                       method = "nearest", ratio = 1)  # 'ratio = 1' for 1:1 matching

matched_df<- match.data(match_model)

saveRDS(matched_df, paste0(path_diabetes_subphenotypes_predictors_folder, "/working/cleaned/dsp01_matched df.RDS"))

glm(dm ~ homa2ir,data=matched_df,family=binomial()) %>% summary()
glm(dm ~ hba1c,data=matched_df,family=binomial()) %>% summary()
glm(dm ~ cluster,data=matched_df,family=binomial()) %>% summary()

ggplot(data=matched_df,aes(x=factor(dm),y=bmi)) +
  geom_boxplot()



prediction model, vars typically in EHR (no fasting stuff, HOMA2; can use slef-reported data)
# multinomial regression, newly diagnosed 3377
# pooled logitic regression, ask shivani
# discrete time event analysis

2y risk of subphenotypes
at this age, what is the probability of falling a subphenotype in how many/2 years from then on
any age, 2y window, the prob of falling into one of the subpenotype
