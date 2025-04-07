rm(list = ls());gc();source(".Rprofile")
# ,"hc","triceps","iliac","abdominal","medial" --> not there in accord
anthro_vars <- c("sbp","dbp","height","weight","bmi","wc")
# "glucose2h","insulinf","ast","apo_a","apo_b","uric_acid" --> not there in accord
lab_vars <- c("hba1c","glucosef","tgl","hdlc","ldlc","totalc","vldlc","alt",
              "serumcreatinine","urinealbumin","urinecreatinine","uacr","egfr","ratio_th")
med_vars <- c("med_chol_use","med_bp_use","med_dm_use")

# no fasting insulin; fasting glucose in mg/gl 
accord_dat_all <- readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/interim/accord_dat_all.RDS")) %>% 
  arrange(study_id,visit) %>% 
  dplyr::mutate(ratio_th=tgl/hdlc,
                bmi = weight/((height/100)^2)) %>% 
  dplyr::filter(!is.na(bsage)) %>% 
  mutate(
    # Extract the numeric part from visits that start with "F"
    visit_month = if_else(
      str_starts(visit, "F"),
      as.numeric(str_extract(visit, "[0-9.]+")),
      NA_real_,
      NA_real_
    ),
    # Add logging to identify problematic entries
    log_visit = if_else(is.na(visit_month) & !visit %in% c("BLR", "EXIT"), visit, NA_character_),
    # Calculate age at visit, handling NA values gracefully
    age = if_else(is.na(visit_month), bsage, bsage + visit_month / 12)
  )

accord_longitudinal = accord_dat_all %>% 
  arrange(study_id,visit) %>% 
  mutate(
    med_chol_use = case_when(
      grepl("Lipid Fibrate", bstreatment) ~ 1,
      TRUE ~ 0),
    med_bp_use = case_when(
      grepl("Intensive BP|Standard BP", bstreatment) ~ 1,
      TRUE ~ 0),
    med_dm_use = case_when(
      grepl("Intensive Glycemia|Standard Glycemia", bstreatment) ~ 1,
      TRUE ~ 0)) %>% 
  mutate(
    available_labs = rowSums(!is.na(.[,lab_vars])),
    available_anthro = rowSums(!is.na(.[,anthro_vars]))) %>% 
  dplyr::select(study_id,visit,bsage,age,dmagediag,dmduration,female,race_eth,alcohol,smoking,
                available_labs,available_anthro,
                one_of(anthro_vars),one_of(lab_vars),one_of(med_vars))



saveRDS(accord_longitudinal,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01g_accord.RDS"))

accord_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01g_accord.RDS"))
write_csv(accord_longitudinal,paste0(path_prediabetes_subphenotypes_folder,"/working/longitudinal/accord.csv"))
