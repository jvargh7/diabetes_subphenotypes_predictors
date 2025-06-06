rm(list = ls());gc();source(".Rprofile")


coxph_coef <- read_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan03_pooled cox ph results.csv")) %>% 
  select(iv, estimate, lci, uci, model) %>% 
  mutate(
    HR = paste0(
      formatC(round(estimate, 2), format = "f", digits = 2), " (",
      formatC(round(lci, 2), format = "f", digits = 2), ", ",
      formatC(round(uci, 2), format = "f", digits = 2), ")"
    )
  ) %>% 
  # dplyr::filter(!iv %in% c("studymesa","studyjhs","studydppos","raceNH Black","raceNH White","raceOther","female","earliest_age")) %>% 
  mutate(term = case_when(
    iv == "bmi" ~ "BMI",
    iv == "sbp_scaled" ~ "SBP",
    iv == "hba1c" ~ "HbA1c",
    iv == "ldlc_scaled" ~ "LDL",
    iv == "homa2b_scaled" ~ "HOMA2-%B",
    iv == "homa2ir" ~ "HOMA2-IR",
    iv == "egfr_ckdepi_2021_scaled" ~ "eGFR",
    TRUE ~ iv  
  )) %>% 
  mutate(model = factor(model,
                        levels = c("Overall", "MOD", "SIDD", "MARD", "SIRD"),
                        labels = c("Overall", "MOD", "SIDD", "MARD", "SIRD"))) %>% 
  select(model, term, HR) %>% 
  mutate(term = factor(term,
                       levels = c("studydppos","studyjhs","studymesa","raceNH Black","raceHispanic","raceOther","earliest_age","female",
                                  "BMI","SBP","HbA1c","LDL","HOMA2-%B","HOMA2-IR","eGFR"),
                       labels = c("Study DPPOS","Study JHS","Study MESA","NH Black","Hispanic","Other","Baseline age (years)","Female",
                                  "BMI (kg/m²)","SBP (per 10 mmHg)","HbA1c (%)","LDL (per 10 mg/dL)","HOMA2-%B (per 10%)","HOMA2-IR (%)",
                                  "eGFR (per 10 mL/min/1.73 m²)"))) %>% 
  arrange(model, term) %>% 
  write.csv(.,paste0("paper/table_cox ph time-varying coefficients.csv"))

