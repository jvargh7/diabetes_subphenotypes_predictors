rm(list = ls());gc();source(".Rprofile")


# predictors  analysis ---------------------------
mean_vars = c("age","dmagediag","bmi","sbp","dbp","ldlc","hdlc",
              "insulinf","glucosef","egfr_ckdepi_2021","tgl","tg_hdl")

median_vars = c("hba1c","homa2b","homa2ir","fu_time","visit_number")

# Baseline
table_df_baseline <- read_csv("analysis/dspan04_predictors descriptive characteristics by subtype at baseline.csv") %>% 
  dplyr::mutate(selected_rows = case_when(variable %in% mean_vars & est %in% c("mean","sd") ~ 1,
                                          variable %in% median_vars & est %in% c("median","q25","q75") ~ 1,
                                          !variable %in% c(mean_vars,median_vars) ~ 1,
                                          TRUE ~ 0
  )) %>% 
  dplyr::filter(selected_rows == 1) %>% 
  dplyr::select(subtype,group,variable,est,value) %>% 
  pivot_wider(names_from=est,values_from=value) %>% 
  mutate(output = case_when(variable %in% mean_vars ~ paste0(round(mean,1)," (",round(sd,1),")"),
                            variable %in% median_vars ~ paste0(round(median,1)," (",round(q25,1),", ",round(q75,1),")"),
                            TRUE ~ paste0(round(freq,0)," (",round(proportion,1),"%)")
  )) %>% 
  dplyr::select(variable,group,subtype,output) %>% 
  pivot_wider(names_from=subtype,values_from=output) %>% 
  dplyr::select(variable,group,Total,NOT2D,MARD,MOD,SIDD,SIRD) %>% 
  write_csv(.,"paper/table_predictors descriptive characteristics by subtype at baseline.csv") 


# last FU
table_df_lastfu <- read_csv("analysis/dspan04_predictors descriptive characteristics by subtype at last follow-up.csv") %>% 
  dplyr::mutate(selected_rows = case_when(variable %in% mean_vars & est %in% c("mean","sd") ~ 1,
                                          variable %in% median_vars & est %in% c("median","q25","q75") ~ 1,
                                          !variable %in% c(mean_vars,median_vars) ~ 1,
                                          TRUE ~ 0
  )) %>% 
  dplyr::filter(selected_rows == 1) %>% 
  dplyr::select(subtype,group,variable,est,value) %>% 
  pivot_wider(names_from=est,values_from=value) %>% 
  mutate(output = case_when(variable %in% mean_vars ~ paste0(round(mean,1)," (",round(sd,1),")"),
                            variable %in% median_vars ~ paste0(round(median,1)," (",round(q25,1),", ",round(q75,1),")"),
                            TRUE ~ paste0(round(freq,0)," (",round(proportion,1),"%)")
  )) %>% 
  dplyr::select(variable,group,subtype,output) %>% 
  pivot_wider(names_from=subtype,values_from=output) %>% 
  dplyr::select(variable,group,Total,NOT2D,MARD,MOD,SIDD,SIRD) %>% 
  write_csv(.,"paper/table_predictors descriptive characteristics by subtype at last follow-up.csv")

