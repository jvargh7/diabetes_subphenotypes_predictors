rm(list = ls());gc();source(".Rprofile")

# trajectory analysis ---------------------------
mean_vars = c("age","dmagediag","bmi","sbp","dbp","ldlc","hdlc",
              "insulinf","glucosef","egfr","tgl","tg_hdl")

median_vars = c("hba1c","homa2b","homa2ir")

table_df <- read_csv("analysis/dspan05_trajectory descriptive characteristics by subtype.csv") %>% 
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
  dplyr::select(variable,group,Total,MARD,MOD,NOT2D,SIDD,SIRD)

write_csv(table_df,"paper/table_trajectory descriptive characteristics by subtype.csv")  


# predictors  analysis ---------------------------
mean_vars = c("age","dmagediag","bmi","sbp","dbp","ldlc","hdlc",
              "insulinf","glucosef","egfr","tgl","tg_hdl")

median_vars = c("hba1c","homa2b","homa2ir")

table_df <- read_csv("analysis/dspan05_predictors descriptive characteristics by subtype.csv") %>% 
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
  dplyr::select(variable,group,Total,MARD,MOD,NOT2D,SIDD,SIRD)

write_csv(table_df,"paper/table_predictors descriptive characteristics by subtype.csv") 


