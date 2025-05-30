rm(list = ls());gc();source(".Rprofile")

library(gridExtra)
library(ggplot2)
library(ggsurvfit)
library(survminer)
library(survival)

tdcm_coef <- read_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan03_tdcm pooled results with multiple imputation.csv")) %>% 
  select(iv, estimate, lci, uci, model) %>% 
  dplyr::filter(model != "Overall") %>% 
  mutate(HR = paste0(format(round(estimate, 2), nsmall = 2), " (",
                        format(round(lci, 2), nsmall = 2), ", ",
                        format(round(uci, 2), nsmall = 2), ")")) %>% 
  dplyr::filter(!iv %in% c("studymesa","studycardia","studyjhs","studydppos","raceNH Black","raceNH White","raceOther","female1","min_age")) %>% 
  mutate(term = case_when(
    iv == "bmi" ~ "BMI",
    iv == "sbp_scaled" ~ "SBP",
    iv == "hba1c" ~ "HbA1c",
    iv == "ldlc_scaled" ~ "LDL",
    iv == "homa2b" ~ "HOMA2-%B",
    iv == "homa2ir" ~ "HOMA2-IR",
    iv == "egfr_ckdepi_2021_scaled" ~ "eGFR",
    TRUE ~ iv  
  ),
  term = factor(term,
                levels = c("eGFR", "HOMA2-IR", "HOMA2-%B", "LDL", "HbA1c", "SBP", "BMI"),
                labels = c("eGFR (per 10 mL/min/1.73 m²)", "HOMA2-IR", "HOMA2-%B", "LDL (per 10 mg/dL)", "HbA1c (%)", "SBP (per 10 mmHg)", "BMI (kg/m²)"))
  ) %>% 
  mutate(model = factor(model,
                        levels = c("MOD", "SIDD", "MARD", "SIRD"),
                        labels = c("MOD", "SIDD", "MARD", "SIRD")))


# forest plot

plot_forest <- ggplot(tdcm_coef, aes(y = term, x = estimate, xmin = lci, xmax = uci, color = model)) + 
  geom_pointrange(position = position_dodge(width = 0.7), size = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(values = cluster_colors_cosmos) +
  scale_x_continuous(limits = c(0, 3.0), breaks = seq(0, 3.0, by = 0.5)) +
  labs(
    x = "Hazard ratio (95% CI)",
    y = NULL,
    # title = "B: Hazard ratio for pathophysiological markers",
    color = "Subtype"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.4)
  ) +
  geom_text(
    aes(x = uci + 0.01, label = HR),
    position = position_dodge(width = 0.7),
    vjust = 0.2,
    hjust = -0.05,
    fontface = "bold",
    size = 5
  ) 


ggsave(plot_forest,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/hazard ratio for pathophysiological markers.png"),width=10,height=10)



#-------------------------------------------------------------------------------------------------------------------
# incidence plot

ipcw_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan02_ipcw dfs.RDS"))

tdcm_fit <- list()

for (i in 1:length(ipcw_dfs)) {
  df <- ipcw_dfs[[i]]
  
  tdcm_df <- df %>%
    arrange(study, study_id, age) %>%
    group_by(study, study_id) %>%
    mutate(
      tstart = case_when(
        row_number() == 1 ~ 0, 
        TRUE ~ age - first(age)
      ), 
      tstop = lead(tstart, default = last(age) - first(age)) 
    ) %>%
    dplyr::filter(tstart < tstop) %>% 
    mutate(
      tstart_age = case_when(row_number() == 1 ~ age, 
                         TRUE ~ dplyr::lag(age, n = 1)), 
      tstop_age = age
    ) %>% 
    dplyr::filter(tstop_age <= censored_age) %>% 
    ungroup() %>% 
    mutate(cluster = factor(cluster,
                            levels = c("MOD", "SIDD", "MARD", "SIRD"),
                            labels = c("MOD", "SIDD", "MARD", "SIRD"))) %>% 
    dplyr::filter(time_to_event <= 10)
  
  
  tdcm_fit[[i]] <- coxph(Surv(tstart, tstop, dmage_filter) ~ strata(cluster), 
                         data = tdcm_df, weights = ipcw_cluster)
  
  
  
}


survfit2(tdcm_fit[[1]])


plot_incidence = survfit2(tdcm_fit[[1]]) %>% 
  ggsurvfit(.,type = "risk") +
  xlab("Time to Diabetes (years)") +
  ylab("") +
  # ggtitle("A: Crude cumulative incidence for T2D subtypes") +
  add_confidence_interval() +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = cluster_not2d_colors, name = "Subtype") +
  scale_fill_manual(values = cluster_not2d_colors, name = "Subtype") +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.4) 
  ) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 3)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))


ggsave(plot_incidence,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/crude cumulative incidence for T2D subtypes.png"),width=12,height=8)

