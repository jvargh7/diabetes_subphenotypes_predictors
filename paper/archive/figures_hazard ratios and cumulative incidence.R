rm(list = ls());gc();source(".Rprofile")

library(gridExtra)
library(ggplot2)
library(survminer)

tdcm_coef <- read_csv("analysis/dspan03_tdcm pooled results with multiple imputation.csv") %>% 
  select(iv, estimate, lci, uci, model) %>% 
  # mutate(HR = paste0(round(estimate, 2), " (", round(lci, 2), ", ", round(uci, 2), ")")) %>% 
  mutate(HR = paste0(format(round(estimate, 2), nsmall = 2), " (",
                        format(round(lci, 2), nsmall = 2), ", ",
                        format(round(uci, 2), nsmall = 2), ")")) %>% 
  dplyr::filter(!iv %in% c("studymesa","studyjhs","race_cleanNH Black","race_cleanNH White","race_cleanOther","female1","min_age"),
                model != "Overall") %>% 
  mutate(term = case_when(
    iv == "bmi" ~ "BMI",
    iv == "sbp" ~ "SBP",
    iv == "hba1c" ~ "HbA1c",
    iv == "ldlc" ~ "LDL",
    iv == "homa2b" ~ "HOMA2-%B",
    iv == "homa2ir" ~ "HOMA2-IR",
    iv == "egfr_ckdepi_2021" ~ "eGFR",
    TRUE ~ iv  
  ),
  term = factor(term,
                levels = c("eGFR", "HOMA2-IR", "HOMA2-%B", "LDL", "HbA1c", "SBP", "BMI"),
                labels = c("eGFR", "HOMA2-IR", "HOMA2-%B", "LDL", "HbA1c", "SBP", "BMI"))
  )


# forest plot

plot_forest <- ggplot(tdcm_coef, aes(y = term, x = estimate, xmin = lci, xmax = uci, color = model)) + 
  geom_pointrange(position = position_dodge(width = 0.7), size = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(values = cluster_colors) +
  scale_x_continuous(limits = c(0, 3.3)) +
  labs(
    x = "Hazard ratio (95% CI)",
    y = NULL,
    title = "B: Hazard ratio of pathophysiological markers",
    color = "Subtype"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.4)
  ) +
  geom_text(
    aes(x = uci + 0.01, label = HR),
    position = position_dodge(width = 0.7),
    vjust = 0.2,
    hjust = -0.05,
    fontface = "bold",
    size = 4
  ) 



#-------------------------------------------------------------------------------------------------------------------
# incidence plot

ipcw_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/8 cohorts/dspan02_ipcw dfs.RDS"))

tdcm_fit <- list()

for (i in 1:length(ipcw_dfs)) {
  df <- ipcw_dfs[[i]]
  
  tdcm_df <- df %>%
    arrange(study, study_id, age) %>%
    group_by(study, study_id) %>%
    mutate(
      # tstart = case_when(row_number() == 1 ~ age, 
      #                    TRUE ~ dplyr::lag(age, n = 1)), 
      # tstop = age
      baseline_age = first(age),  # Assuming 'age' at first observation is baseline
      tstart = case_when(
        row_number() == 1 ~ 0, 
        TRUE ~ age - first(age)
      ), 
      tstop = lead(tstart, default = last(age) - first(age)) 
    ) %>%
    ungroup() %>% 
    dplyr::filter((tstart < tstop) & (tstop <= censored_age)) %>% 
    mutate(cluster = factor(cluster,
                            levels = c("MOD", "SIDD", "MARD", "SIRD"),
                            labels = c("MOD", "SIDD", "MARD", "SIRD"))) %>% 
    dplyr::filter(time_to_event <= 10)
  
  
  tdcm_fit[[i]] <- coxph(Surv(tstart, tstop, event) ~ strata(cluster), 
                         data = tdcm_df, weights = ipcw_cluster)
  
  
  
}


survfit2(tdcm_fit[[1]])

plot_incidence = survfit2(tdcm_fit[[1]]) %>% 
  ggsurvfit(.,type = "risk") +
  xlab("Time to Diabetes (years)") +
  ylab("") +
  ggtitle("A: Crude cumulative incidence for T2D subtypes") +
  add_confidence_interval() +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = cluster_colors, name = "Subtype") +
  scale_fill_manual(values = cluster_colors, name = "Subtype") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.4) 
  ) +
  scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))




final_plot <- grid.arrange(plot_incidence, plot_forest, ncol = 2)

ggsave(final_plot,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/hazard ratios and incidence by subtype.png"),width=16,height=8)



