rm(list = ls());gc();source(".Rprofile")

library(gridExtra)
library(ggplot2)
library(ggsurvfit)
library(survminer)
library(survival)

cluster_all_colors = c(cluster_colors_cosmos,"#CD5C5C")
names(cluster_all_colors) = c(names(cluster_colors_cosmos),"New T2D")


tdcm_coef <- read_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan03_pooled tdcm results.csv")) %>% 
  select(iv, estimate, lci, uci, model) %>% 
  # dplyr::filter(model != "Overall") %>% 
  # mutate(HR = paste0(format(round(estimate, 2), nsmall = 2), " (",
  #                       format(round(lci, 2), nsmall = 2), ", ",
  #                       format(round(uci, 2), nsmall = 2), ")")) %>% 
  mutate(
    HR = paste0(
      formatC(round(estimate, 2), format = "f", digits = 2), " (",
      formatC(round(lci, 2), format = "f", digits = 2), ", ",
      formatC(round(uci, 2), format = "f", digits = 2), ")"
    )
  ) %>%
  dplyr::filter(iv %in% c("bmi","egfr_ckdepi_2021_scaled","hba1c","homa2b_scaled","homa2ir","ldlc_scaled","sbp_scaled")) %>% 
  mutate(term = case_when(
    iv == "bmi" ~ "BMI",
    iv == "sbp_scaled" ~ "SBP",
    iv == "hba1c" ~ "HbA1c",
    iv == "ldlc_scaled" ~ "LDL",
    iv == "homa2b_scaled" ~ "HOMA2-%B",
    iv == "homa2ir" ~ "HOMA2-IR",
    iv == "egfr_ckdepi_2021_scaled" ~ "eGFR",
    TRUE ~ iv  
  ),
  term = factor(term,
                levels = c("eGFR", "LDL", "SBP", "HOMA2-IR", "HOMA2-%B", "HbA1c", "BMI"),
                labels = c("eGFR (per 10 mL/\nmin/1.73 m²)", "LDL \n(per 10 mg/dL)", "SBP \n(per 10 mmHg)", 
                           "HOMA2-IR", "HOMA2-%B \n(per 10 units)", "HbA1c (%)", "BMI (kg/m²)"))
  ) %>% 
  mutate(model = case_when(model == "Overall" ~ "New T2D",
                           TRUE ~ model)) %>% 
  mutate(model = factor(model,levels = names(cluster_all_colors))) %>% 
  dplyr::filter(model == "New T2D") 


# forest plot

plot_forest <- ggplot(tdcm_coef, aes(y = term, x = estimate, xmin = lci, xmax = uci)) + 
  geom_pointrange(color = cluster_all_colors["New T2D"], size = 1.1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgrey") +
  scale_x_continuous(limits = c(0, 2), breaks = seq(0, 2, by = 0.5)) +
  labs(
    x = "Hazard ratio (95% CI)",
    y = NULL
  ) +
  theme_bw(base_size = 20) +  # Increased font size
  theme(
    legend.position = "none",  # Remove legend since only one model
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.6)
  )


ggsave(plot_forest,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/tdcm hazard ratio for pathophysiological markers newT2D only.png"),width=7,height=8)

