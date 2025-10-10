rm(list = ls());gc();source(".Rprofile")

library(gridExtra)
library(ggplot2)
library(ggsurvfit)
library(survminer)
library(survival)

cluster_all_colors = c(cluster_colors_cosmos,"#CD5C5C")
names(cluster_all_colors) = c(names(cluster_colors_cosmos),"New T2D")


tdcm_coef <- read_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspse03_pooled tdcm results without dpp intervention.csv")) %>% 
  select(iv, estimate, lci, uci, model) %>% 
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
  mutate(model = factor(model,levels = names(cluster_all_colors)))


# forest plot

plot_forest <- ggplot(tdcm_coef, aes(y = term, x = estimate, xmin = lci, xmax = uci, color = model)) + 
  geom_pointrange(position = position_dodge(width = 0.7), size = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(values = cluster_all_colors) +
  scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 1)) +
  labs(
    x = "Hazard ratio (95% CI)",
    y = NULL,
    # title = "B: Hazard ratio for pathophysiological markers",
    color = "Subtype"
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  ) +
  geom_text(
    aes(x = uci + 0.01, label = HR),
    position = position_dodge(width = 0.7),
    vjust = 0.2,
    hjust = -0.05,
    fontface = "bold",
    size = 6
  ) 


ggsave(plot_forest,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/tdcm hazard ratio for pathophysiological markers without dpp intervention.png"),width=10,height=13)


