rm(list = ls());gc();source(".Rprofile")

library(gridExtra)
library(ggplot2)
library(ggsurvfit)
library(survminer)
library(survival)

cluster_all_colors = c(cluster_colors_cosmos,"#CD5C5C")
names(cluster_all_colors) = c(names(cluster_colors_cosmos),"Overall")


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
  dplyr::filter(!iv %in% c("studymesa","studyjhs","studydppos","raceNH Black","raceNH White","raceOther","female","earliest_age","dpp_intervention")) %>% 
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
                levels = c("eGFR", "HOMA2-IR", "HOMA2-%B", "LDL", "HbA1c", "SBP", "BMI"),
                labels = c("eGFR (per 10 mL/min/1.73 m²)", "HOMA2-IR (%)", "HOMA2-%B (per 10%)", "LDL (per 10 mg/dL)", "HbA1c (%)", "SBP (per 10 mmHg)", "BMI (kg/m²)"))
  ) %>% 
  mutate(model = factor(model,levels = names(cluster_all_colors)))


# forest plot

plot_forest <- ggplot(tdcm_coef, aes(y = term, x = estimate, xmin = lci, xmax = uci, color = model)) + 
  geom_pointrange(position = position_dodge(width = 0.7), size = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(values = cluster_all_colors) +
  scale_x_continuous(limits = c(0, 5.5), breaks = seq(0, 5.5, by = 0.5)) +
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


ggsave(plot_forest,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/hazard ratio for pathophysiological markers.png"),width=10,height=11.5)


