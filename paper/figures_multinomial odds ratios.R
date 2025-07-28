rm(list = ls());gc();source(".Rprofile")

library(gridExtra)
library(ggplot2)
library(ggsurvfit)
library(survminer)
library(survival)

cluster_all_colors = c(cluster_colors_cosmos,"#CD5C5C")
names(cluster_all_colors) = c(names(cluster_colors_cosmos),"New T2D")


multinom_coef <- read_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan05_pooled multinomial results.csv")) %>% 
  select(iv, estimate, lci, uci, `y.level`, reference) %>% 
  rename(model = `y.level`) %>% 
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
  mutate(model = factor(model,levels = names(cluster_all_colors)))


# forest plot

multinom_all <- multinom_coef %>% dplyr::filter(reference == "NOT2D")

plot_forest <- ggplot(multinom_all, aes(y = term, x = estimate, xmin = lci, xmax = uci, color = model)) + 
  geom_pointrange(position = position_dodge(width = 0.7), size = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(values = cluster_all_colors) +
  scale_x_continuous(limits = c(0, 11), breaks = seq(0, 11, by = 1)) +
  labs(
    x = "Odds ratio (95% CI)",
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


multinom_ref <- multinom_coef %>% dplyr::filter(reference == "MARD")

plot_forest_ref <- ggplot(multinom_ref, aes(y = term, x = estimate, xmin = lci, xmax = uci, color = model)) + 
  geom_pointrange(position = position_dodge(width = 0.7), size = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(values = cluster_all_colors) +
  scale_x_continuous(limits = c(0, 11), breaks = seq(0, 11, by = 1)) +
  labs(
    x = "Odds ratio (95% CI)",
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


library(patchwork)

plot_forest <- plot_forest + ggtitle("A. Reference group = 'NOT2D'")
plot_forest_ref <- plot_forest_ref + 
  ggtitle("B. Reference group = 'MARD'") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

combined_plot <- plot_forest + plot_forest_ref +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(combined_plot,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/multinomial odds ratio for pathophysiological markers.png"),width=17,height=13)


