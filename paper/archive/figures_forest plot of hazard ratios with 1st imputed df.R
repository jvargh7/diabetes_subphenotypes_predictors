rm(list = ls());gc();source(".Rprofile")

coxph_coef <- read_csv("analysis/dspan03_cox ph with multiple imputation.csv") %>% 
  dplyr::filter(model == "m1",
                !term %in% c("raceNH Black","raceNH White","raceOther","female1","age")) %>% 
  select(-model)

library(ggplot2)
library(forestplot)
library(ggpubr)

coxph_long <- coxph_coef %>%
  gather(key = "outcome", value = "HR", Overall, MARD, MOD, SIDD, SIRD) %>%
  mutate(outcome = factor(outcome, levels = c("Overall", "MARD", "MOD", "SIDD", "SIRD"))) %>% 
  separate(HR, into = c("estimate", "ci_range"), sep = " \\(") %>%
  separate(ci_range, into = c("ci_low", "ci_high"), sep = ", ") %>%
  mutate(ci_high = sub("\\)", "", ci_high),
         estimate = as.numeric(estimate),
         ci_low = as.numeric(ci_low),
         ci_high = as.numeric(ci_high))


coxph_forest_plot <- ggplot(coxph_long, aes(x = estimate, y = term)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), 
                 position = position_dodge(width = 0.5), height = 0.2) +
  facet_grid(. ~ outcome, scales = "free_x", switch = "x") +
  theme_bw() +
  labs(x = "Hazard Ratio (95% CI)", y = "Covariates") +
  scale_x_continuous(limits = c(0, 4), breaks = 0:4) +  # Set x-axis limits and breaks
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") + # Adds a vertical line at HR = 1
  theme(
    strip.text.x = element_text(angle = 0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) 

ggsave(coxph_forest_plot,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/forest plot of cox ph hazard ratio by covariates.jpg"),width=12,height =5.5)


#--------------------------------------------------------------------------------------------------------------------

tdcm_coef <- read_csv("analysis/dspan03_tdcm with multiple imputation.csv") %>% 
  dplyr::filter(model == "m1",
                !term %in% c("raceNH Black","raceNH White","raceOther","female1","min_age")) %>% 
  select(-model)

tdcm_long <- tdcm_coef %>%
  gather(key = "Outcome", value = "HR", Overall, MARD, MOD, SIDD, SIRD) %>%
  mutate(Outcome = factor(Outcome, levels = c("Overall", "MARD", "MOD", "SIDD", "SIRD"))) %>% 
  separate(HR, into = c("estimate", "ci_range"), sep = " \\(") %>%
  separate(ci_range, into = c("ci_low", "ci_high"), sep = ", ") %>%
  mutate(ci_high = sub("\\)", "", ci_high),
         estimate = as.numeric(estimate),
         ci_low = as.numeric(ci_low),
         ci_high = as.numeric(ci_high)) %>% 
  mutate(term = case_when(
    term == "sbp" ~ "SBP",
    term == "ldlc" ~ "LDL-C",
    term == "homa2b" ~ "HOMA2-%B",
    term == "homa2ir" ~ "HOMA2-IR",
    term == "hba1c" ~ "HbA1c",
    term == "bmi" ~ "BMI",
    term == "egfr_ckdepi_2021" ~ "EGFR",
    TRUE ~ term  
  ),
  Outcome = case_when(
    Outcome == "Overall" ~ "Diabetes",
    TRUE ~ Outcome
  ))


tdcm_forest_plot <- ggplot(tdcm_long, aes(x = estimate, y = term, color = Outcome, shape = Outcome)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), 
                 position = position_dodge(width = 0.5), height = 0.2) +
  facet_grid(. ~ Outcome, scales = "free_x", switch = "x") +
  theme_bw() +
  labs(x = "Hazard Ratio (95% CI)", y = "Biomarkers") +
  scale_x_continuous(limits = c(0, 2), breaks = 0:2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange")) + # Customize colors for 5 outcomes
  scale_shape_manual(values = c(16, 17, 18, 19, 15)) + # Customize shapes for 5 outcomes
  theme(
    strip.text.x = element_text(angle = 0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "bottom"  # Adjust legend position if needed
  )

ggsave(tdcm_forest_plot,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/forest plot of tdcm hazard ratio by covariates.jpg"),width=12,height =5.5)
