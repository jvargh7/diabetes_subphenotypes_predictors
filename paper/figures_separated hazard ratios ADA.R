rm(list = ls());gc();source(".Rprofile")

library(ggplot2)
library(readr)


cluster_all_colors = c(cluster_colors_ada,"#CD5C5C")
names(cluster_all_colors) = c(names(cluster_colors_ada),"New T2D")

tdcm_coef <- read_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan03_pooled tdcm results.csv")) %>% 
  select(iv, estimate, lci, uci, model) %>% 
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
                labels = c("eGFR (per 10 mL/min/1.73 m²)", "HOMA2-IR (%)", "HOMA2-%B (per 10%)", "LDL Cholesterol (per 10 mg/dL)", "HbA1c (%)", "SBP (per 10 mmHg)", "BMI (kg/m²)"))
  ) %>% 
  mutate(model = case_when(model == "Overall" ~ "New T2D",
                           TRUE ~ model)) %>% 
  mutate(model = factor(model,levels = names(cluster_all_colors)))



# Use your original tdcm_coef
unique_terms <- unique(tdcm_coef$term)

# Create and save/display a plot for each term
for (term_i in unique_terms) {
  
  df_term <- tdcm_coef %>% dplyr::filter(term == term_i)
  
  plot_term <- ggplot(df_term, aes(y = model, x = estimate, xmin = lci, xmax = uci, color = model)) + 
    geom_pointrange(position = position_dodge(width = 0.5), size = 1) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "darkgrey") +
    scale_color_manual(values = cluster_all_colors) +
    scale_x_continuous(limits = c(0, 4), breaks = seq(0, 4, 0.5)) +
    labs(
      title = paste0(term_i),
      x = "Hazard ratio (95% CI)",
      y = NULL,
      color = "Subtype"
    ) +
    theme_bw(base_size = 16) +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15),
      panel.grid = element_blank()
    ) +
    geom_text(
      aes(x = uci + 0.05, label = HR),
      hjust = 0,
      size = 6
    )
  
  # Save the plot
  ggsave(
    plot = plot_term,
    filename = paste0(
      path_diabetes_subphenotypes_predictors_folder,
      "/figures/ADA plots/hazard ratio - ", gsub("[:/ ]", "_", as.character(term_i)), ".png"
    ),
    width = 11,
    height = 4,
    dpi = 300
  )
}

