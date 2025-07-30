rm(list = ls());gc();source(".Rprofile")


lastfu_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>%
  group_by(joint_id) %>% 
  mutate(fu_time = first(censored_age) - first(earliest_age)) %>% 
  ungroup() %>% 
  mutate(
    event_type = case_when(
      subtype == "NOT2D" ~ "censored",
      subtype == "MARD"  ~ "MARD",
      subtype == "MOD"   ~ "MOD",
      subtype == "SIDD"  ~ "SIDD",
      subtype == "SIRD"  ~ "SIRD",
      TRUE ~ NA_character_
    ),
    event_type = factor(event_type, levels = c("censored", "MARD", "MOD", "SIDD", "SIRD")),
    event_status = as.numeric(event_type) - 1 
  ) %>% 
  # only last follow-up
  dplyr::filter(age == censored_age) %>% 
  select(joint_id, fu_time, event_status, event_type)


library(prodlim)

fit <- prodlim(Hist(fu_time, event_status, cens.code = 0) ~ 1,
               data = lastfu_df)


# Get the CIF summary for all causes
fit_summary <- summary(fit, times = seq(0, max(lastfu_df$fu_time), by = 0.5))

cif_df <- as_tibble(fit_summary) %>%
  dplyr::filter(cause %in% 1:4) %>%  # Keep only subtypes
  mutate(
    cause = factor(cause,
                   levels = c(1, 2, 3, 4),
                   labels = c("MARD", "MOD", "SIDD", "SIRD"))
  )


cif_plot <- ggplot(cif_df, aes(x = time, y = cuminc, color = cause)) +
  geom_step(size = 1) +
  labs(
    x = "Time from first included visit (years)",
    y = "Cumulative incidence",
    color = "Subtype"
  ) +
  scale_color_manual(values = cluster_colors_ada) +
  scale_y_continuous(limits = c(0, 0.50)) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    text = element_text(size = 16)
  )


ggsave(
  filename = paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/Aalen-Johansen cumulative incidence curves.png"), 
  plot = cif_plot,
  width = 6,      
  height = 6,      
  dpi = 300      
)




