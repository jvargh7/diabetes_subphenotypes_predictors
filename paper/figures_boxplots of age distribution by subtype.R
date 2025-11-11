rm(list = ls()); gc(); source(".Rprofile")

library(ggplot2)
library(patchwork)

# Load data
df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder, "/working/processed/dspan01_analytic sample.RDS"))

# Filter data by subtype
df_not2d <- df %>% dplyr::filter(subtype == "NOT2D")
df_mard <- df %>% dplyr::filter(subtype == "MARD")
df_mod <- df %>% dplyr::filter(subtype == "MOD")
df_sidd <- df %>% dplyr::filter(subtype == "SIDD")
df_sird <- df %>% dplyr::filter(subtype == "SIRD")

# Plot 1: NOT2D
p_not2d <- ggplot(df_not2d, aes(x = earliest_age, y = censored_age)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  labs(
    title = "NOT2D",
    x = NULL,
    y = "Age at Last Visit/T2D Diagnosis"
  ) +
  xlim(20, 100) +
  ylim(20, 100) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold",),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Plot 2: MARD
p_mard <- ggplot(df_mard, aes(x = earliest_age, y = censored_age)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  labs(
    title = "MARD",
    x = NULL,
    y = NULL
  ) +
  xlim(20, 100) +
  ylim(20, 100) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Plot 3: MOD
p_mod <- ggplot(df_mod, aes(x = earliest_age, y = censored_age)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  labs(
    title = "MOD",
    x = NULL,
    y = NULL
  ) +
  xlim(20, 100) +
  ylim(20, 100) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Plot 4: SIDD
p_sidd <- ggplot(df_sidd, aes(x = earliest_age, y = censored_age)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  labs(
    title = "SIDD",
    x = NULL,
    y = NULL
  ) +
  xlim(20, 100) +
  ylim(20, 100) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Plot 5: SIRD
p_sird <- ggplot(df_sird, aes(x = earliest_age, y = censored_age)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  labs(
    title = "SIRD",
    x = NULL,
    y = NULL
  ) +
  xlim(20, 100) +
  ylim(20, 100) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white", color = NA)
  )


# Combine all 5 plots into a single 5-column figure
combined_plot <- (p_not2d | p_mard | p_mod | p_sidd | p_sird) +
  plot_layout(guides = "collect") +
  plot_annotation(caption = "Age at First Selected Visit") &
  theme(plot.caption = element_text(hjust = 0.5, size = 11, face = "bold", margin = margin(t = 10)))


ggsave(filename = paste0(path_diabetes_subphenotypes_predictors_folder, "/figures/scatterplot of age by subtype combined.png"), 
       combined_plot, width = 13, height = 5.5, dpi = 300)

