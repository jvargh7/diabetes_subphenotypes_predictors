rm(list = ls());gc();source(".Rprofile")

# restrict to 15y of follow-up time for trajectory analysis
analytic_df = readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>% 
  dplyr::filter(t <= 0 & t >= -15) %>% 
  arrange(joint_id,t) %>% 
  distinct(joint_id,t,.keep_all=TRUE) %>% 
  mutate(across(one_of(c("subtype","race","study")),.fns=~as.factor(.)))


cluster_not2d_colors = c(cluster_colors_cosmos,"#5C4033")
names(cluster_not2d_colors) = c(names(cluster_colors_cosmos),"NOT2D")


df_long <- analytic_df %>%
  select(t, subtype, hba1c, bmi, homa2b, homa2ir) %>%
  pivot_longer(cols = c(hba1c, bmi, homa2b, homa2ir), names_to = "variable", values_to = "value") %>%
  dplyr::filter(!is.na(value)) %>%
  mutate(
    variable = recode(variable,
                      "bmi" = "BMI",
                      "hba1c" = "HbA1c",
                      "homa2b" = "HOMA2-%B",
                      "homa2ir" = "HOMA2-IR"),
    variable = factor(variable, levels = c("BMI", "HbA1c", "HOMA2-%B", "HOMA2-IR")),
    subtype = factor(subtype, levels = names(cluster_not2d_colors))  # make sure order matches color map
  ) 

df_counts <- df_long %>%
  mutate(t_bin = floor(t)) %>%
  count(t_bin, variable, subtype, name = "n")


fig_data <- ggplot(df_counts, aes(x = t_bin, y = n, fill = subtype)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ variable, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = cluster_not2d_colors) +
  labs(
    x = "Time (years)",
    y = "Number of Observations",
    fill = "Subtype"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

ggsave(fig_data,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/histogram data availability of biomarkers before diagnosis.jpg"),width = 12,height = 6)




# by visit year and study -------------------------------------------

df_hba1c_visits <- analytic_df %>%
  dplyr::filter(!is.na(hba1c)) %>%
  distinct(study, visit, visit_year_start, visit_year_end) %>%
  mutate(study = factor(study, levels = rev(sort(unique(study))),labels = c("CARDIA","DPP/OS","JHS","MESA")))  # optional: reverse for better order

# Create long format for visit points
df_hba1c_points <- df_hba1c_visits %>%
  pivot_longer(cols = c(visit_year_start, visit_year_end), names_to = "point_type", values_to = "year")

# Plot
fig_year <- ggplot() +
  geom_segment(
    data = df_hba1c_visits,
    aes(x = visit_year_start, xend = visit_year_end, y = study, yend = study),
    color = "#2C77B1", linewidth = 1
  ) +
  geom_point(
    data = df_hba1c_points,
    aes(x = year, y = study),
    color = "#D23E42", size = 2.5
  ) +
  labs(
    x = "Calendar Year",
    y = "Study"
  ) +
  scale_x_continuous(breaks = seq(1995, 2020, 5)) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(face = "bold")
  )

ggsave(fig_year,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/data availability of biomarkers before diagnosis by calendar year and study.jpg"),width = 10,height = 6)

