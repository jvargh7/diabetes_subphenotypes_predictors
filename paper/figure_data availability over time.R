rm(list = ls());gc();source(".Rprofile")

# restrict to 15y of follow-up time for trajectory analysis
analytic_df = readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>% 
  dplyr::filter(t <= 0 & t >= -15) %>% 
  arrange(joint_id,t) %>% 
  distinct(joint_id,t,.keep_all=TRUE) %>% 
  mutate(across(one_of(c("subtype","race","study")),.fns=~as.factor(.))) %>% 
  dplyr::filter(subtype != "NOT2D")



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
    subtype = factor(subtype, levels = names(cluster_colors_ada))  # make sure order matches color map
  ) 

df_counts <- df_long %>%
  mutate(t_bin = floor(t)) %>%
  count(t_bin, variable, subtype, name = "n")


fig_data <- ggplot(df_counts, aes(x = t_bin, y = n, fill = subtype)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ variable, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = cluster_colors_ada) +
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
library(ggplot2)
library(xlsx)

a1c_df <- read.xlsx("data/table_hba1c availability by study over time.xlsx") %>%
  mutate(study = factor(study, levels = c("mesa","jhs","dppos","cardia"),
                        labels = c("MESA","JHS","DPP/OS","CARDIA"))) %>% 
  rowwise() %>%
  mutate(year = list(seq(visit_start, visit_stop))) %>%
  unnest(cols = c(year))  # this ensures `year` becomes a flat column

fig_year <- ggplot() +
  geom_segment(data = a1c_df %>% distinct(study, visit, visit_start, visit_stop),
               aes(x = visit_start, xend = visit_stop,
                   y = study, yend = study),
               color = "blue", size = 0.7) +
  geom_point(data = a1c_df,
             aes(x = year, y = study,
                 shape = factor(a1c_ava),
                 fill = factor(a1c_ava)),
             color = "red", size = 2) +
  scale_shape_manual(values = c(1, 21), labels = c("No HbA1c", "Has HbA1c")) +
  scale_fill_manual(values = c(NA, "red"), guide = "none") +
  labs(x = "Calendar Year", y = "Study", shape = "HbA1c Availability") +
  scale_x_continuous(breaks = seq(1980, 2025, 5), limits = c(1980, 2025)) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

ggsave(fig_year,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/data availability of hba1c by calendar year and study.jpg"),width = 10,height = 6)

