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
data <- data.frame(
  Study = rep(c("CARDIA", "DPP/OS", "JHS", "MESA"), each = 21),
  Year = rep(2000:2020, times = 4),
  HbA1c = c(
    # CARDIA: collected at 2000, 2005, 2010
    rep(0, 0), rep(1, 1), rep(0, 4), rep(1, 1), rep(0, 4), rep(1, 1), rep(0, 10),
    
    # DPP/OS: annual measurement
    rep(1, 21),
    
    # JHS: Visit 1 (2000–2004), Visit 2 (2005–2008), Visit 3 (2009–2013)
    rep(1, 5), rep(1, 4), rep(1, 5), rep(0, 7),
    
    # MESA: Visit 1 (2000–2002), Visit 5 (2010–2012), Visit 6 (2016–2018)
    rep(1, 3), rep(0, 7), rep(1, 3), rep(0, 3), rep(1, 3), rep(0, 2)
  )
)

data$Study <- factor(data$Study, levels = c("MESA",'JHS','DPP/OS','CARDIA'))

fig_year <- ggplot(data, aes(x = Year, y = Study)) +
  geom_point(aes(shape = factor(HbA1c), color = factor(HbA1c), fill = factor(HbA1c)), size = 3, stroke = 1.2) +
  scale_shape_manual(values = c(1, 16), labels = c("Not Collected", "Collected")) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("white", "red")) +
  labs(x = "Calendar Year", y = NULL, shape = "HbA1c Status", color = NULL) +
  theme_minimal(base_size = 15) +
  theme(
    plot.background = element_rect(color = "black", fill = NA, size = 1),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(1985, 2020), breaks = seq(1985, 2020, by = 5))


ggsave(fig_year,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/data availability of hba1c by calendar year and study.jpg"),width = 10,height = 6)

