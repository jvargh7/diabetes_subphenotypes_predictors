rm(list = ls());gc();source(".Rprofile")

out_combined <- read_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan02_modeled trajectories of biomarkers.csv"))


cluster_not2d_colors = c(cluster_colors_ada,"#5C4033")
names(cluster_not2d_colors) = c(names(cluster_colors_ada),"No T2D")

fig_homa2b = out_combined %>% 
  dplyr::filter(outcome == "HOMA2B") %>% 
  mutate(cluster = factor(group,levels=c("NOT2D","MOD","SIRD","SIDD","MARD"),
                          labels=c("No T2D","MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=., aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, col = cluster)) +
  geom_path(size = 1.2) +  # thicker lines for main curves
  geom_ribbon(aes(fill = cluster), alpha = 0.1, color = NA) +  # optional soft fill for CIs
  theme_bw(base_size = 15) +
  theme(
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15)
  ) +
  xlab("Time (years)") +
  ylab("HOMA2-%B (%)") +
  scale_color_manual(name = "", values = cluster_not2d_colors) +
  scale_fill_manual(name = "", values = cluster_not2d_colors)  # if using ribbon fill


fig_homa2b

fig_bmi = out_combined %>%
  dplyr::filter(outcome == "BMI") %>% 
  mutate(cluster = factor(group,levels=c("NOT2D","MOD","SIRD","SIDD","MARD"),
                          labels=c("No T2D","MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=., aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, col = cluster)) +
  geom_path(size = 1.2) +  # thicker lines for main curves
  geom_ribbon(aes(fill = cluster), alpha = 0.1, color = NA) +  # optional soft fill for CIs
  theme_bw(base_size = 15) +
  theme(
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15)
  ) +
  xlab("Time (years)") +
  ylab("BMI (kg/m²)") +
  scale_color_manual(name = "", values = cluster_not2d_colors) +
  scale_fill_manual(name = "", values = cluster_not2d_colors) + # if using ribbon fill
  coord_cartesian(ylim = c(20, 35))

fig_hba1c = out_combined %>%
  dplyr::filter(outcome == "HbA1c") %>% 
  mutate(cluster = factor(group,levels=c("NOT2D","MOD","SIRD","SIDD","MARD"),
                          labels=c("No T2D","MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=., aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, col = cluster)) +
  geom_path(size = 1.2) +  # thicker lines for main curves
  geom_ribbon(aes(fill = cluster), alpha = 0.1, color = NA) +  # optional soft fill for CIs
  theme_bw(base_size = 15) +
  theme(
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15)
  ) +
  xlab("Time (years)") +
  ylab("HbA1c (%)") +
  scale_color_manual(name = "", values = cluster_not2d_colors) +
  scale_fill_manual(name = "", values = cluster_not2d_colors)  # if using ribbon fill

fig_homa2ir = out_combined %>%
  dplyr::filter(outcome == "HOMA2IR") %>% 
  mutate(cluster = factor(group,levels=c("NOT2D","MOD","SIRD","SIDD","MARD"),
                          labels=c("No T2D","MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=., aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, col = cluster)) +
  geom_path(size = 1.2) +  # thicker lines for main curves
  geom_ribbon(aes(fill = cluster), alpha = 0.1, color = NA) +  # optional soft fill for CIs
  theme_bw(base_size = 15) +
  theme(
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15)
  ) +
  xlab("Time (years)") +
  ylab("HOMA2-IR (%)") +
  scale_color_manual(name = "", values = cluster_not2d_colors) +
  scale_fill_manual(name = "", values = cluster_not2d_colors)  # if using ribbon fill



library(ggpubr)
ggarrange(fig_bmi,
          fig_hba1c,
          fig_homa2b,
          fig_homa2ir,
          nrow=2,ncol=2,labels=LETTERS[1:4],common.legend = TRUE,legend = "right") %>% 
  ggsave(.,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/ADA plots/trajectory of biomarkers before diagnosis.jpg"),width = 10,height = 8)
