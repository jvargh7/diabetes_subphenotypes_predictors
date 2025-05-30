rm(list = ls());gc();source(".Rprofile")

source("functions/egfr_ckdepi_2021.R")



# new dm
final_dataset_temp = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/final_dataset_temp.RDS")) %>%
  mutate(joint_id = paste(study, original_study_id, sep = "_"))

analytic_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dsppre01_analytic df.RDS")) %>% 
  mutate(egfr_ckdepi_2021 = egfr_ckdepi_2021(scr = serumcreatinine,female = female,age = age),
         joint_id = paste(study, study_id, sep = "_")) %>%
  dplyr::filter(joint_id %in% final_dataset_temp$joint_id) %>% 
  # Restrict to newly diagnosed and those with cluster data
  dplyr::filter(!is.na(dmagediag),!is.na(cluster)) %>% 
  # Restrict to those with eGFR
  dplyr::filter(!is.na(egfr_ckdepi_2021)) %>% 
  # Restrict to values within +/- 1 year
  dplyr::filter(abs(age-dmagediag)<=1) %>% 
  group_by(joint_id) %>% 
  # Take latest value 
  dplyr::filter(age == max(age)) %>% 
  ungroup() %>% 
  distinct(joint_id,.keep_all=TRUE)


cluster_colors_cosmos

fig_egfr = analytic_df %>% 
  mutate(cluster = factor(cluster,levels=c("MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=.,aes(x=cluster,y=egfr_ckdepi_2021,fill=cluster)) +
  geom_boxplot() +
  theme_bw() + 
  xlab("") +
  ylab(bquote('eGFR 2021 ml/min/1.73' /m^2~'')) +
  scale_fill_manual(name="",values=cluster_colors_cosmos) +
  theme(legend.position = "bottom",
        legend.margin=margin(c(-3,5,5,5))) + 
  guides(fill = guide_legend(nrow=2,ncol=2))

fig_serumcreatinine = analytic_df %>% 
  mutate(cluster = factor(cluster,levels=c("MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=.,aes(x=cluster,y=serumcreatinine,fill=cluster)) +
  geom_boxplot() +
  theme_bw() + 
  xlab("") +
  ylab("Serum Creatinine (ng/dL)") +
  scale_fill_manual(name="",values=cluster_colors_cosmos) +
  guides(fill = guide_legend(nrow=2,ncol=2)) +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits=c(0,1.5),breaks=seq(0,1.5,0.3))


fig_scatter_A <- analytic_df %>% 
  mutate(cluster = factor(cluster,levels=c("MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=.,aes(col=cluster,x=dmagediag,y=egfr_ckdepi_2021)) +
  geom_density2d()+
  # geom_point() +
  theme_bw() + 
  xlab("Age at diagnosis (y)") +
  ylab(bquote('eGFR 2021 ml/min/1.73' /m^2~'')) +
  scale_color_manual(name="",values=cluster_colors_cosmos) +
  theme(legend.position = "bottom",
        legend.margin=margin(c(-3,5,5,5))) + 
  guides(fill = guide_legend(nrow=2,ncol=2))

fig_scatter_B <- analytic_df %>% 
  mutate(cluster = factor(cluster,levels=c("MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=.,aes(col=cluster,x=hba1c,y=egfr_ckdepi_2021)) +
  # geom_point() +
  geom_density2d()+
  theme_bw() + 
  xlab(bquote('BMI kg' /m^2~'')) +
  ylab(bquote('eGFR 2021 ml/min/1.73' /m^2~'')) +
  scale_color_manual(name="",values=cluster_colors_cosmos) +
  theme(legend.position = "bottom",
        legend.margin=margin(c(-3,5,5,5))) + 
  guides(fill = guide_legend(nrow=2,ncol=2))

alibrary(ggpubr)

ggarrange(fig_serumcreatinine,
          fig_egfr,
          nrow=2,ncol=1,labels=LETTERS[1:2],common.legend = TRUE,legend = "bottom") %>% 
  ggsave(.,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/creatinine and egfr.jpg"),width = 3,height = 6)



ggsave(fig_egfr,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/cross sectional egfr.jpg"),width = 3,height = 3)
