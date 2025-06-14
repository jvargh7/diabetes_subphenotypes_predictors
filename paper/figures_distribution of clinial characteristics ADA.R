rm(list=ls());gc();source(".Rprofile")

library(ggplot2)

boxplot_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>%
  dplyr::filter(age == censored_age) %>% 
  dplyr::select("joint_id","subtype","hba1c", "sbp", "dbp", "ldlc", "hdlc", "bmi","homa2b","homa2ir","egfr_ckdepi_2021") %>% 
  mutate(subtype = factor(subtype,levels = c("NOT2D",'MARD','SIDD','SIRD','MOD'),labels = c("No T2D",'MARD','SIDD','SIRD','MOD')))


cluster_not2d_colors = c(cluster_colors_ada,"#5C4033")
names(cluster_not2d_colors) = c(names(cluster_colors_ada),"No T2D")

fig_A = boxplot_df %>% 
  ggplot(data=.,aes(x=subtype,y=bmi,fill=subtype)) +
  geom_boxplot(position = position_dodge(width=0.9)) +
  xlab("") +
  ylab(bquote('BMI ( kg' /m^2~')')) +
  scale_y_continuous(limits=c(10,70),breaks=seq(10,70,by=10)) +
  theme_bw() +
  scale_fill_manual(name="",values=cluster_not2d_colors) +
  theme_bw(base_size = 14) +
  theme(legend.text = element_text(size = 14)) 

fig_B = boxplot_df %>% 
  ggplot(data=.,aes(x=subtype,y=sbp,fill=subtype)) +
  geom_boxplot(position = position_dodge(width=0.9)) +
  xlab("") +
  ylab("SBP (mmHg)") +
  scale_y_continuous(limits=c(50,250),breaks=seq(50,250,by=50)) +
  theme_bw() +
  scale_fill_manual(name="",values=cluster_not2d_colors) +
  theme_bw(base_size = 14) +
  theme(legend.text = element_text(size = 14))

fig_C = boxplot_df %>% 
  ggplot(data=.,aes(x=subtype,y=hba1c,fill=subtype)) +
  geom_boxplot(position = position_dodge(width=0.9)) +
  xlab("") +
  ylab("HbA1c (%)") +
  scale_y_continuous(limits=c(3,15),breaks=seq(3,15,by=3)) +
  theme_bw() +
  scale_fill_manual(name="",values=cluster_not2d_colors) +
  theme_bw(base_size = 14) +
  theme(legend.text = element_text(size = 14)) 

fig_D = boxplot_df %>% 
  ggplot(data=.,aes(x=subtype,y=egfr_ckdepi_2021,fill=subtype)) +
  geom_boxplot(position = position_dodge(width=0.9)) +
  xlab("") +
  ylab("eGFR (per 10 mL/min/1.73 m²)") +
  scale_y_continuous(limits=c(0,180),breaks=seq(0,180,by=60)) +
  theme_bw() +
  scale_fill_manual(name="",values=cluster_not2d_colors) +
  theme_bw(base_size = 14) +
  theme(legend.text = element_text(size = 14)) 

fig_E = boxplot_df %>% 
  ggplot(data=.,aes(x=subtype,y=ldlc,fill=subtype)) +
  geom_boxplot(position = position_dodge(width=0.9)) +
  xlab("") +
  ylab("LDL Cholesterol (mg/dL)") +
  scale_y_continuous(limits=c(0,250),breaks=seq(0,250,by=50)) +
  theme_bw() +
  scale_fill_manual(name="",values=cluster_not2d_colors) +
  theme_bw(base_size = 14) +
  theme(legend.text = element_text(size = 14)) 

fig_F = boxplot_df %>% 
  ggplot(data=.,aes(x=subtype,y=hdlc,fill=subtype)) +
  geom_boxplot(position = position_dodge(width=0.9)) +
  xlab("") +
  ylab("HDL Cholesterol (mg/dL)") +
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,by=20)) +
  theme_bw() +
  scale_fill_manual(name="",values=cluster_not2d_colors) +
  theme_bw(base_size = 14) +
  theme(legend.text = element_text(size = 14)) 

fig_G = boxplot_df %>% 
  ggplot(data=.,aes(x=subtype,y=homa2b,fill=subtype)) +
  geom_boxplot(position = position_dodge(width=0.9)) +
  xlab("") +
  ylab("HOMA2-%B (%)") +
  scale_y_continuous(limits=c(0,300),breaks=seq(0,300,by=100)) +
  theme_bw() +
  scale_fill_manual(name="",values=cluster_not2d_colors) +
  theme_bw(base_size = 14) +
  theme(legend.text = element_text(size = 14)) 

fig_H = boxplot_df %>% 
  ggplot(data=.,aes(x=subtype,y=homa2ir,fill=subtype)) +
  geom_boxplot(position = position_dodge(width=0.9)) +
  xlab("") +
  ylab("HOMA2-IR (%)") +
  scale_y_continuous(limits=c(0,10),breaks=seq(0,10,by=2.5)) +
  theme_bw() +
  scale_fill_manual(name="",values=cluster_not2d_colors) +
  theme_bw(base_size = 14) +
  theme(legend.text = element_text(size = 14)) 




library(ggpubr)

ggarrange(
          fig_A,
          fig_B,
          fig_C,
          fig_D,
          fig_E,
          fig_F,
          fig_G,
          fig_H,
          nrow=2,
          ncol=4,
          common.legend = TRUE,legend = "bottom") %>% 
  ggsave(.,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/ADA plots/distribution of clinical characteristics by subtype at last follow-up.png"),width=14,height =5.5)


