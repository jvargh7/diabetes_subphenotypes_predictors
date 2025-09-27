rm(list = ls());gc();source(".Rprofile")

# restrict to 15y of follow-up time for trajectory analysis
analytic_df = readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan01_analytic sample.RDS")) %>% 
  dplyr::filter(t <= 0 & t >= -15) %>% 
  arrange(joint_id,t) %>% 
  distinct(joint_id,t,.keep_all=TRUE) %>% 
  mutate(across(one_of(c("subtype","race","study")),.fns=~as.factor(.))) %>% 
  mutate(dpp_intervention = case_when(dpp_intervention == 1 ~ "intervention arm",
                                      TRUE ~ "others"))


summary(analytic_df[,c("homa2b","homa2ir","t","subtype","censored_age","female","study","race")])
table(analytic_df$subtype,useNA="always")


analytic_df %>% 
  group_by(subtype) %>% 
  distinct(study,study_id,joint_id) %>% 
  tally() %>% 
  mutate(p = n/sum(n))

ind_df = analytic_df %>% 
  group_by(study,study_id,joint_id) %>% 
  dplyr::filter(t == max(t)) %>% 
  ungroup()

ind_df %>% 
  group_by(subtype == "NOT2D") %>% 
  summarize(age = mean(censored_age),
            sd_age = sd(censored_age),
            female = mean(female))

library(splines)
library(emmeans)
library(ggeffects)
emm_options(pbkrtest.limit = 999999)


library(lme4)

m1 = lmer(ldlc ~ subtype*ns(t,df = 3) + censored_age + female + study + race  + dpp_intervention + (1|joint_id), data = analytic_df)
m2 = lmer(hdlc ~ subtype*ns(t,df = 3) + censored_age + female + study + race  + dpp_intervention + (1|joint_id), data = analytic_df)
m3 = lmer(tgl ~ subtype*ns(t,df = 3) + censored_age + female + study + race  + dpp_intervention + (1|joint_id), data = analytic_df)
m4 = lmer(sbp ~ subtype*ns(t,df = 3) + censored_age + female + study + race  + dpp_intervention + (1|joint_id), data = analytic_df)
m5 = lmer(egfr_ckdepi_2021 ~ subtype*ns(t,df = 3) + censored_age + female + study + race  + dpp_intervention + (1|joint_id), data = analytic_df)
m6 = lmer(glucosef ~ subtype*ns(t,df = 3) + censored_age + female + study + race  + dpp_intervention + (1|joint_id), data = analytic_df)


out1 = ggpredict(m1, c("t [all]", "subtype"))
out2 = ggpredict(m2, c("t [all]", "subtype"))
out3 = ggpredict(m3, c("t [all]", "subtype"))
out4 = ggpredict(m4, c("t [all]", "subtype"))
out5 = ggpredict(m5, c("t [all]", "subtype"))
out6 = ggpredict(m6, c("t [all]", "subtype"))

bind_rows(out1 %>% as.data.frame() %>%  mutate(outcome = "LDL"),
          out2 %>% as.data.frame() %>% mutate(outcome = "HDL"),
          out3 %>% as.data.frame() %>% mutate(outcome = "TG"),
          out4 %>% as.data.frame() %>% mutate(outcome = "SBP"),
          out5 %>% as.data.frame() %>% mutate(outcome = "eGFR"),
          out6 %>% as.data.frame() %>% mutate(outcome = "FPG")
) %>% 
  write_csv(.,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan09_modeled trajectories of additional biomarkers.csv"))


out_combined <- read_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan09_modeled trajectories of additional biomarkers.csv"))


# cluster_not2d_colors = c(cluster_colors,"#AC94F4")
cluster_not2d_colors = c(cluster_colors_cosmos,"#5C4033")
names(cluster_not2d_colors) = c(names(cluster_colors_cosmos),"No T2D")


fig_ldl = out_combined %>% 
  dplyr::filter(outcome == "LDL") %>% 
  mutate(cluster = factor(group,levels=c("NOT2D","MOD","SIRD","SIDD","MARD"),
                          labels=c("No T2D","MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=.,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,col=cluster)) +
  geom_path() +
  geom_ribbon(fill=NA,linetype = "dotted") +
  theme_bw() + 
  xlab("Time (years)") +
  ylab("LDL (mg/dL)") +
  scale_color_manual(name="",values=cluster_not2d_colors) 

fig_ldl

fig_hdl = out_combined %>% 
  dplyr::filter(outcome == "HDL") %>% 
  mutate(cluster = factor(group,levels=c("NOT2D","MOD","SIRD","SIDD","MARD"),
                          labels=c("No T2D","MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=.,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,col=cluster)) +
  geom_path() +
  geom_ribbon(fill=NA,linetype = "dotted") +
  theme_bw() + 
  xlab("Time (years)") +
  ylab("HDL (mg/dL)") +
  scale_color_manual(name="",values=cluster_not2d_colors) 

fig_tgl = out_combined %>% 
  dplyr::filter(outcome == "TG") %>% 
  mutate(cluster = factor(group,levels=c("NOT2D","MOD","SIRD","SIDD","MARD"),
                          labels=c("No T2D","MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=.,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,col=cluster)) +
  geom_path() +
  geom_ribbon(fill=NA,linetype = "dotted") +
  theme_bw() + 
  xlab("Time (years)") +
  ylab("Triglycerides (mg/dL)") +
  scale_color_manual(name="",values=cluster_not2d_colors) 

fig_sbp = out_combined %>% 
  dplyr::filter(outcome == "SBP") %>% 
  mutate(cluster = factor(group,levels=c("NOT2D","MOD","SIRD","SIDD","MARD"),
                          labels=c("No T2D","MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=.,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,col=cluster)) +
  geom_path() +
  geom_ribbon(fill=NA,linetype = "dotted") +
  theme_bw() + 
  xlab("Time (years)") +
  ylab("SBP (mmHg)") +
  scale_color_manual(name="",values=cluster_not2d_colors) 

fig_egfr = out_combined %>% 
  dplyr::filter(outcome == "eGFR") %>% 
  mutate(cluster = factor(group,levels=c("NOT2D","MOD","SIRD","SIDD","MARD"),
                          labels=c("No T2D","MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=.,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,col=cluster)) +
  geom_path() +
  geom_ribbon(fill=NA,linetype = "dotted") +
  theme_bw() + 
  xlab("Time (years)") +
  ylab("eGFR (mL/min/1.73 m²)") +
  scale_color_manual(name="",values=cluster_not2d_colors) 

fig_fpg = out_combined %>% 
  dplyr::filter(outcome == "FPG") %>% 
  mutate(cluster = factor(group,levels=c("NOT2D","MOD","SIRD","SIDD","MARD"),
                          labels=c("No T2D","MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=.,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,col=cluster)) +
  geom_path() +
  geom_ribbon(fill=NA,linetype = "dotted") +
  theme_bw() + 
  xlab("Time (years)") +
  ylab("FPG (mg/dL)") +
  scale_color_manual(name="",values=cluster_not2d_colors) 


library(ggpubr)
ggarrange(fig_ldl,
          fig_hdl,
          fig_tgl,
          fig_sbp,
          fig_egfr,
          fig_fpg,
          nrow=2,ncol=3,labels=LETTERS[1:6],common.legend = TRUE,legend = "bottom") %>% 
  ggsave(.,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/trajectory of additional biomarkers before diagnosis.jpg"),width = 9,height = 6)


