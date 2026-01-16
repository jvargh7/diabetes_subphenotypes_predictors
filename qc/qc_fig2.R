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
m1 = lmer(homa2b ~ subtype*ns(t, df = 3) + censored_age + female + study + race + dpp_intervention + (1|joint_id), data = analytic_df)
m2 = lmer(bmi ~ subtype*ns(t,df = 3) + censored_age + female + study + race  + dpp_intervention + (1|joint_id), data = analytic_df)
m3 = lmer(hba1c ~ subtype*ns(t,df = 3) + censored_age + female + study + race  + dpp_intervention + (1|joint_id), data = analytic_df)
m4 = lmer(homa2ir ~ subtype*ns(t,df = 3) + censored_age + female + study + race  + dpp_intervention + (1|joint_id), data = analytic_df)



# https://stackoverflow.com/questions/77811491/how-to-perform-piecewise-linear-mixed-regression-with-multiple-breakpoints-in-r

# The below take a long time to fit -- ~ 5 mins per row -------------------
out1 = ggpredict(m1, c("t [all]", "subtype"))
out2 = ggpredict(m2, c("t [all]", "subtype"))
out3 = ggpredict(m3, c("t [all]", "subtype"))
out4 = ggpredict(m4, c("t [all]", "subtype"))

bind_rows(out1 %>% as.data.frame() %>%  mutate(outcome = "HOMA2B"),
          out2 %>% as.data.frame() %>% mutate(outcome = "BMI"),
          out3 %>% as.data.frame() %>% mutate(outcome = "HbA1c"),
          out4 %>% as.data.frame() %>% mutate(outcome = "HOMA2IR")
) %>% 
  write_csv(.,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/qc/dspan02_modeled trajectories of biomarkers.csv"))


out_combined <- read_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/qc/dspan02_modeled trajectories of biomarkers.csv"))


# cluster_not2d_colors = c(cluster_colors,"#AC94F4")
cluster_not2d_colors = c(cluster_colors_cosmos,"#5C4033")
names(cluster_not2d_colors) = c(names(cluster_colors_cosmos),"No T2D")

fig_homa2b = out_combined %>% 
  dplyr::filter(outcome == "HOMA2B") %>% 
  mutate(cluster = factor(group,levels=c("NOT2D","MOD","SIRD","SIDD","MARD"),
                          labels=c("No T2D","MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=.,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,col=cluster)) +
  geom_path() +
  geom_ribbon(fill=NA,linetype = "dotted") +
  theme_bw() + 
  xlab("Time (years)") +
  ylab("HOMA2-%B") +
  scale_color_manual(name="",values=cluster_not2d_colors) 

fig_homa2b

fig_bmi = out_combined %>%
  dplyr::filter(outcome == "BMI") %>% 
  
  mutate(cluster = factor(group,levels=c("NOT2D","MOD","SIRD","SIDD","MARD"),
                          labels=c("No T2D","MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=.,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,col=cluster)) +
  geom_path() +
  geom_ribbon(fill=NA,linetype = "dotted") +
  theme_bw() + 
  xlab("Time (years)") +
  ylab("BMI (kg/m2)") +
  scale_color_manual(name="",values=cluster_not2d_colors) +
  coord_cartesian(ylim = c(20, 35))

fig_hba1c = out_combined %>%
  dplyr::filter(outcome == "HbA1c") %>% 
  
  mutate(cluster = factor(group,levels=c("NOT2D","MOD","SIRD","SIDD","MARD"),
                          labels=c("No T2D","MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=.,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,col=cluster)) +
  geom_path() +
  geom_ribbon(fill=NA,linetype = "dotted") +
  theme_bw() + 
  xlab("Time (years)") +
  ylab("HbA1c (%)") +
  scale_color_manual(name="",values=cluster_not2d_colors)

fig_homa2ir = out_combined %>%
  dplyr::filter(outcome == "HOMA2IR") %>% 
  mutate(cluster = factor(group,levels=c("NOT2D","MOD","SIRD","SIDD","MARD"),
                          labels=c("No T2D","MOD","SIRD","SIDD","MARD"))) %>% 
  ggplot(data=.,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,col=cluster)) +
  geom_path() +
  geom_ribbon(fill=NA,linetype = "dotted") +
  theme_bw() + 
  xlab("Time (years)") +
  ylab("HOMA2-IR") +
  scale_color_manual(name="",values=cluster_not2d_colors) 



library(ggpubr)
ggarrange(fig_bmi,
          fig_hba1c,
          fig_homa2b,
          fig_homa2ir,
          nrow=2,ncol=2,labels=LETTERS[1:4],common.legend = TRUE,legend = "bottom") %>% 
  ggsave(.,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/qc/trajectory of biomarkers before diagnosis.jpg"),width = 6,height = 6)

ggarrange(fig_bmi,
          fig_hba1c,
          fig_homa2b,
          fig_homa2ir,
          nrow=2,ncol=2,labels=LETTERS[1:4],common.legend = TRUE,legend = "bottom") %>% 
  ggsave(.,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/qc/trajectory of biomarkers before diagnosis.pdf"),width = 8,height = 8)




out_combined <- read_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/qc/dspan04_modeled trajectories of biomarkers.csv")) %>% 
  dplyr::filter(x %in% c(-5,0)) %>% 
  dplyr::select(group,outcome,x,predicted,std.error) %>% 
  mutate(t = case_when(x == -5 ~ "t5",
                       x == 0 ~ "t0")) %>% 
  dplyr::select(-x) %>% 
  pivot_wider(values_from = c("predicted","std.error"),names_from=c("t")) %>% 
  mutate(diff = (predicted_t0 - predicted_t5),
         se_diff = sqrt(std.error_t5^2 + std.error_t0^2)) %>% 
  mutate(diff_lci = diff - 1.96*se_diff,
         diff_uci = diff + 1.96*se_diff)

not2d_combined = out_combined %>% 
  dplyr::filter(group == "NOT2D") %>% 
  dplyr::select(group,outcome,diff,se_diff) %>% 
  rename(not2d_diff = diff,
         not2d_se_diff = se_diff)


out_combined_not2d = out_combined %>% 
  left_join(not2d_combined %>% dplyr::select(-group),
            by = c("outcome")) %>% 
  mutate(diff_not2d_diff = diff - not2d_diff,
         se_diff_not2d_diff = sqrt(se_diff^2 + not2d_se_diff^2)) %>% 
  mutate(diff_not2d_diff_lci = diff_not2d_diff - 1.96*se_diff_not2d_diff,
         diff_not2d_diff_uci = diff_not2d_diff + 1.96*se_diff_not2d_diff)


