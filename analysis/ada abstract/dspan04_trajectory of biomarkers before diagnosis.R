rm(list = ls());gc();source(".Rprofile")

source("functions/egfr_ckdepi_2021.R")

# 1. Restrict to individuals with newly diagnosed T2D or no T2D. Set max_age as dmagediag if new T2D or max(age) if no T2D

# 2. Take their data retrospectively up to 15 years (i.e. max_age - age <= 15)
# Create a 't' variable: max_age - age
# Retain the wave of max_age (t = 0)

# 3. Include individuals who have 
# - at least one wave before T2D diagnosis (if newly diagnosed) AND have cluster info  
# - OR have at least two waves (if no T2D)
# Create an exposure variable called 'subtype' with 5 categories: SIDD, SIRD, MOD, MARD and NoT2D


# 4. Fit a model that accounts for clustering
# study_id: Variable that is a unique ID for each individual

# 5. Estimate predicted trajectories for each subtype over time along with confidence intervals
# https://stats.stackexchange.com/questions/109369/how-can-i-estimate-model-predicted-means-a-k-a-marginal-means-lsmeans-or-em

# new dm
final_dataset_temp = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/final_dataset_temp.RDS")) %>%
  mutate(joint_id = paste(study, original_study_id, sep = "_"))

analytic_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dsppre01_analytic df.RDS")) %>% 
  mutate(egfr_ckdepi_2021 = egfr_ckdepi_2021(scr = serumcreatinine,female = female,age = age),
         joint_id = paste(study, study_id, sep = "_"),
         # newly diagnosed T2D or no T2D
         newdm_event = case_when(joint_id %in% final_dataset_temp$joint_id ~ 1,
                           is.na(dmagediag) ~ 0,
                           TRUE ~ NA)) %>%
  dplyr::filter(!is.na(newdm_event)) %>% 
  group_by(study,study_id,joint_id) %>% 
  mutate(max_age = case_when(newdm_event == 1 ~ dmagediag,
                             TRUE ~ max(age)),
         t = age - max_age) %>% 
  dplyr::filter(t <= 0 & t >= -15) %>% 
  ungroup()

wave_df <- analytic_df %>% 
  group_by(study,study_id,joint_id) %>%
  mutate(has_age_before_max = any(age < max_age)) %>% 
  dplyr::filter((newdm_event == 1 & has_age_before_max & !is.na(cluster)) 
                | (newdm_event == 0 & has_age_before_max)) %>% 
  mutate(subtype = case_when(is.na(dmagediag) ~ "NOT2D",
                             !is.na(cluster) ~ cluster,
                             TRUE ~ NA_character_)) %>% 
  ungroup() %>% 
  arrange(joint_id,t) %>% 
  distinct(joint_id,t,.keep_all=TRUE) %>% 
  mutate(across(one_of(c("subtype","race","study")),.fns=~as.factor(.)))
  
summary(wave_df[,c("homa2b","t","subtype","max_age","female","study","race")])
table(wave_df$subtype,useNA="always")

wave_df %>% 
  group_by(subtype) %>% 
  distinct(study,study_id,joint_id) %>% 
  tally() %>% 
  mutate(p = n/sum(n))

ind_df = wave_df %>% 
  group_by(study,study_id) %>% 
  dplyr::filter(t == max(t)) %>% 
  ungroup()

ind_df %>% 
  group_by(subtype == "NOT2D") %>% 
  summarize(age = mean(max_age),
            sd_age = sd(max_age),
            female = mean(female))

library(splines)
library(emmeans)
library(ggeffects)
emm_options(pbkrtest.limit = 999999)

# library(geepack)
# m1a = geeglm(homa2b ~ subtype*ns(t) + max_age + female + study + race, family = gaussian(),data = wave_df,id = joint_id,corstr = "independence")
# m2 = geeglm(bmi ~ subtype*ns(t) + max_age + female + study + race, family = gaussian(),data = wave_df,id = joint_id,corstr = "independence")
# m3 = geeglm(hba1c ~ subtype*ns(t) + max_age + female + study + race, family = gaussian(),data = wave_df,id = joint_id,corstr = "independence")
 

library(lme4)
m1 = lmer(log(homa2b) ~ subtype*ns(t, df = 3) + max_age + female + study + race + (1|joint_id), data = wave_df)
m2 = lmer(bmi ~ subtype*ns(t,df = 3) + max_age + female + study + race  + (1|joint_id), data = wave_df)
m3 = lmer(hba1c ~ subtype*ns(t,df = 3) + max_age + female + study + race  + (1|joint_id), data = wave_df)
m4 = lmer(log(homa2ir) ~ subtype*ns(t,df = 3) + max_age + female + study + race  + (1|joint_id), data = wave_df)



# https://stackoverflow.com/questions/77811491/how-to-perform-piecewise-linear-mixed-regression-with-multiple-breakpoints-in-r

# The below take a long time to fit -- ~ 5 mins per row -------------------
out1 = ggpredict(m1, c("t [all]", "subtype"))
# Message: Model has log-transformed response. Back-transforming predictions to original response scale. Standard errors are still on the transformed scale.
out2 = ggpredict(m2, c("t [all]", "subtype"))
out3 = ggpredict(m3, c("t [all]", "subtype"))
out4 = ggpredict(m4, c("t [all]", "subtype"))

bind_rows(out1 %>% as.data.frame() %>%  mutate(outcome = "HOMA2B"),
          out2 %>% as.data.frame() %>% mutate(outcome = "BMI"),
          out3 %>% as.data.frame() %>% mutate(outcome = "HbA1c"),
          out4 %>% as.data.frame() %>% mutate(outcome = "HOMA2IR")
          ) %>% 
  write_csv(.,paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan04_modeled trajectories of biomarkers.csv"))


out_combined <- read_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan04_modeled trajectories of biomarkers.csv"))

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
  scale_color_manual(name="",values=cluster_not2d_colors)

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
  ggsave(.,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/trajectory of biomarkers before diagnosis.jpg"),width = 6,height = 6)







out_combined <- read_csv(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan04_modeled trajectories of biomarkers.csv")) %>% 
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


















































