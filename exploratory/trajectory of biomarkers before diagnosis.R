rm(list = ls());gc();source(".Rprofile")


final_dataset_temp = readRDS(paste0(path_diabetes_subphenotypes_adults_folder,"/working/cleaned/final_dataset_temp.RDS"))

clusters = read_csv(paste0(path_diabetes_subphenotypes_adults_folder,"/working/processed/dec_an02_clean_kmeans_5var_mi_knn_cluster.csv")) %>% 
  dplyr::select(-one_of("...1")) %>% 
  left_join(final_dataset_temp %>% 
              dplyr::select(study_id,original_study_id),
            by=c("study_id")) %>% 
  rename(cluster_study_id = study_id)

aric_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01a_aric.RDS")) %>% 
  # dplyr::filter((age < dmagediag))
  dplyr::filter(!is.na(dmagediag))

cardia_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01b_cardia.RDS")) %>% 
  # dplyr::filter((age < dmagediag))
  dplyr::filter(!is.na(dmagediag))

jhs_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01e_jhs.RDS")) %>% 
  # dplyr::filter((age < dmagediag))
  dplyr::filter(!is.na(dmagediag))

dppos_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01c_dppos.RDS")) %>% 
  # dplyr::filter((age < dmagediag))
  dplyr::filter(!is.na(dmagediag))

mesa_longitudinal <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01f_mesa.RDS")) %>% 
  # dplyr::filter((age < dmagediag))
  dplyr::filter(!is.na(dmagediag))


longitudinal_df = bind_rows(aric_longitudinal %>% mutate(study = "aric") %>% mutate(study_id = as.numeric(str_replace(study_id,"C",""))),
                            cardia_longitudinal %>% mutate(study = "cardia"),
                            jhs_longitudinal %>% mutate(study = "jhs"),
                            dppos_longitudinal %>% mutate(study = "dppos"),
                            mesa_longitudinal %>% mutate(study = "mesa")) %>%
  
  bind_rows(final_dataset_temp %>% dplyr::select(-study_id) %>% rename(study_id = original_study_id) %>% 
              dplyr::select(one_of(unique(c(
                                          names(aric_longitudinal),
                                          names(cardia_longitudinal),
                                          names(jhs_longitudinal),
                                          names(dppos_longitudinal),
                                          names(mesa_longitudinal))
                                     
                                     )
                                   ))) %>% 
  
  mutate(glucosef2 = case_when(is.na(glucosef2) ~ glucosef*0.0555,
                               TRUE ~ glucosef2),
         # 1 μIU/mL = 6.00 pmol/L
         insulinf2 = case_when(is.na(insulinf2) ~ insulinf*6,
                               TRUE ~ insulinf2)) %>% 
  mutate(glucosef = case_when(is.na(glucosef) ~ glucosef2/0.0555,
                              TRUE ~ glucosef),
         insulinf = case_when(is.na(insulinf) ~ insulinf2/6,
                              TRUE ~ insulinf)
         
         ) %>% 

  left_join(clusters %>% 
              dplyr::select(cluster_study_id,original_study_id,cluster,study,female),
            by=c("study"="study","study_id" = "original_study_id")) %>% 
  dplyr::filter(!is.na(cluster_study_id)) %>% 
  mutate(t = age - dmagediag) %>% 
  dplyr::filter(t >= -15, t < 0)

# Not vectorized yet -- 
source("C:/code/external/functions/cgm/egfr_ckdepi_2021.R")
longitudinal_df$egfr_2021 = egfr_ckdepi_2021(longitudinal_df$serumcreatinine,longitudinal_df$female,longitudinal_df$age)

clusters_in_longitudinal <- clusters %>% 
  dplyr::filter(cluster_study_id %in% longitudinal_df$cluster_study_id)

table(clusters_in_longitudinal$cluster)
table(clusters_in_longitudinal$cluster) %>% prop.table()

table(longitudinal_df$cluster)
table(longitudinal_df$cluster) %>% prop.table()

# EXPLORATORY -----------
ggplot(data=longitudinal_df,
       aes(x=t,y=bmi,col=cluster,group=cluster)) +
  # geom_point() +
  geom_smooth()

 ggplot(data=longitudinal_df,
       aes(x=t,y=hba1c,col=cluster,group=cluster)) +
  # geom_point() +
  geom_smooth()

ggplot(data=longitudinal_df,
       aes(x=t,y=glucosef,col=cluster,group=cluster)) +
  # geom_point() +
  geom_smooth(method = "loess")


ggplot(data=longitudinal_df,
       aes(x=t,y=insulinf,col=cluster,group=cluster)) +
  # geom_point() +
  geom_smooth()

# ADJUSTED MODELS -------------
library(geepack)
library(marginaleffects)
library(splines)
m1_bmi = geeglm(bmi ~ cluster*ns(t,1) + dmagediag + female + study,data = longitudinal_df,id = cluster_study_id)
m1_hba1c = geeglm(hba1c ~ cluster*ns(t,1) + dmagediag + female + study,data = longitudinal_df,id = cluster_study_id)
m1_glucosef = geeglm(glucosef ~ cluster*ns(t,2) + dmagediag + female + study,data = longitudinal_df,id = cluster_study_id)
m1_insulinf = geeglm(insulinf ~ cluster*ns(t,3) + dmagediag + female + study,data = longitudinal_df,id = cluster_study_id)

m1_ldlc = geeglm(ldlc ~ cluster*ns(t,2) + dmagediag + female + study,data = longitudinal_df,id = cluster_study_id)
m1_hdlc = geeglm(hdlc ~ cluster*ns(t,2) + dmagediag + female + study,data = longitudinal_df,id = cluster_study_id)
m1_tgl = geeglm(tgl ~ cluster*ns(t,2) + dmagediag + female + study,data = longitudinal_df,id = cluster_study_id)
m1_egfr_2021 = geeglm(egfr_2021 ~ cluster*ns(t,2) + dmagediag + female + study,data = longitudinal_df,id = cluster_study_id)




predictions_trajectory = function(model){
  t_predict_range = seq(-15,-1,by=0.5)
  t_display_range = c(-15,0)
  # t_predict_range = seq(-15,6,by=0.5) 
  # t_display_range = c(-15,6)
  t_display_breaks = seq(t_display_range[1],t_display_range[2],by=3)
  
  predictions(
    model,
    newdata = datagrid(cluster = unique,t=t_predict_range)) %>% 
    return()
  
  
}

plot_predictions <- function(pred_out){
  
  ggplot(data=pred_out,aes(x=t,ymin=conf.low,ymax=conf.high,y=estimate,col=cluster,fill=cluster,group=cluster)) +
    geom_path() + 
    geom_ribbon(alpha=0.5) +
    xlab("Years before diagnosis") +
    theme_bw() + 
    scale_x_continuous(limits=t_display_range,breaks=t_display_breaks) +
    scale_color_manual("",values=cluster_colors) +
    scale_fill_manual("",values=cluster_colors) +
    geom_vline(xintercept = c(0,1),col="red",linetype = 2)
  
  
}


## BMI ---------
# https://marginaleffects.com/vignettes/marginalmeans.html
p1_bmi <- predictions_trajectory(m1_bmi)
p1_hba1c <- predictions_trajectory(m1_hba1c)
p1_glucosef <- predictions_trajectory(m1_glucosef)
p1_insulinf <- predictions_trajectory(m1_insulinf)
p1_ldlc <- predictions_trajectory(m1_ldlc)
p1_hdlc <- predictions_trajectory(m1_hdlc)
p1_tgl <- predictions_trajectory(m1_tgl)
p1_egfr_2021 <- predictions_trajectory(m1_egfr_2021)

f1_bmi <- plot_predictions(p1_bmi) +
  ylab(bquote('BMI ( kg' /m^2~')'))

f1_hba1c <- plot_predictions(p1_hba1c) +
  ylab("HbA1c (%)")

f1_glucosef <- plot_predictions(p1_glucosef) +
  ylab("Fasting Glucose (mg/dL)")

f1_insulinf <- plot_predictions(p1_insulinf) +
  ylab("Fasting Insulin (μIU/mL)")

f1_ldlc <- plot_predictions(p1_ldlc) +
  ylab("LDL cholesterol (mg/dL)")

f1_hdlc <- plot_predictions(p1_hdlc) +
  ylab("HDL cholesterol (mg/dL)")

f1_tgl <- plot_predictions(p1_tgl) +
  ylab("Triglyceride (mg/dL)")

f1_egfr_2021 <- plot_predictions(p1_egfr_2021) +
  ylab(bquote('eGFR ( mL/min/1.73' /m^2~')'))




## Combined Figure ----------
library(ggpubr)
ggarrange(f1_bmi,
          f1_hba1c,
          f1_glucosef,
          f1_insulinf,
          f1_ldlc,
          f1_hdlc,
          f1_tgl,
          f1_egfr_2021,
          common.legend = TRUE,
          legend = "bottom",
          nrow = 2,
          ncol = 4) %>% 
  ggsave(.,filename=paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/trajectory of biomarkers before diagnosis.jpg"),width=12,height = 8)

