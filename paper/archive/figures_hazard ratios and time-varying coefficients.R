rm(list = ls());gc();source(".Rprofile")

library(gridExtra)
library(ggplot2)
library(survival)

# time-varying coefficients plot
# Ref: https://pmc.ncbi.nlm.nih.gov/articles/PMC6015946/#sec2


ipcw_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan02_ipcw dfs.RDS"))

for (i in 1:length(ipcw_dfs)) {
  df <- ipcw_dfs[[1]]
  
  tdcm_df <- df %>%
    arrange(study, study_id, age) %>%
    group_by(study, study_id) %>%
    mutate(
      # tstart = case_when(row_number() == 1 ~ age, 
      #                   TRUE ~ dplyr::lag(age, n = 1)), 
      # tstop = age
      baseline_age = first(age),  # Assuming 'age' at first observation is baseline
      tstart = case_when(
        row_number() == 1 ~ 0, 
        TRUE ~ age - first(age)
      ), 
      tstop = lead(tstart, default = last(age) - first(age)) 
    ) %>%
    ungroup() %>% 
    mutate(across(c(bmi, hba1c, homa2b, homa2ir, ldlc, sbp, egfr_ckdepi_2021), ~replace(., is.infinite(.), NA))) %>% 
    dplyr::filter((tstart < tstop) & (tstop <= censored_age)) %>% 
    # error due to 0 ppl in NH Other (sidd == 1), ignore this category
    mutate(race = case_when(race == "NH Other" ~ "Other", 
                            TRUE ~ race)) %>% 
    mutate(cluster = factor(cluster,
                            levels = c("MOD", "SIDD", "MARD", "SIRD"),
                            labels = c("MOD", "SIDD", "MARD", "SIRD"))
    ) 
  

}


bmi_fit<- coxph(Surv(tstart, tstop, event) ~ bmi + study + female + race + min_age, 
                      data = tdcm_df, weights = ipcw_cluster)

hba1c_fit <- coxph(Surv(tstart, tstop, event) ~ hba1c + study + female + race + min_age, 
                        data = tdcm_df, weights = ipcw_cluster)

homa2b_fit <- coxph(Surv(tstart, tstop, event) ~ homa2b + study + female + race + min_age, 
                         data = tdcm_df, weights = ipcw_cluster)

homa2ir_fit <- coxph(Surv(tstart, tstop, event) ~ homa2ir + study + female + race + min_age, 
                          data = tdcm_df, weights = ipcw_cluster)

ldlc_fit <- coxph(Surv(tstart, tstop, event) ~ ldlc + study + female + race + min_age, 
                       data = tdcm_df, weights = ipcw_cluster)

sbp_fit <- coxph(Surv(tstart, tstop, event) ~ sbp + study + female + race + min_age, 
                      data = tdcm_df, weights = ipcw_cluster)

egfr_fit <- coxph(Surv(tstart, tstop, event) ~ egfr_ckdepi_2021 + study + female + race + min_age, 
                       data = tdcm_df, weights = ipcw_cluster)


zph_bmi <- cox.zph(bmi_fit)
zph_hba1c <- cox.zph(hba1c_fit)
zph_homa2b <- cox.zph(homa2b_fit)
zph_homa2ir <- cox.zph(homa2ir_fit)
zph_ldlc <- cox.zph(ldlc_fit)
zph_sbp <- cox.zph(sbp_fit)
zph_egfr <- cox.zph(egfr_fit)



jpeg(paste0(path_diabetes_subphenotypes_predictors_folder,"/figures/time-varying coefficients.jpg"), width = 1600, height = 800)

# Setting up the plotting area for 2 rows and 4 columns
par(mfrow = c(2, 4), mar = c(4, 4, 2, 1) + 0.1)

plot(zph_bmi[1], lwd = 2)
abline(0, 0, col=1, lty=3, lwd=2)
abline(h=zph_bmi$coef[1],col=3,lwd=2,lty=2)

plot(zph_hba1c[1], lwd = 2)
abline(0, 0, col=1, lty=3, lwd=2)
abline(h=zph_hba1c$coef[1],col=3,lwd=2,lty=2)

plot(zph_homa2b[1], lwd = 2)
abline(0, 0, col=1, lty=3, lwd=2)
abline(h=zph_homa2b$coef[1],col=3,lwd=2,lty=2)

plot(zph_homa2ir[1], lwd = 2)
abline(0, 0, col=1, lty=3, lwd=2)
abline(h=zph_homa2ir$coef[1],col=3,lwd=2,lty=2)

plot(zph_ldlc[1], lwd = 2)
abline(0, 0, col=1, lty=3, lwd=2)
abline(h=zph_ldlc$coef[1],col=3,lwd=2,lty=2)

plot(zph_sbp[1], lwd = 2)
abline(0, 0, col=1, lty=3, lwd=2)
abline(h=zph_sbp$coef[1],col=3,lwd=2,lty=2)

plot(zph_egfr[1], lwd = 2)
abline(0, 0, col=1, lty=3, lwd=2)
abline(h=zph_egfr$coef[1],col=3,lwd=2,lty=2)

dev.off()

# Reset the plotting layout if you're going to continue using the R session
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)







