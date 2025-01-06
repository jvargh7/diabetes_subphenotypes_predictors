rm(list = ls());gc();source(".Rprofile")

library(gridExtra)
library(ggplot2)
library(survminer)

tdcm_coef <- read_csv("analysis/dspan03_tdcm pooled results with multiple imputation.csv") %>% 
  select(iv, estimate, lci, uci, model) %>% 
  mutate(HR = paste0(round(estimate, 2), " (", round(lci, 2), ", ", round(uci, 2), ")")) %>% 
  dplyr::filter(!iv %in% c("studymesa","studyjhs","raceNH Black","raceNH White","raceOther","female1","min_age"),
                model != "Overall") %>% 
  mutate(term = case_when(
    iv == "bmi" ~ "BMI",
    iv == "sbp" ~ "SBP",
    iv == "hba1c" ~ "HbA1c",
    iv == "ldlc" ~ "LDL",
    iv == "homa2b" ~ "HOMA2-%B",
    iv == "homa2ir" ~ "HOMA2-IR",
    iv == "egfr_ckdepi_2021" ~ "eGFR",
    TRUE ~ iv  
  ),
  term = factor(term,
                levels = c("eGFR", "HOMA2-IR", "HOMA2-%B", "LDL", "HbA1c", "SBP", "BMI"),
                labels = c("eGFR", "HOMA2-IR", "HOMA2-%B", "LDL", "HbA1c", "SBP", "BMI"))
  )


# forest plot

plot_forest <- ggplot(tdcm_coef, aes(y = term, x = estimate, xmin = lci, xmax = uci, color = model)) + 
  geom_pointrange(position = position_dodge(width = 0.7), size = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(values = cluster_colors) +
  scale_x_continuous(limits = c(0, 2.5)) +
  labs(
    x = "Hazard ratio (95% CI)",
    y = NULL,
    title = "A: Hazard Ratios by Pathophysiological Markers",
    color = "Subtype"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.4)
  ) +
  geom_text(
    aes(x = uci + 0.05, label = HR),
    position = position_dodge(width = 0.7),
    vjust = 0.2,
    hjust = -0.05,
    fontface = "bold",
    size = 4
  ) 



#-------------------------------------------------------------------------------------------------------------------
# time-varying coefficients plot

ipcw_dfs <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/processed/dspan02_ipcw dfs.RDS"))

for (i in 1:length(ipcw_dfs)) {
  df <- ipcw_dfs[[i]]
  
  tdcm_df <- df %>%
    arrange(study, study_id, age) %>%
    group_by(study, study_id) %>%
    mutate(
      tstart = case_when(row_number() == 1 ~ age, 
                        TRUE ~ dplyr::lag(age, n = 1)), 
      tstop = age
      # baseline_age = first(age),  # Assuming 'age' at first observation is baseline
      # tstart = case_when(
      #   row_number() == 1 ~ 0, 
      #   TRUE ~ age - first(age)
      # ), 
      # tstop = lead(tstart, default = last(age) - first(age)) 
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
  
  
  bmi_fit[[i]] <- coxph(Surv(tstart, tstop, event) ~ bmi + study + female + race + min_age, 
                         data = tdcm_df, weights = ipcw_cluster)
  
  hba1c_fit[[i]] <- coxph(Surv(tstart, tstop, event) ~ hba1c + study + female + race + min_age, 
                        data = tdcm_df, weights = ipcw_cluster)
  
  homa2b_fit[[i]] <- coxph(Surv(tstart, tstop, event) ~ homa2b + study + female + race + min_age, 
                        data = tdcm_df, weights = ipcw_cluster)
  
  homa2ir_fit[[i]] <- coxph(Surv(tstart, tstop, event) ~ homa2ir + study + female + race + min_age, 
                        data = tdcm_df, weights = ipcw_cluster)
  
  ldlc_fit[[i]] <- coxph(Surv(tstart, tstop, event) ~ ldlc + study + female + race + min_age, 
                        data = tdcm_df, weights = ipcw_cluster)
  
  sbp_fit[[i]] <- coxph(Surv(tstart, tstop, event) ~ sbp + study + female + race + min_age, 
                        data = tdcm_df, weights = ipcw_cluster)
  
  egfr_fit[[i]] <- coxph(Surv(tstart, tstop, event) ~ egfr_ckdepi_2021 + study + female + race + min_age, 
                        data = tdcm_df, weights = ipcw_cluster)
  
  
  
}

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







