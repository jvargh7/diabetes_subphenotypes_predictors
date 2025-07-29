
get_calibration_slope <- function(data, event_name) {
  # Step 1: Fit original Cox model
  formula_str <- paste0("Surv(tstart, tstop, ", event_name, ") ~ study + female + race + earliest_age + bmi + hba1c + ",
                        "homa2b_scaled + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention")
  surv_formula <- as.formula(formula_str)
  
  fit <- try(coxph(surv_formula, data = data, ties = "breslow", control = coxph.control(iter.max = 20)), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA)
  
  # Step 2: Get linear predictor
  data$lp <- predict(fit, type = "lp")
  
  # Step 3: Fit calibration slope model
  slope_fit <- try(coxph(Surv(tstart, tstop, data[[event_name]]) ~ lp, data = data), silent = TRUE)
  if (inherits(slope_fit, "try-error")) return(NA)
  
  # Step 4: Return slope coefficient
  slope <- coef(slope_fit)["lp"]
  return(slope)
}
