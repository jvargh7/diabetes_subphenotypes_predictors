evaluate_coxph_quantile <- function(model, data, time_horizon, quantiles = c(0.10, 0.25, 0.50, 0.75, 0.90)) {
  # Predicted linear predictor
  lp <- predict(model, newdata = data, type = "lp")
  
  # Get baseline survival
  base_surv <- survival::basehaz(model, centered = FALSE)
  S0 <- approx(base_surv$time, exp(-base_surv$hazard), xout = time_horizon, rule = 2)$y
  
  # Predicted survival and risk
  surv_prob <- S0^exp(lp)
  risk_prob <- 1 - surv_prob
  
  # Observed event within time_horizon
  obs_event <- with(data, ifelse(!is.na(time_to_event) & !is.na(event) &
                                   time_to_event <= time_horizon & event == 1, 1, 0))
  
  # Valid for evaluation: enough follow-up OR had the event
  valid <- with(data, !is.na(time_to_event) & !is.na(event) &
                  (time_to_event >= time_horizon | event == 1))
  
  # Guard: if nothing valid, return NAs
  if (!any(valid)) {
    results <- data.frame(
      threshold   = quantiles,
      sensitivity = NA_real_,
      specificity = NA_real_,
      F1          = NA_real_,
      c_index     = NA_real_,
      cal_slope   = NA_real_
    )
    return(results)
  }
  
  # Use quantile-based thresholds instead of fixed thresholds
  risk_quantiles <- quantile(risk_prob[valid], probs = 1 - quantiles, na.rm = TRUE)
  
  results <- lapply(seq_along(quantiles), function(idx) {
    thresh <- risk_quantiles[idx]
    pred <- ifelse(risk_prob >= thresh, 1, 0)
    
    # Only use valid observations for metrics
    TP <- sum(pred[valid] == 1 & obs_event[valid] == 1, na.rm = TRUE)
    TN <- sum(pred[valid] == 0 & obs_event[valid] == 0, na.rm = TRUE)
    FP <- sum(pred[valid] == 1 & obs_event[valid] == 0, na.rm = TRUE)
    FN <- sum(pred[valid] == 0 & obs_event[valid] == 1, na.rm = TRUE)
    
    sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
    f1 <- if ((2 * TP + FP + FN) > 0) 2 * TP / (2 * TP + FP + FN) else NA_real_
    
    return(c(threshold = quantiles[idx], sensitivity = sensitivity, specificity = specificity, F1 = f1))
  }) %>% bind_rows()
  
  # C-index
  cidx <- tryCatch(survival::concordance(model)$concordance, error = function(e) NA_real_)
  
  # Calibration slope - only use valid observations
  cal_slope <- tryCatch({
    if (sum(valid) > 1) {
      cal_model <- glm(obs_event[valid] ~ lp[valid], family = binomial)
      coef(cal_model)[2]
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)
  
  results$c_index <- cidx
  results$cal_slope <- cal_slope
  return(results)
}