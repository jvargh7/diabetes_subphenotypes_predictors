
evaluate_coxph <- function(model, data, time_horizon, cutpoints = c(0.10, 0.25, 0.50, 0.75, 0.90)) {
  # Predicted linear predictor
  lp <- predict(model, newdata = data, type = "lp")
  
  # Get baseline survival
  base_surv <- basehaz(model, centered = FALSE)
  S0 <- approx(base_surv$time, exp(-base_surv$hazard), xout = time_horizon, rule = 2)$y
  
  # Predicted survival and risk
  surv_prob <- S0^exp(lp)
  risk_prob <- 1 - surv_prob
  
  # Observed event within time_horizon
  obs_event <- with(data, ifelse(time_to_event <= time_horizon & event == 1, 1, 0))
  
  # Filter only those with follow-up time >= time_horizon
  valid <- data$time_to_event >= time_horizon | data$event == 1
  
  results <- lapply(cutpoints, function(thresh) {
    pred <- ifelse(risk_prob >= thresh, 1, 0)
    
    TP <- sum(pred == 1 & obs_event == 1 & valid)
    TN <- sum(pred == 0 & obs_event == 0 & valid)
    FP <- sum(pred == 1 & obs_event == 0 & valid)
    FN <- sum(pred == 0 & obs_event == 1 & valid)
    
    sensitivity <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    f1 <- ifelse((2 * TP + FP + FN) == 0, NA, 2 * TP / (2 * TP + FP + FN))
    
    return(c(threshold = thresh, sensitivity = sensitivity, specificity = specificity, F1 = f1))
  }) %>% bind_rows()
  
  # C-index
  cidx <- concordance(model)$concordance
  
  # Calibration slope
  cal_model <- glm(obs_event[valid] ~ lp[valid], family = binomial)
  cal_slope <- coef(cal_model)[2]
  
  results$c_index <- cidx
  results$cal_slope <- cal_slope
  return(results)
}
