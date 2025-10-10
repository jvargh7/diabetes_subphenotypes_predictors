evaluate_coxph <- function(model, data, time_horizon,
                           cutpoints = c(0.10, 0.25, 0.50, 0.75, 0.90)) {
  stopifnot(inherits(model, "coxph"))
  
  # LP on eval data
  lp <- predict(model, newdata = data, type = "lp")
  
  # Baseline survival S0(t) at horizon
  bh <- basehaz(model, centered = FALSE)
  S0 <- approx(bh$time, exp(-bh$hazard), xout = time_horizon, rule = 2)$y
  
  # Predicted survival & risk at horizon
  surv_prob <- S0^exp(lp)
  risk_prob <- 1 - surv_prob
  
  # Observed event by horizon
  obs_event <- with(data, ifelse(!is.na(time_to_event) & !is.na(event) &
                                   time_to_event <= time_horizon & event == 1, 1L, 0L))
  
  # Valid for evaluation: enough follow-up OR had the event before horizon
  valid <- with(data, !is.na(time_to_event) & !is.na(event) &
                  (time_to_event >= time_horizon | event == 1))
  
  # Guard: if nothing valid, return NAs
  if (!any(valid)) {
    out <- tibble::tibble(
      threshold   = cutpoints,
      sensitivity = NA_real_,
      specificity = NA_real_,
      F1          = NA_real_,
      c_index     = NA_real_,
      cal_slope   = NA_real_
    )
    return(out)
  }
  
  # Threshold metrics
  res <- lapply(cutpoints, function(thresh) {
    pred <- as.integer(risk_prob >= thresh)
    TP <- sum(pred[valid] == 1 & obs_event[valid] == 1, na.rm = TRUE)
    TN <- sum(pred[valid] == 0 & obs_event[valid] == 0, na.rm = TRUE)
    FP <- sum(pred[valid] == 1 & obs_event[valid] == 0, na.rm = TRUE)
    FN <- sum(pred[valid] == 0 & obs_event[valid] == 1, na.rm = TRUE)
    
    sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
    f1   <- if ((2*TP + FP + FN) > 0) 2 * TP / (2*TP + FP + FN) else NA_real_
    
    c(threshold = thresh, sensitivity = sens, specificity = spec, F1 = f1)
  }) %>% dplyr::bind_rows()
  
  # Concordance on model (same as your original)
  cidx <- tryCatch(concordance(model)$concordance, error = function(e) NA_real_)
  
  # Calibration slope (with intercept, as you had)
  cal_slope <- tryCatch({
    coef(glm(obs_event[valid] ~ lp[valid], family = binomial()))[2]
  }, error = function(e) NA_real_)
  
  res$c_index   <- cidx
  res$cal_slope <- cal_slope
  res
}
