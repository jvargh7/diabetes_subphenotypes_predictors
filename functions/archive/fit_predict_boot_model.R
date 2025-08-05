
fit_predict_boot_model <- function(data, indices, event_name, m, b) {
  unique_ids <- unique(data$joint_id)
  sampled_ids <- unique_ids[indices]
  d <- data[data$joint_id %in% sampled_ids, ]
  
  # Skip if only 1 class present
  if (length(unique(d[[event_name]])) < 2 || length(unique(d$race)) == 1 || length(unique(d$female)) == 1) {
    return(NULL)
  }
  
  # Build survival formula
  formula_str <- paste0("Surv(tstart, tstop, ", event_name, ") ~ study + female + race + earliest_age + bmi + hba1c + ",
                        "homa2b_scaled + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention")
  surv_formula <- as.formula(formula_str)
  
  # Fit model
  fit <- try(coxph(surv_formula, data = d, cluster = joint_id,
                   ties = "breslow", control = coxph.control(iter.max = 20)), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  
  # C-index
  c_index <- try(concordance(fit)$concordance, silent = TRUE)
  if (inherits(c_index, "try-error")) return(NULL)
  
  # Predict probability at max time
  lp <- predict(fit, type = "lp", newdata = d)
  basehaz <- basehaz(fit, centered = FALSE)
  pred_time <- max(d$tstop)
  h0_t <- basehaz[which.min(abs(basehaz$time - pred_time)), "hazard"]
  pred_prob <- 1 - exp(-h0_t * exp(lp))
  
  # Prepare output
  out_df <- data.frame(
    m = m,
    b = b,
    joint_id = d$joint_id,
    event = d[[event_name]],
    pred_prob = pred_prob,
    concordance = c_index
  )
  
  return(out_df)
}

