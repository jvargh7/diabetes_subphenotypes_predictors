
############  C-INDEX ############

c_index_function <- function(data, indices, event_name) {
  unique_ids <- unique(data$joint_id)
  sampled_ids <- unique_ids[indices]
  d <- data[data$joint_id %in% sampled_ids, ]
  
  if (length(unique(d$race)) == 1 || length(unique(d$female)) == 1) return(NA)
  
  formula_str <- paste0("Surv(tstart, tstop, ", event_name, ") ~ study + female + race + earliest_age + bmi + hba1c + ",
                        "homa2b_scaled + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention")
  surv_formula <- as.formula(formula_str)
  
  fit <- try(coxph(surv_formula, data = d, cluster = joint_id,
                   ties = "breslow", control = coxph.control(iter.max = 20)), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA)
  
  return(concordance(fit)$concordance)
}


########### SENSITIVITY & SPECIFICITY ########### 

sens_spec_function <- function(data, indices, event_name) {
  unique_ids <- unique(data$joint_id)
  sampled_ids <- unique_ids[indices]
  d <- data[data$joint_id %in% sampled_ids, ]
  
  if (length(unique(d[[event_name]])) < 2) return(c(NA, NA, NA))  # ensure variability
  
  formula_str <- paste0("Surv(tstart, tstop, ", event_name, ") ~ study + female + race + earliest_age + bmi + hba1c + ",
                        "homa2b_scaled + homa2ir + ldlc_scaled + sbp_scaled + egfr_ckdepi_2021_scaled + dpp_intervention")
  surv_formula <- as.formula(formula_str)
  
  fit <- try(coxph(surv_formula, data = d, ties = "breslow", control = coxph.control(iter.max = 20)), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(NA, NA, NA))
  
  lp <- predict(fit, type = "lp", newdata = d)
  basehaz <- basehaz(fit, centered = FALSE)
  
  pred_time <- max(d$tstop)
  h0_t <- basehaz[which.min(abs(basehaz$time - pred_time)), "hazard"]
  
  d$pred_prob <- 1 - exp(-h0_t * exp(lp))
  d$predicted_class <- ifelse(d$pred_prob >= 0.5, 1, 0)
  
  actual <- d[[event_name]]
  predicted <- d$predicted_class
  
  TP <- sum(actual == 1 & predicted == 1)
  TN <- sum(actual == 0 & predicted == 0)
  FP <- sum(actual == 0 & predicted == 1)
  FN <- sum(actual == 1 & predicted == 0)
  
  sensitivity <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  specificity <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
  precision   <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
  recall      <- sensitivity
  f1 <- ifelse(is.na(precision) || is.na(recall) || (precision + recall) == 0,
               NA,
               2 * precision * recall / (precision + recall))
  
  return(c(sensitivity, specificity, f1))
}


########### CALIBRATION SLOPE ########### 

