
evaluate_metrics_with_ci <- function(df, event_label, cutoffs = c(0.10, 0.25, 0.5, 0.75, 0.9)) {
  results <- list()
  
  for (cut in cutoffs) {
    df <- df %>% 
      mutate(pred_class = ifelse(pred_prob >= cut, 1, 0))
    
    metrics_by_boot <- df %>%
      group_by(m, b) %>%
      summarise(
        sensitivity = {
          TP <- sum(event == 1 & pred_class == 1)
          FN <- sum(event == 1 & pred_class == 0)
          ifelse((TP + FN) == 0, NA, TP / (TP + FN))
        },
        specificity = {
          TN <- sum(event == 0 & pred_class == 0)
          FP <- sum(event == 0 & pred_class == 1)
          ifelse((TN + FP) == 0, NA, TN / (TN + FP))
        },
        precision = {
          TP <- sum(event == 1 & pred_class == 1)
          FP <- sum(event == 0 & pred_class == 1)
          ifelse((TP + FP) == 0, NA, TP / (TP + FP))
        },
        recall = sensitivity,
        f1 = {
          p <- precision
          r <- recall
          ifelse(is.na(p) || is.na(r) || (p + r) == 0, NA, 2 * p * r / (p + r))
        },
        c_index = if (cut == 0.5) mean(unique(concordance), na.rm = TRUE) else NA,
        .groups = "drop"
      )
    
    # Summarize mean + 95% CI
    result_row <- data.frame(
      event = event_label,
      cutoff = cut,
      c_index = round(mean(metrics_by_boot$c_index, na.rm = TRUE), 3),
      c_index_lower = round(quantile(metrics_by_boot$c_index, 0.025, na.rm = TRUE), 3),
      c_index_upper = round(quantile(metrics_by_boot$c_index, 0.975, na.rm = TRUE), 3),
      sensitivity = round(mean(metrics_by_boot$sensitivity, na.rm = TRUE), 3),
      sens_lower = round(quantile(metrics_by_boot$sensitivity, 0.025, na.rm = TRUE), 3),
      sens_upper = round(quantile(metrics_by_boot$sensitivity, 0.975, na.rm = TRUE), 3),
      specificity = round(mean(metrics_by_boot$specificity, na.rm = TRUE), 3),
      spec_lower = round(quantile(metrics_by_boot$specificity, 0.025, na.rm = TRUE), 3),
      spec_upper = round(quantile(metrics_by_boot$specificity, 0.975, na.rm = TRUE), 3),
      f1 = round(mean(metrics_by_boot$f1, na.rm = TRUE), 3),
      f1_lower = round(quantile(metrics_by_boot$f1, 0.025, na.rm = TRUE), 3),
      f1_upper = round(quantile(metrics_by_boot$f1, 0.975, na.rm = TRUE), 3)
    )
    
    results[[as.character(cut)]] <- result_row
  }
  
  return(bind_rows(results))
}
