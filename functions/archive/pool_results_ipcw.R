# pool results following Rubin's rules

pool_results <- function(model) {
  
  D = length(model)
  
  # Bind rows from all tidied models
  results <- bind_rows(lapply(model, broom::tidy))
  
  # Prepare output dataframe with calculations
  output <- data.frame(Estimate = results$estimate,
                       SE = results$robust.se,
                       term = results$term
                       # dfcom = dfcom_coxph  # Uncomment and define dfcom_coxph if needed
  ) %>% 
    mutate(W_d = SE^2) %>% 
    group_by(term) %>% 
    mutate(B_D = var(Estimate)) %>% 
    dplyr::summarize(B_D = mean(B_D), # B: Variance of estimates (between imputation variance)
                     W_D = mean(W_d), # \bar{V}: average of V_d over D imputed datasets
                     theta_D = mean(Estimate) # \bar{\theta}: mean of estimates,
                     # dfcom = mean(dfcom)  # Uncomment if dfcom_coxph is defined
    ) %>% 
    mutate(T_D = W_D + (1 + 1/D)*B_D, # Var(\theta|Y_{0}) ~ improved approximation of posterior variance [\bar{V} + B] 
           gamma_D = (1 + 1/D)*(B_D/T_D), # \hat{\gamma}_D = between imputation : total variance --> fraction of missing information
           nu = (D-1)*((1+ (1/(D+1))*(W_D/B_D))^2), # degrees of freedom of t-distribution
           nu2 = (D-1)/(gamma_D)^2 # equivalent to mice:::pool.fitlist >> mice:::barnard.rubin()'s dfold; (D/(D+1)) and not (1/(D+1))
           # nu_improved = mice:::barnard.rubin(D,B_D,T_D,dfcom = dfcom)  # Uncomment if applicable
    ) %>% 
    mutate(L = theta_D + qt(p = 0.025, df = nu2) * sqrt(T_D),
           U = theta_D + qt(p = 0.975, df = nu2) * sqrt(T_D),
           sqrt_T_D = sqrt(T_D)
    ) %>% 
    mutate(HR = paste0(round(exp(theta_D), 2), " \t (",
                       round(exp(L), 2), ", ",
                       round(exp(U), 2), ")"),
           estimate = exp(theta_D),
           lci = exp(L),
           uci = exp(U)
    ) %>% 
    rename(iv = term)
  
  # Return the final output
  return(output)
}
