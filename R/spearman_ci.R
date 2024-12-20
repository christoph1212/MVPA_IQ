spearman_ci <- function(x, y, cor_test) {
  
  # Spearman Rho 
  rho <- cor_test$estimate
  
  # Sample Size
  n <- sum(!is.na(x) & !is.na(y))
  
  # Z-Transformation
  z <- 0.5 * log((1 + rho) / (1 - rho))
  
  # Standard Error
  se <- 1 / sqrt(n - 3)
  
  # Z-Value for CI
  alpha = 0.05
  z_alpha <- qnorm(1 - alpha)
  
  # CI
  ci_lower_z <- z - z_alpha * se
  ci_upper_z <- z + z_alpha * se
  
  ci_lower <- (exp(2 * ci_lower_z) - 1) / (exp(2 * ci_lower_z) + 1)
  ci_upper <- (exp(2 * ci_upper_z) - 1) / (exp(2 * ci_upper_z) + 1)
  
  output <- sprintf("r(%d) = %.2f, 95%% CI [%.2f, %.2f]", n, rho, ci_lower, ci_upper)
  
  return(output)
  
}