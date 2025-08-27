#' @export
step_fista_scad <- function(alpha, beta,
                            value_old,
                            opt,
                            step_size, lambda, t_old,
                            intercept, va, vb, x, y,
                            prob_fun = getProbRR.org,
                            a = 3.7) {  # SCAD hyperparameter
  
  if (!(opt %in% c("alpha","beta"))) {
    cli::cli_abort("Option 'opt' must be either 'alpha' or 'beta'.")
  }
  
  if (opt == "alpha") value <- alpha
  if (opt == "beta") value <- beta
  
  # Momentum update
  t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
  # damping_factor <- 0.8  # Reduce momentum effect
  # a_new <- damping_factor * (t_old - 1) / t_new
  a_new <- (t_old - 1) / t_new
  y_value_new <- value + a_new * (value - value_old)
  
  # Compute gradient
  if (opt == "alpha") gradient <- grad_nll(alpha = y_value_new, beta = beta,
                                           y = y, x = x, va = va, vb = vb,
                                           prob_fun, opt = "alpha"
  )
  
  if (opt == "beta") {
    gradient <- grad_nll(alpha = alpha, beta = y_value_new,
                         y = y, x = x, va = va, vb = vb,
                         prob_fun, opt = "beta"
    )
    
  }
  
  
  # Proximal gradient update with SCAD thresholding
  input <- y_value_new - step_size * gradient
  value_new <- scad_thres(input, lambda * step_size, a)
  
  # Maintain intercept term if specified
  if (intercept) value_new[1] <- input[1]
  
  # Return updated values in a structured list
  return(list(value_new = value_new, t_value = t_new, y_value = y_value_new))
}

#' @export
scad_thres <- function(entry, lambda, a) {
  # size safety of equivalence between SCAD and hard thresholding
  if (!(a >= 2)) cli::cli_abort("a < 2")
  
  # Vectorized Version
  e1 <- abs(entry) <= 2 * lambda & abs(entry) - lambda > 0
  
  e1.5 <- abs(entry) <= 2 * lambda & abs(entry) - lambda <= 0
  
  e2 <- abs(entry) > 2 * lambda & abs(entry) <= a * lambda
  
  entry[which(e1)] <-  sign(entry[which(e1)]) * (abs(entry[which(e1)]) - lambda)
  
  entry[which(e1.5)] <- 0
  
  entry[which(e2)] <- ((a - 1) * entry[which(e2)] - sign(entry[which(e2)]) * a * lambda) / (a - 2)
  
  return(entry)
}

#' @export
step_fista_adaptive_lasso <- function(alpha, beta,
                                      value_old,
                                      opt,
                                      step_size, lambda, t_old,
                                      intercept, va, vb, x, y,
                                      prob_fun = getProbRR.org,
                                      weights, gamma = 1) {  # Adaptive Lasso parameters
  
  if (!(opt %in% c("alpha","beta"))) {
    cli::cli_abort("Option 'opt' must be either 'alpha' or 'beta'.")
  }
  
  if (opt == "alpha") value <- alpha
  if (opt == "beta") value <- beta
  
  # Momentum update
  t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
  a_new <- (t_old - 1) / t_new
  y_value_new <- value + a_new * (value - value_old)
  
  # Compute gradient
  if (opt == "alpha") gradient <- grad_nll(y_value_new, beta,
                                           x, y, va, vb,
                                           prob_fun, opt = "alpha")
  
  
  if (opt == "beta") gradient <- grad_nll(alpha, y_value_new, 
                                          x, y, va, vb,
                                          prob_fun, opt = "beta")
  
  # Clean any NA gradients to prevent issues during computation
  gradient[is.na(gradient)] <- 0
  
  # Proximal gradient update with adaptive soft-thresholding
  input <- y_value_new - step_size * gradient
  value_new <- adaptive_soft_thres(input, lambda * step_size, weights)
  
  # Maintain intercept term if specified
  if (intercept) value_new[1] <- input[1]
  
  # Return updated values in a structured list
  return(list(value_new = value_new, t_value = t_new, y_value = y_value_new))
}

#' @export
# Adaptive Soft-Thresholding Function
# adaptive_soft_thres <- function(z, lambda, weights) {
#   sign_z <- sign(z)
#   abs_z <- abs(z)
#   # Adaptive soft-thresholding rule
#   result <- sign_z * pmax(abs_z - lambda * weights, 0)
#   return(result)
# }
adaptive_soft_thres <- function(entry, lambda, n) {
  # define thresholding multiplier
  s <- abs(entry) - (lambda^(n + 1)) * abs(entry)^(-n)
  
  # apply regularization
  s[s < 0] <- 0
  s[s > 0] <- sign(entry[s > 0]) * s[s > 0]
  
  # output regularized entry
  return(s)
}

