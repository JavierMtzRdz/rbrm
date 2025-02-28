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
  damping_factor <- 0.8  # Reduce momentum effect
  a_new <- damping_factor * (t_old - 1) / t_new
  # a_new <- (t_old - 1) / t_new
  y_value_new <- value + a_new * (value - value_old)
  
  # Compute gradient
  if (opt == "alpha") {
    gradient <- numDeriv::grad(function(.x) {
      nllh(.x, beta, va, vb, x, y, prob_fun = prob_fun)
    }, y_value_new, method = "simple")
  }
  
  if (opt == "beta") {
    gradient <- numDeriv::grad(function(.x) {
      nllh(alpha, .x, va, vb, x, y, prob_fun = prob_fun)
    }, y_value_new, method = "simple")
  }
  
  
  # Clean any NA gradients to prevent issues during computation
  if (any(is.na(gradient))) {
    cli::cli_alert_danger("NaN in gradient, replacing with 0.")
    print(gradient)
    gradient[is.na(gradient)] <- 0
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
  lambda <- lambda*10
  
  # Vectorized Version
  e1 <- abs(entry) <= 2 * lambda
  e2 <- abs(entry) > 2 * lambda & abs(entry) <= a * lambda
  
  entry[e1] <- ifelse(
    abs(entry[e1]) - lambda > 0, sign(entry[e1]) * (abs(entry[e1]) - lambda), 0
  )
  
  entry[e2] <- ((a - 1) * entry[e2] - sign(entry[e2]) * a * lambda) / (a - 2)
  
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
  if (opt == "alpha") {
    gradient <- numDeriv::grad(function(.x) {
      nllh(.x, beta, va, vb, x, y, prob_fun = prob_fun)
    }, y_value_new, method = "simple")
  }
  
  if (opt == "beta") {
    gradient <- numDeriv::grad(function(.x) {
      nllh(alpha, .x, va, vb, x, y, prob_fun = prob_fun)
    }, y_value_new, method = "simple")
  }
  
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
adaptive_soft_thres <- function(z, lambda, weights) {
  sign_z <- sign(z)
  abs_z <- abs(z)
  # Adaptive soft-thresholding rule
  result <- sign_z * pmax(abs_z - lambda * weights, 0)
  return(result)
}

