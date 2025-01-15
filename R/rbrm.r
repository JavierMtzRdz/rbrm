

#' FISTA Proximal Gradient Descent for Alpha and beta
#' 
#' Applies the FISTA algorithm for optimizing alpha and beta with 
#' proximal gradient descent.
#' 
#' @param alpha A numeric vector of alpha coefficients.
#' @param step_size A numeric value for the step size.
#' @param lambda A numeric value for the L1 regularization parameter.
#' @param t_old A numeric value for the previous t parameter in FISTA.
#' @param last_alpha A numeric vector representing the alpha coefficients from the previous iteration.
#' @return A list with the updated alpha, t, and y_alpha values.
#' @export
proximal.gd.fista <- function(alpha, beta, last_y,
                              opt,
                              step_size, lambda, t_old,
                                    intercept, va, vb, x, y,
                              prob_fun = getProbRR.org) {
  
  if (!(opt %in% c("alpha","beta"))) {
    cli::cli_abort("Option 'opt' must be either 'alpha' or 'beta'.")
  }
  
  if (opt == "alpha") {
    
    value <- alpha
    
    gradient <- numDeriv::grad(function(.x){nllh(.x, beta, va, vb, x, y,
                                                 prob_fun = prob_fun)}, 
                               last_y, method = "simple")
  }
  if (opt == "beta") {
    
    value <- beta
    
    gradient <- numDeriv::grad(function(.x) {nllh(alpha, .x, va, vb, x, y,
                                                  prob_fun = prob_fun)},
                               last_y, method = "simple")
  }
  
  
  # Clean any NA gradients to prevent issues during computation
  gradient[is.na(gradient)] <- 0
  
  # Proximal gradient update with soft-thresholding
  input <- last_y - step_size * gradient
  value_new <- soft_thres(input, lambda * step_size)
  
  # Maintain intercept term if specified
  if (intercept) value_new[1] <- input[1]
  
  
  # FISTA momentum update
  t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
  y_value_new <- value_new + (t_old - 1) / t_new * (value_new - value)
  
  # Return updated values in a structured list
  return(list(value_new = value_new, t_value = t_new, y_value = y_value_new))
}

#' Adaptive Step-Size FISTA Proximal Gradient Descent for Alpha and Beta
#' 
#' This function applies the FISTA (Fast Iterative Shrinkage-Thresholding Algorithm) 
#' for optimizing either `alpha` or `beta` parameters using proximal gradient descent 
#' with adaptive step size via backtracking line search.
#' 
#' The adaptive step size dynamically adjusts during each iteration to improve convergence 
#' by using backtracking line search. When a step size produces a sufficient objective decrease, 
#' it is increased in future iterations to accelerate convergence.
#' 
#' @param alpha A numeric vector of alpha coefficients.
#' @param beta A numeric vector of beta coefficients.
#' @param last_value A numeric vector representing the coefficients from the previous iteration.
#' @param opt A character string indicating the parameter to optimize; either "alpha" or "beta".
#' @param step_size A numeric value for the initial step size. This will be adjusted adaptively.
#' @param lambda A numeric value for the L1 regularization parameter.
#' @param t_old A numeric value for the previous t parameter in FISTA (used for momentum).
#' @param intercept A logical indicating whether to maintain the intercept term.
#' @param va, vb, x, y Data or model components passed to the objective function.
#' @param max_backtrack Integer, maximum number of backtracking iterations for adaptive step size.
#' @param backtrack_factor A numeric factor by which the step size is reduced when backtracking.
#' @param beta_increase A numeric factor by which the step size is increased if backtracking succeeds.
#' 
#' @return A list with the following components:
#'   \item{value_new}{A numeric vector of updated coefficients for the specified parameter (alpha or beta).}
#'   \item{t_value}{Updated t parameter for the next FISTA iteration.}
#'   \item{y_value}{A numeric vector of coefficients adjusted with FISTA momentum for the next iteration.}
#'   \item{step_size}{The adapted step size for use in subsequent iterations.}
#' 
#' @details
#' The algorithm performs a backtracking line search to find an adaptive step size that 
#' satisfies the Armijo condition (sufficient decrease condition), ensuring stable convergence.
#' If the condition is met, the step size is increased by a factor (`beta_increase`) 
#' to speed up future iterations. If not, the step size is decreased progressively 
#' (by `backtrack_factor`) until a satisfactory objective decrease is achieved.
#' 
#' @export
proximal.gd.asfista <- function(alpha, beta, last_y, 
                                opt, 
                                step_size, lambda, t_old,
                                intercept, va, vb, x, y,
                                max_backtrack = 10,  # Max backtracking iterations
                                backtrack_factor = 0.8,  # Step size reduction factor
                                beta_increase = 1.1,  # Factor to increase step size
                                prob_fun = getProbRR.org) {
  
  # Ensure the optimization option is valid
  if (!(opt %in% c("alpha", "beta"))) {
    cli::cli_abort("Option 'opt' must be either 'alpha' or 'beta'.")
  }
  
  # Initialize based on the optimization target ('alpha' or 'beta')
  if (opt == "alpha") {
    value <- alpha
    gradient <- numDeriv::grad(function(.x) {
      nllh(.x, beta, va, vb, x, y, prob_fun = prob_fun)
    }, last_y, method = "simple")
  } else if (opt == "beta") {
    value <- beta
    gradient <- numDeriv::grad(function(.x) {
      nllh(alpha, .x, va, vb, x, y, prob_fun = prob_fun)
    }, last_y, method = "simple")
  }
  
  prev_step_size <- step_size
  # Handle NA gradients
  gradient[is.na(gradient)] <- 0
  
  # Proximal gradient update
  input <- last_y - step_size * gradient
  
  value_new <- soft_thres(input, lambda * step_size)
  
  # Preserve intercept if specified
  if (intercept) {
    value_new[1] <- input[1]
  }
  
  # FISTA momentum update
  t_new <- (1 + sqrt(1 + (step_size/prev_step_size)*4 * t_old^2)) / 2
  
  y_value_new <- value_new + ((t_old - 1) / t_new )* (value_new - value)
  
  # Compute the initial objective value
  # obj_old <- penalized_nllh(alpha, beta, va, vb, x, y, 
  #                           lambda, intercept,
  #                           prob_fun = prob_fun)
  
  # Backtracking line search
  for (bt_iter in 1:max_backtrack) {
    if (opt == "alpha") {
      obj_new <- nllh(value_new, beta, va, vb, x, y,
                      # lambda, intercept,
                      prob_fun = prob_fun)
    } else if (opt == "beta") {
      
      obj_new <- nllh(alpha, value_new, va, vb, x, y,
                      # lambda, intercept,
                      prob_fun = prob_fun)
    }
    
    
    obj_old <- nllh(alpha, beta, va, vb, x, y,
                    # lambda, intercept,
                    prob_fun = prob_fun)
    
    
    # Check Armijo condition for sufficient decrease
    # purrr::walk(value_new, cli::cli_li)
    # cli::cli_ol(step_size)
    # cli::cli_ol(obj_new)
    # cli::cli_ol(obj_old - 0.5 * step_size * sum(gradient^2))
    
    if (obj_new <= obj_old - 0.5 * step_size * sum(gradient^2)) {
      # Successful line search: increase step size for the next iteration
      step_size <- step_size * backtrack_factor
      break
    } else {
      # Reduce step size and recompute the proximal gradient step
      step_size <- step_size * beta_increase
      # Proximal gradient update
      input <- last_y - step_size * gradient
      
      value_new <- soft_thres(input, lambda * step_size)
      
      # Preserve intercept if specified
      if (intercept) {
        value_new[1] <- input[1]
      }
      
      # FISTA momentum update
      t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
      
      y_value_new <- value_new + ((t_old - 1) / t_new )*(value_new - value)
    }
  }
  
  
  
  # Return updated values and adapted step size
  return(list(value_new = value_new, 
              t_value = t_new, 
              y_value = y_value_new, 
              step_size = step_size,
              gradient = gradient))
}


#' FISTA
#' 
#' Applies the FISTA algorithm for optimizing alpha and beta with 
#' proximal gradient descent.
#' 
#' @param alpha A numeric vector of alpha coefficients.
#' @param step_size A numeric value for the step size.
#' @param lambda A numeric value for the L1 regularization parameter.
#' @param t_old A numeric value for the previous t parameter in FISTA.
#' @param last_alpha A numeric vector representing the alpha coefficients from the previous iteration.
#' @return A list with the updated alpha, t, and y_alpha values.
#' @export
fista <- function(alpha, beta,
                  y_alpha, y_beta,
                  alpha.start, beta.start,
                  step_size_alpha, step_size_beta,
                  lambda, 
                  t_alpha, t_beta,
                  intercept, 
                  max.step, thres,
                  va, vb, x, y,
                  prob_fun = getProbRR.org,
                  opt_step = proximal.gd.fista){
  ## Optimization
  step <- 0
  for (iter in 1:max.step) {
    step <- step + 1
    # FISTA update for alpha
    last_alpha <- alpha
    res_alpha <- opt_step(alpha = alpha, beta = beta, 
                         last_y = y_alpha,
                         opt = "alpha",
                         step_size = step_size_alpha, 
                         lambda = lambda, 
                         t_old = t_alpha,
                         intercept = intercept,
                         va = va, vb = vb, x = x, y = y,
                         prob_fun = prob_fun)
    
    step_size_alpha <- ifelse(is.null(res_alpha$step_size), 
                              step_size_alpha, res_alpha$step_size)
    alpha <- res_alpha$value_new
    t_alpha <- res_alpha$t_value
    y_alpha <- res_alpha$y_value
    
    
    # FISTA update for beta
    last_beta <- beta
    res_beta <- opt_step(alpha = alpha, beta = beta, 
                        last_y = y_beta,
                        opt = "beta",
                        step_size = step_size_beta, 
                        lambda = lambda, t_old = t_beta,
                        intercept = intercept,
                        va = va, vb = vb, x = x, y = y,
                        prob_fun = prob_fun)
    
    
    step_size_beta <- ifelse(is.null(res_beta$step_size), 
                             step_size_beta, res_beta$step_size)
    beta <- res_beta$value_new
    t_beta <- res_beta$t_value
    y_beta <- res_beta$y_value
    
    # Update alpha and beta for the next iteration
    # beta <- beta_new
    
    grad_alpha <- numDeriv::grad(function(.x) {nllh(.x, beta, va, vb, x, y,
                                                    prob_fun = prob_fun)}, 
                                 alpha, method = "simple")
    
    grad_beta <- numDeriv::grad(function(.x) {nllh(alpha, .x, va, vb, x, y,
                                                   prob_fun = prob_fun)}, 
                                beta, method = "simple")
    
    
    # Stopping criteria based on gradient norm, relative change in objective, and relative change in parameters
    grad_norm_alpha <- max(abs(grad_alpha))
    grad_norm_beta <- max(abs(grad_beta))
    
    # rel_param_change_alpha <- max(abs(alpha_new - last_beta)) / max(1, max(abs(last_beta)))
    # rel_param_change_beta <- max(abs(beta_new - last_beta)) / max(1, max(abs(last_beta)))
    rel_param_change_alpha <- max(abs(alpha - last_alpha) / 
                                    pmax(1e-15, abs(last_alpha)))
    rel_param_change_beta <- max(abs(beta - last_beta) / 
                                   pmax(1e-15, abs(last_beta)))
    
    # pmax()
    # Break if all conditions are met
    
    if ((grad_norm_alpha < thres && grad_norm_beta < thres) #||
        # (rel_param_change_alpha < tol_ch && rel_param_change_beta < tol_ch)#||
        # (step_size_alpha < min_ss && step_size_beta < min_ss) 
        # rel_obj_change < tol && 
    ) {
      break
    }
    
  }
  
  return(list(alpha = alpha,
              beta = beta,
              step = step))
}


#' Greedy FISTA
#' 
#' Applies the Greedy FISTA algorithm for optimizing alpha and beta with 
#' proximal gradient descent.
#' 
#' @param alpha A numeric vector of alpha coefficients.
#' @param step_size A numeric value for the step size.
#' @param lambda A numeric value for the L1 regularization parameter.
#' @param t_old A numeric value for the previous t parameter in FISTA.
#' @param last_alpha A numeric vector representing the alpha coefficients from the previous iteration.
#' @return A list with the updated alpha, t, and y_alpha values.
#' @export
greedy_fista <- function(
    alpha, beta,
    y_alpha, y_beta,
    alpha.start, beta.start,
    step_size_alpha, step_size_beta,
    lambda, 
    t_alpha, t_beta,
    intercept, 
    max.step, thres,
    safeguard_factor = 10, gamma_shrink = 0.9,
    va, vb, x, y,
    prob_fun = getProbRR.org,
    opt_step = proximal.gd.fista
) {
  # Initialize variables
  step <- 0
  initial_diff <- NULL
  
  for (iter in 1:max.step) {
    step <- step + 1
    
    # FISTA update for alpha
    last_alpha <- alpha
    res_alpha <- opt_step(
      alpha = alpha, beta = beta, 
      last_y = y_alpha,
      opt = "alpha",
      step_size = step_size_alpha, 
      lambda = lambda, 
      t_old = t_alpha,
      intercept = intercept,
      va = va, vb = vb, x = x, y = y,
      prob_fun = prob_fun
    )
    step_size_alpha <- ifelse(is.null(res_alpha$step_size), step_size_alpha, res_alpha$step_size)
    alpha <- res_alpha$value_new
    t_alpha <- res_alpha$t_value
    y_alpha <- res_alpha$y_value
    
    # FISTA update for beta
    last_beta <- beta
    res_beta <- opt_step(
      alpha = alpha, beta = beta, 
      last_y = y_beta,
      opt = "beta",
      step_size = step_size_beta, 
      lambda = lambda, 
      t_old = t_beta,
      intercept = intercept,
      va = va, vb = vb, x = x, y = y,
      prob_fun = prob_fun
    )
    step_size_beta <- ifelse(is.null(res_beta$step_size), step_size_beta, res_beta$step_size)
    beta <- res_beta$value_new
    t_beta <- res_beta$t_value
    y_beta <- res_beta$y_value
    
    # Safeguard mechanism: Check if update diverges
    current_diff <- sqrt(sum((alpha - last_alpha)^2 + (beta - last_beta)^2))
    if (is.null(initial_diff)) initial_diff <- current_diff
    if (current_diff > safeguard_factor * initial_diff) {
      step_size_alpha <- max(step_size_alpha * gamma_shrink, 1e-10)
      step_size_beta <- max(step_size_beta * gamma_shrink, 1e-10)
    }
    
    # Restarting mechanism
    alpha_diff <- alpha - last_alpha
    beta_diff <- beta - last_beta
    if ((y_alpha - alpha) %*% alpha_diff >= 0) y_alpha <- alpha
    if ((y_beta - beta) %*% beta_diff >= 0) y_beta <- beta
    
    # Gradient calculations for stopping criteria
    grad_alpha <- numDeriv::grad(function(.x) {
      nllh(.x, beta, va, vb, x, y, prob_fun = prob_fun)
    }, alpha, method = "simple")
    
    grad_beta <- numDeriv::grad(function(.x) {
      nllh(alpha, .x, va, vb, x, y, prob_fun = prob_fun)
    }, beta, method = "simple")
    
    grad_norm_alpha <- max(abs(grad_alpha))
    grad_norm_beta <- max(abs(grad_beta))
    
    # Relative change for alpha and beta
    rel_param_change_alpha <- max(abs(alpha - last_alpha) / pmax(1e-15, abs(last_alpha)))
    rel_param_change_beta <- max(abs(beta - last_beta) / pmax(1e-15, abs(last_beta)))
    
    # Convergence check
    if (grad_norm_alpha < thres && grad_norm_beta < thres) {
      break
    }
  }
  
  return(list(alpha = alpha, beta = beta, step = step))
}

#' Regularized Binary Regression Model experimental (RBRM)
#'
#' Performs a regularized binary regression model (RBRM) using FISTA proximal gradient descent. 
#' The model penalizes the negative log-likelihood and includes optional early stopping.
#'
#' @param va A matrix of independent variables (without an intercept) for alpha.
#' @param vb A matrix of independent variables (without an intercept) for beta. If NULL, va is used.
#' @param y A vector of binary outcomes (0/1).
#' @param x A vector indicating the group assignment (1 or 0) for each observation.
#' @param alpha.start Initial values for alpha coefficients. Defaults to a zero vector.
#' @param beta.start Initial values for beta coefficients. Defaults to a vector of 0.01.
#' @param max.step Maximum number of optimization steps. Default is 3000.
#' @param thres Threshold for convergence. Default is 1e-04.
#' @param lambda Regularization parameter for L1 penalty. Default is 0 (no regularization).
#' @param lr.alpha Learning rate for alpha parameters. Default is 0.06.
#' @param lr.beta Learning rate for beta parameters. Default is 0.02.
#' @param intercept Logical. If TRUE, the intercept is included and not penalized. Default is TRUE.
#' @param early_stopping_rounds The number of rounds without improvement before early stopping occurs.
#'
#' @return A list with the following components:
#' \item{point.est}{The estimated alpha and beta coefficients.}
#' \item{convergence}{Logical indicating whether convergence was achieved within the max steps.}
#' \item{value}{The penalized negative log-likelihood value at convergence.}
#' \item{step}{The number of steps taken to converge.}
#'
#' @examples
#' # Example usage:
#' #va <- matrix(c(1, 1, 1, 1, 0, 0, 0, 0), ncol = 2)
#' #vb <- matrix(c(1, 1, 0, 0, 1, 1, 0, 0), ncol = 2)
#' #y <- c(1, 0, 1, 0)
#' #x <- c(1, 1, 0, 0)
#' #result <- rbrm(va, vb, y, x)
#'
#' @export

rbrm.experimental <- function(va, vb, x, y,
                              alpha.start = NULL, beta.start = NULL,
                 max.step = 3000, thres = 1e-04, lambda = 0,
                 lr.alpha = 0.1, lr.beta = 0.1,
                 intercept = TRUE, early_stopping_rounds = 10,
                 min_ss = 1e-13, tol_ch = 1e-5,
                 prob_fun = getProbRR.org,
                 opt_fun = fista) {
  
  # va <- v; vb <- v; alpha.start = NULL; beta.start = NULL;
  # max.step = 3000; thres = 1e-04; lambda = 0;
  # lr.alpha = 0.06; lr.beta = 0.02;
  # intercept = TRUE; early_stopping_rounds = 10
  tictoc::tic("Total time")
  
  if (is.null(vb)) {
    vb <- va
  }
  pa <- dim(va)[2]
  pb <- dim(vb)[2]
  
  # sanity check for the intercept term in va, vb
  if (all(va[, 1] == 1) & all(vb[, 1] == 1)) {
    intercept <- TRUE
  }
  
  ## starting values for parameter optimization
  if (is.null(alpha.start)) {
    alpha.start <- c(rep(0, pa))
  }
  if (is.null(beta.start)) {
    beta.start <- c(rep(0.01, pb))
  }
  
  ## Optimization
  alpha <- y_alpha <- alpha.start
  beta <- y_beta <- beta.start
  t_alpha <- t_beta <- 1
  step_size_alpha <- lr.alpha
  step_size_beta <- lr.beta
  
  opt_result <- opt_fun(alpha, beta,
                        y_alpha, y_beta,
                        alpha.start, beta.start,
                        step_size_alpha, step_size_beta,
                        lambda, 
                        t_alpha, t_beta,
                        intercept, 
                        max.step, thres,
                        va, vb, x, y,
                        prob_fun)
  
  step <- opt_result$step
  alpha <- opt_result$alpha
  beta <- opt_result$beta
  
  time <- tictoc::toc(quiet = TRUE)
  
  opt <- list(
    point.est = c(alpha, beta), 
    convergence = (step < max.step),
    value = penalized.neg.log.likelihood(c(alpha, beta), 
                                         pa, pb, intercept, lambda,
                                         alpha, beta, va, vb, x, y,
                                         prob_fun = prob_fun),
    step = step,
    time = round(time$toc - time$tic, 4)
  )
  
  return(structure(opt, class = c("rbrm")))
}

#' Regularized Binary Regression Model experimental (RBRM) 2
#'
#' Performs a regularized binary regression model (RBRM) using FISTA proximal gradient descent. 
#' The model penalizes the negative log-likelihood and includes optional early stopping.
#'
#' @param va A matrix of independent variables (without an intercept) for alpha.
#' @param vb A matrix of independent variables (without an intercept) for beta. If NULL, va is used.
#' @param y A vector of binary outcomes (0/1).
#' @param x A vector indicating the group assignment (1 or 0) for each observation.
#' @param alpha.start Initial values for alpha coefficients. Defaults to a zero vector.
#' @param beta.start Initial values for beta coefficients. Defaults to a vector of 0.01.
#' @param max.step Maximum number of optimization steps. Default is 3000.
#' @param thres Threshold for convergence. Default is 1e-04.
#' @param lambda Regularization parameter for L1 penalty. Default is 0 (no regularization).
#' @param lr.alpha Learning rate for alpha parameters. Default is 0.06.
#' @param lr.beta Learning rate for beta parameters. Default is 0.02.
#' @param intercept Logical. If TRUE, the intercept is included and not penalized. Default is TRUE.
#' @param early_stopping_rounds The number of rounds without improvement before early stopping occurs.
#'
#' @return A list with the following components:
#' \item{point.est}{The estimated alpha and beta coefficients.}
#' \item{convergence}{Logical indicating whether convergence was achieved within the max steps.}
#' \item{value}{The penalized negative log-likelihood value at convergence.}
#' \item{step}{The number of steps taken to converge.}
#'
#' @examples
#' # Example usage:
#' #va <- matrix(c(1, 1, 1, 1, 0, 0, 0, 0), ncol = 2)
#' #vb <- matrix(c(1, 1, 0, 0, 1, 1, 0, 0), ncol = 2)
#' #y <- c(1, 0, 1, 0)
#' #x <- c(1, 1, 0, 0)
#' #result <- rbrm(va, vb, y, x)
#'
#' @export

rbrm.experimental2 <- function(va, vb, x, y,
                              alpha.start = NULL, beta.start = NULL,
                              max.step = 3000, thres = 1e-04, lambda = 0,
                              lr.alpha = 0.1, lr.beta = 0.1,
                              intercept = TRUE, early_stopping_rounds = 10,
                              min_ss = 1e-13, tol_ch = 1e-3,
                              prob_fun = getProbRR.alt,
                              opt_fun = proximal.gd.fista) {
  

  # va <- v; vb <- v; alpha.start = NULL; beta.start = NULL;
  # max.step = 3000; thres = 1e-04; lambda = 0;
  # lr.alpha = 0.06; lr.beta = 0.02;
  # intercept = TRUE; early_stopping_rounds = 10
  
  tictoc::tic("Total time")
  
  if (is.null(vb)) {
    vb <- va
  }
  pa <- dim(va)[2]
  pb <- dim(vb)[2]
  
  # sanity check for the intercept term in va, vb
  if (all(va[, 1] == 1) & all(vb[, 1] == 1)) {
    intercept <- TRUE
  }
  
  ## starting values for parameter optimization
  if (is.null(alpha.start)) {
    alpha.start <- c(rep(0, pa))
  }
  if (is.null(beta.start)) {
    beta.start <- c(rep(0.01, pb))
  }
  
  ## Optimization
  alpha <- alpha.start
  beta <- beta.start
  last_alpha <- last_beta <- alpha
  t_alpha <- t_beta <- 1
  y_alpha <- y_beta <- alpha.start
  step_size_alpha <- lr.alpha
  step_size_beta <- lr.beta
  step <- 0
  for (iter in 1:max.step) {
    # print step size every 100 steps
    # if (step %% 100 == 0) {
    #     print(paste0("step ", step))
    #     # print("alpha step size")
    #     # print(t_alpha)
    #     # print("beta step size")
    #     # print(t_beta)
    # }
    step <- step + 1
    # FISTA update for alpha
    last_alpha <- alpha
    res_alpha <- opt_fun(alpha = y_alpha,
                         beta = beta,
                         last_value = last_alpha,
                         opt = "alpha", 
                         step_size = step_size_alpha,
                         lambda = lambda, 
                         t_old = t_alpha, 
                         intercept = intercept,
                         va = va, vb = vb, x = x, y = y,
                         prob_fun = prob_fun)
    
    step_size_alpha <- ifelse(is.null(res_alpha$step_size), 
                              step_size_alpha, res_alpha$step_size)
    
    alpha_new <- res_alpha$value_new
    
    t_alpha <- res_alpha$t_value
    y_alpha <- res_alpha$y_value
    
    # FISTA update for beta
    last_beta <- beta
    res_beta <- opt_fun(alpha = alpha,
                        beta = y_beta,
                        last_value = last_beta,
                        opt = "beta", 
                        step_size = step_size_beta,
                        lambda = lambda, 
                        t_old = t_beta, 
                        intercept = intercept,
                        va = va, vb = vb, x = x, y = y,
                        prob_fun = prob_fun)
    
    step_size_beta <- ifelse(is.null(res_beta$step_size), 
                             step_size_beta, res_beta$step_size)
    
    beta_new <- res_beta$value_new
    
    t_beta <- res_beta$t_value
    y_beta <- res_beta$y_value
    
    grad_alpha <- numDeriv::grad(function(.x) {nllh(.x, beta, va, vb, x, y,
                                                    prob_fun = prob_fun)}, 
                                 y_alpha, method = "simple")
    grad_beta <- numDeriv::grad(function(.x) {nllh(alpha, .x, va, vb, x, y,
                                                   prob_fun = prob_fun)}, 
                                y_beta, method = "simple")
    
    if ((max(abs(grad_alpha)) < thres) & (max(abs(grad_beta)) < thres)) {
      break
    }
    
    # Update alpha and beta for the next iteration
    alpha <- alpha_new
    beta <- beta_new
  }
  
  time <- tictoc::toc(quiet = TRUE)
  
  opt <- list(
    point.est = c(alpha, beta), convergence = (step < max.step),
    value = penalized.neg.log.likelihood(c(alpha, beta), 
                                         pa, pb, intercept, lambda,
                                         alpha, beta, va, vb, x, y,
                                         prob_fun = prob_fun),
    step = step,
    time = round(time$toc - time$tic, 4)
  )
  
  return(structure(opt, class = c("rbrm")))
}

