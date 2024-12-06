#' Soft Thresholding Operator
#' 
#' Applies the soft thresholding operator to the input.
#' 
#' @param x A numeric vector to apply the soft threshold.
#' @param lambda A positive numeric value representing the threshold parameter.
#' @return A numeric vector where the soft thresholding has been applied.
#' @examples
#' #soft_thres(c(3, -1.5, 0.2), 0.5)
#' @export
soft_thres <- function(x, lambda) {
  sx <- abs(x) - lambda
  sx[sx < 0] <- 0
  return(sx * sign(x))
}



#' Negative Log-Likelihood for Alpha
#' 
#' Computes the negative log-likelihood for the alpha coefficients.
#' 
#' @param alpha Coefficients.
#' @return The negative log-likelihood for the given alpha.
#' @examples
#' #alpha <- c(1, 2, 3)
#' #nllh.alpha(alpha)
#' @export
nllh <- function(alpha, beta, va, vb, x, y,
                 prob_fun = brm::getProbRR) {
  logrr <- (va %*% alpha)
  logop <- (vb %*% beta)
  
  ps <- prob_fun(logrr, logop)
  
  p0 <- ps[, 1]
  p1 <- ps[, 2]
  
  # Clipping probabilities
  p0 <- pmin(pmax(p0, 1e-15), 1 - 1e-15)
  p1 <- pmin(pmax(p1, 1e-15), 1 - 1e-15)
  
  return(-sum((1 - y[x == 0]) * log1p(-p0[x == 0]) + 
                (y[x == 0]) * log(p0[x == 0])) -
           sum((1 - y[x == 1]) * log1p(-p1[x == 1]) + 
                 (y[x == 1]) * log(p1[x == 1])))
}


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
proximal.gd.fista <- function(alpha, beta, last_value,
                              opt,
                              step_size, lambda, t_old,
                                    intercept, va, vb, x, y,
                              prob_fun = brm::getProbRR) {
  
  if (!(opt %in% c("alpha","beta"))) {
    cli::cli_abort("Option 'opt' must be either 'alpha' or 'beta'.")
  }
  
  if (opt == "alpha") {
    
    value <- alpha
    
    gradient <- numDeriv::grad(function(.x){nllh(.x, beta, va, vb, x, y,
                                                 prob_fun = prob_fun)}, 
                               alpha, method = "simple")
  }
  if (opt == "beta") {
    
    value <- beta
    
    gradient <- numDeriv::grad(function(.x) {nllh(alpha, .x, va, vb, x, y,
                                                  prob_fun = prob_fun)},
                               beta, method = "simple")
  }
  
  
  # Clean any NA gradients to prevent issues during computation
  if (any(is.na(gradient))) gradient[is.na(gradient)] <- 0
  
  # Proximal gradient update with soft-thresholding
  input <- value - step_size * gradient
  value_new <- soft_thres(input, lambda * step_size)
  
  # Maintain intercept term if specified
  if (intercept) value_new[1] <- input[1]
  
  
  # FISTA momentum update
  t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
  y_value_new <- value_new + (t_old - 1) / t_new * (value_new - last_value)
  
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
proximal.gd.asfista <- function(alpha, beta, last_value,
                                opt,
                                step_size, lambda, t_old,
                                intercept, va, vb, x, y,
                                max_backtrack = 15,  # Max backtracking iterations
                                backtrack_factor = 0.5,  # Step size reduction factor
                                beta_increase = 1.2,
                                prob_fun = brm::getProbRR) {  # Factor to increase step size
  
  # Set initial values based on 'opt' parameter
  if (opt == "alpha") {
    value <- alpha
    gradient <- numDeriv::grad(function(.x) {nllh(.x, beta, va, vb, x, y,
                                                  prob_fun = prob_fun)}, 
                               alpha, method = "simple")
  } else if (opt == "beta") {
    value <- beta
    gradient <- numDeriv::grad(function(.x) {nllh(alpha, .x, va, vb, x, y,
                                                  prob_fun = prob_fun)}, 
                               beta, method = "simple")
  } else {
    stop("Option 'opt' must be either 'alpha' or 'beta'.")
  }
  
  # Clean up any NA values in the gradient
  gradient[is.na(gradient)] <- 0
  
  # Adaptive step-size adjustment with backtracking line search
  input <- value - step_size * gradient
  value_new <- soft_thres(input, lambda * step_size)
  
  # Objective function value after update (for backtracking line search)
  obj_old <- nllh(alpha, beta, va, vb, x, y,
                  prob_fun = prob_fun)  # Calculate at current parameters
  
  
  
  
  # Backtracking line search to adapt step size
  for (bt_iter in 1:max_backtrack) {
    if (opt == "alpha") {
      obj_new <- nllh(value_new, beta, va, vb, x, y,
                      prob_fun = prob_fun)
    } else {
      obj_new <- nllh(alpha, value_new, va, vb, x, y,
                      prob_fun = prob_fun)
    }
    
    # Check Armijo condition for sufficient decrease
    # logrr <- (va %*% alpha)
    # logop <- (vb %*% beta)
    # 
    # ps <- prob_fun(logrr, logop)
    # num <- length(ps)/2
    # 
    # cli::cli_alert("p0: {ps[1:num]}")
    # cli::cli_alert("p1: {ps[(num+1):(num*2)]}")
    
    
    if (obj_new <= obj_old - 0.5 * step_size * sum(gradient^2)) {
      # If line search was successful, increase step size for future iterations
      step_size <- step_size * beta_increase
      break
    } else {
      # Reduce step size and try again
      step_size <- step_size * backtrack_factor
      input <- value - step_size * gradient
      value_new <- soft_thres(input, lambda * step_size)
    }
  }
  
  
  # FISTA momentum update
  t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
  y_value_new <- value_new + (t_old - 1) / t_new * (value_new - last_value)
  
  # Return updated values and adapted step size
  return(list(value_new = value_new, 
              t_value = t_new, 
              y_value = y_value_new, 
              step_size = step_size,
              gradient = gradient))
}


#' FISTA Proximal Gradient Descent for Alpha
#' 
#' Applies the FISTA algorithm for optimizing alpha with proximal gradient descent.
#' 
#' @param alpha A numeric vector of alpha coefficients.
#' @param step_size A numeric value for the step size.
#' @param lambda A numeric value for the L1 regularization parameter.
#' @param t_old A numeric value for the previous t parameter in FISTA.
#' @param last_alpha A numeric vector representing the alpha coefficients from the previous iteration.
#' @return A list with the updated alpha, t, and y_alpha values.
#' @export
proximal.gd.alpha.fista <- function(alpha, step_size, lambda, t_old, last_alpha,
                                    intercept,
                                    beta, va, vb, x, y,
                                    prob_fun = brm::getProbRR) {
  gradient <- numDeriv::grad(function(.x){nllh(.x, beta, va, vb, x, y,
                                               prob_fun = prob_fun)}, alpha, method = "simple")
  gradient[is.na(gradient)] <- 0
  input <- alpha - step_size * gradient
  alpha_new <- soft_thres(input, lambda * step_size)
  if (intercept == TRUE) {
    alpha_new[1] <- input[1]
  }
  t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
  y_alpha_new <- alpha_new + (t_old - 1) / t_new * (alpha_new - last_alpha)
  return(list(alpha_new = alpha_new, t_alpha = t_new, y_alpha = y_alpha_new))
}




#' FISTA Proximal Gradient Descent for Beta
#' 
#' Applies the FISTA algorithm for optimizing beta with proximal gradient descent.
#' 
#' @param beta A numeric vector of beta coefficients.
#' @param step_size A numeric value for the step size.
#' @param lambda A numeric value for the L1 regularization parameter.
#' @param t_old A numeric value for the previous t parameter in FISTA.
#' @param last_beta A numeric vector representing the beta coefficients from the previous iteration.
#' @return A list with the updated beta, t, and y_beta values.
#' @export
proximal.gd.beta.fista <- function(beta, step_size, lambda, t_old, last_beta,
                                   intercept,
                                   alpha, va, vb, x, y,
                                   prob_fun = brm::getProbRR) {
  gradient <- numDeriv::grad(function(.x) {nllh(alpha, .x, va, vb, x, y,
                                                prob_fun = prob_fun)},
                             beta, method = "simple")
  gradient[is.na(gradient)] <- 0
  input <- beta - step_size * gradient
  beta_new <- soft_thres(input, lambda * step_size)
  if (intercept == TRUE) {
    beta_new[1] <- input[1]
  }
  t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
  y_beta_new <- beta_new + (t_old - 1) / t_new * (beta_new - last_beta)
  return(list(beta_new = beta_new, t_beta = t_new, y_beta = y_beta_new))
}





#' Penalized Negative Log-Likelihood
#' 
#' Computes the penalized negative log-likelihood for alpha and beta.
#' 
#' @param pars A numeric vector of the concatenated alpha and beta coefficients.
#' @return The penalized negative log-likelihood.
#' @examples
#' #pars <- c(alpha = 1, beta = 2)
#' #penalized.neg.log.likelihood(pars)
#' @export
penalized.neg.log.likelihood <- function(pars, pa, pb, intercept, lambda,
                                         alpha, beta, va, vb, x, y,
                                         prob_fun = brm::getProbRR) {
  alpha <- pars[1:pa]
  beta <- pars[(pa + 1):(pa + pb)]
  
  unpenalized.nllh <- nllh(alpha, beta, va, vb, x, y,
                           prob_fun = prob_fun)
  
  # Applying the penalty term
  if (intercept == TRUE) {
    penalty <- lambda * (sum(abs(alpha[-1])) + sum(abs(beta[-1]))) # Exclude intercept
  } else {
    penalty <- lambda * (sum(abs(alpha)) + sum(abs(beta)))
  }
  return(unpenalized.nllh + penalty)
}

#' Penalized Negative Log-Likelihood
#' 
#' Computes the penalized negative log-likelihood for alpha and beta.
#' 
#' @param pars A numeric vector of the concatenated alpha and beta coefficients.
#' @return The penalized negative log-likelihood.
#' @examples
#' #pars <- c(alpha = 1, beta = 2)
#' #penalized.neg.log.likelihood(pars)
#' @export
penalized_nllh <- function(alpha, beta, va, vb, x, y,
                           lambda, intercept,
                           prob_fun = brm::getProbRR) {
  logrr <- (va %*% alpha)
  logop <- (vb %*% beta)
  
  ps <- prob_fun(logrr, logop)
  
  p0 <- ps[, 1]
  p1 <- ps[, 2]
  
  # Clipping probabilities
  p0 <- pmin(pmax(p0, 1e-15), 1 - 1e-15)
  p1 <- pmin(pmax(p1, 1e-15), 1 - 1e-15)
  
  unpenalized.nllh <- -sum((1 - y[x == 0]) * log1p(-p0[x == 0]) +
                             (y[x == 0]) * log(p0[x == 0])) - 
    sum((1 - y[x == 1]) * log1p(-p1[x == 1]) + 
          (y[x == 1]) * log(p1[x == 1]))
  
  # Applying the penalty term
  if (intercept == TRUE) {
    penalty <- lambda * (sum(abs(alpha[-1])) + sum(abs(beta[-1]))) # Exclude intercept
  } else {
    penalty <- lambda * (sum(abs(alpha)) + sum(abs(beta)))
  }
  return(unpenalized.nllh + penalty)
}






#' Regularized Binary Regression Model (RBRM)
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
rbrm <- function(va, vb, x, y,
                 alpha.start = NULL, beta.start = NULL,
                 max.step = 3000, thres = 1e-04, lambda = 0,
                 lr.alpha = 0.06, lr.beta = 0.02,
                 intercept = TRUE, early_stopping_rounds = 10) {
  # 
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
        res_alpha <- proximal.gd.alpha.fista(y_alpha, step_size_alpha, 
                                             lambda, t_alpha, last_alpha,
                                             intercept,
                                             beta, va, vb, x, y)
        
        
        alpha_new <- res_alpha$alpha_new
        t_alpha <- res_alpha$t_alpha
        y_alpha <- res_alpha$y_alpha

        # FISTA update for beta
        last_beta <- beta
        res_beta <- proximal.gd.beta.fista(y_beta, step_size_beta, 
                                           lambda, t_beta, last_beta,
                                           intercept,
                                           alpha, va, vb, x, y)
        beta_new <- res_beta$beta_new
        t_beta <- res_beta$t_beta
        y_beta <- res_beta$y_beta

        grad_alpha <- numDeriv::grad(function(.x) {nllh(.x, beta, va, vb, x, y)}, 
                                     y_alpha, method = "simple")
        grad_beta <- numDeriv::grad(function(.x) {nllh(alpha, .x, va, vb, x, y)}, 
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
                                           alpha, beta, va, vb, x, y),
      step = step,
      time = round(time$toc - time$tic, 4)
    )

    return(structure(opt, class = c("rbrm")))
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
                 min_ss = 1e-13, tol_ch = 1e-3) {
  
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
    res_alpha <- proximal.gd.asfista(y_alpha, beta, last_alpha,
                                     opt = "alpha", step_size_alpha, 
                                     lambda, t_alpha,
                                     intercept,
                                     va, vb, x, y)
    
    step_size_alpha <- res_alpha$step_size
    alpha_new <- res_alpha$value_new
    t_alpha <- res_alpha$t_value
    y_alpha <- res_alpha$y_value
    
    # FISTA update for beta
    last_beta <- beta
    res_beta <- proximal.gd.asfista(alpha, y_beta, last_beta,
                                  opt = "beta", step_size_beta, 
                                  lambda, t_beta, 
                                  intercept,
                                  va, vb, x, y)
    
    step_size_beta <- res_beta$step_size
    beta_new <- res_beta$value_new
    t_beta <- res_beta$t_value
    y_beta <- res_beta$y_value
    
    # Update alpha and beta for the next iteration
    alpha <- alpha_new
    beta <- beta_new
    
    grad_alpha <- numDeriv::grad(function(.x) {nllh(.x, beta, va, vb, x, y)}, 
                                 y_alpha, method = "simple")
    
    grad_beta <- numDeriv::grad(function(.x) {nllh(alpha, .x, va, vb, x, y)}, 
                                y_beta, method = "simple")
    
    
    # Stopping criteria based on gradient norm, relative change in objective, and relative change in parameters
    grad_norm_alpha <- max(abs(grad_alpha))
    grad_norm_beta <- max(abs(grad_beta))

    # rel_param_change_alpha <- max(abs(alpha_new - last_beta)) / max(1, max(abs(last_beta)))
    # rel_param_change_beta <- max(abs(beta_new - last_beta)) / max(1, max(abs(last_beta)))
    rel_param_change_alpha <-max(abs(alpha_new - last_alpha) / 
                                   max(1e-5, abs(last_alpha)))
    rel_param_change_beta <- max(abs(beta_new - last_beta) / 
                                   max(1e-5, abs(last_beta)))
    
    
    # Break if all conditions are met
    
    if ((grad_norm_alpha < thres && grad_norm_beta < thres) ||
        (rel_param_change_alpha < tol_ch && rel_param_change_beta < tol_ch)||
        (step_size_alpha < min_ss && step_size_beta < min_ss) 
        # rel_obj_change < tol && 
        ) {
      break
    }
    
  }
  
  time <- tictoc::toc(quiet = TRUE)
  
  opt <- list(
    point.est = c(alpha, beta), 
    convergence = (step < max.step),
    value = penalized.neg.log.likelihood(c(alpha, beta), 
                                         pa, pb, intercept, lambda,
                                         alpha, beta, va, vb, x, y),
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
                              prob_fun = getProbRR.alt) {
  
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
    res_alpha <- proximal.gd.asfista(y_alpha, beta, last_alpha,
                                     opt = "alpha", step_size_alpha, 
                                     lambda, t_alpha,
                                     intercept,
                                     va, vb, x, y,
                                     prob_fun = prob_fun)
    
    step_size_alpha <- res_alpha$step_size
    alpha_new <- res_alpha$value_new
    t_alpha <- res_alpha$t_value
    y_alpha <- res_alpha$y_value
    
    # FISTA update for beta
    last_beta <- beta
    res_beta <- proximal.gd.asfista(alpha, y_beta, last_beta,
                                    opt = "beta", step_size_beta, 
                                    lambda, t_beta, 
                                    intercept,
                                    va, vb, x, y,
                                    prob_fun = prob_fun)
    
    step_size_beta <- res_beta$step_size
    beta_new <- res_beta$value_new
    t_beta <- res_beta$t_value
    y_beta <- res_beta$y_value
    
    # Update alpha and beta for the next iteration
    alpha <- alpha_new
    beta <- beta_new
    
    grad_alpha <- numDeriv::grad(function(.x) {nllh(.x, beta, va, vb, x, y,
                                                    prob_fun = prob_fun)}, 
                                 y_alpha, method = "simple")
    
    grad_beta <- numDeriv::grad(function(.x) {nllh(alpha, .x, va, vb, x, y,
                                                   prob_fun = prob_fun)}, 
                                y_beta, method = "simple")
    
    
    # Stopping criteria based on gradient norm, relative change in objective, and relative change in parameters
    grad_norm_alpha <- max(abs(grad_alpha))
    grad_norm_beta <- max(abs(grad_beta))
    
    # rel_param_change_alpha <- max(abs(alpha_new - last_beta)) / max(1, max(abs(last_beta)))
    # rel_param_change_beta <- max(abs(beta_new - last_beta)) / max(1, max(abs(last_beta)))
    rel_param_change_alpha <-max(abs(alpha_new - last_alpha) / 
                                   max(1e-5, abs(last_alpha)))
    rel_param_change_beta <- max(abs(beta_new - last_beta) / 
                                   max(1e-5, abs(last_beta)))
    
    
    # Break if all conditions are met
    
    if ((grad_norm_alpha < thres && grad_norm_beta < thres) ||
        (rel_param_change_alpha < tol_ch && rel_param_change_beta < tol_ch)||
        (step_size_alpha < min_ss && step_size_beta < min_ss) 
        # rel_obj_change < tol && 
    ) {
      break
    }
    
  }
  
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


#' Regularized Binary Regression Model (RBRM)
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
rbrm.original <- function(va, vb, x, y, # Order corrected
                          alpha.start = NULL, beta.start = NULL,
                 max.step = 3000, thres = 1e-04, lambda = 0,
                 lr.alpha = 0.06, lr.beta = 0.02,
                 intercept = TRUE, early_stopping_rounds = 10) {
  soft_thres <- function(x, lambda) {
    sx <- abs(x) - lambda
    sx[sx < 0] <- 0
    return(sx * sign(x))
  }
  
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
  
  nllh.alpha <- function(alpha) {
    logrr <- (va %*% alpha)
    logop <- (vb %*% beta)
    
    ps <- brm::getProbRR(logrr, logop)
    
    p0 <- ps[, 1]
    p1 <- ps[, 2]
    
    p0[p0 == 1] <- 1 - 1e-15
    p0[p0 <= 0] <- 1e-15
    
    p1[p1 == 1] <- 1 - 1e-15
    p1[p1 <= 0] <- 1e-15
    
    
    if (any(is.nan(log1p(-p1[x == 1])))) {
      p1[which(is.nan(log1p(-p1[x == 1])))] <- 1 - 1e-15
    }
    if (any(is.nan(log1p(-p0[x == 0])))) {
      p1[which(is.nan(log1p(-p0[x == 0])))] <- 1 - 1e-15
    }
    
    return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) +
                  (y[x == 0]) * log(p0[x == 0])) -
             sum((1 - y[x == 1]) * log(1 - p1[x == 1]) +
                   (y[x == 1]) * log(p1[x == 1])))
  }
  
  nllh.beta <- function(beta) {
    logrr <- (va %*% alpha)
    logop <- (vb %*% beta)
    
    ps <- brm::getProbRR(logrr, logop)
    
    p0 <- ps[, 1]
    p1 <- ps[, 2]
    
    p0[p0 == 1] <- 1 - 1e-15
    p0[p0 <= 0] <- 1e-15
    
    p1[p1 == 1] <- 1 - 1e-15
    p1[p1 <= 0] <- 1e-15
    
    if (any(is.nan(log1p(-p1[x == 1])))) {
      p1[which(is.nan(log1p(-p1[x == 1])))] <- 1 - 1e-15
    }
    if (any(is.nan(log1p(-p0[x == 0])))) {
      p1[which(is.nan(log1p(-p0[x == 0])))] <- 1 - 1e-15
    }
    
    return(-sum((1 - y[x == 0]) * log1p(-p0[x == 0]) +
                  (y[x == 0]) * log(p0[x == 0])) -
             sum((1 - y[x == 1]) * log1p(-p1[x == 1]) +
                   (y[x == 1]) * log(p1[x == 1])))
  }
  
  ## FISTA proximal gradient descent functions
  proximal.gd.alpha.fista <- function(alpha, step_size, lambda, t_old, last_alpha) {
    gradient <- numDeriv::grad(nllh.alpha, alpha, method = "simple")
    gradient[is.na(gradient)] <- 0
    input <- alpha - step_size * gradient
    alpha_new <- soft_thres(input, lambda * step_size)
    if (intercept == TRUE) {
      alpha_new[1] <- input[1]
    }
    t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
    y_alpha_new <- alpha_new + (t_old - 1) / t_new * (alpha_new - last_alpha)
    return(list(alpha_new = alpha_new, t_alpha = t_new, y_alpha = y_alpha_new))
  }
  
  proximal.gd.beta.fista <- function(beta, step_size, lambda, t_old, last_beta) {
    gradient <- numDeriv::grad(nllh.beta, beta, method = "simple")
    gradient[is.na(gradient)] <- 0
    input <- beta - step_size * gradient
    beta_new <- soft_thres(input, lambda * step_size)
    if (intercept == TRUE) {
      beta_new[1] <- input[1]
    }
    t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
    y_beta_new <- beta_new + (t_old - 1) / t_new * (beta_new - last_beta)
    return(list(beta_new = beta_new, t_beta = t_new, y_beta = y_beta_new))
  }
  
  penalized.neg.log.likelihood <- function(pars) {
    alpha <- pars[1:pa]
    beta <- pars[(pa + 1):(pa + pb)]
    
    logrr <- (va %*% alpha)
    logop <- (vb %*% beta)
    
    ps <- brm::getProbRR(logrr, logop)
    
    p0 <- ps[, 1]
    p1 <- ps[, 2]
    
    p0[p0 == 1] <- 1 - 1e-15
    p0[p0 <= 0] <- 1e-15
    
    p1[p1 == 1] <- 1 - 1e-15
    p1[p1 <= 0] <- 1e-15
    
    if (any(is.nan(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) +
                        (y[x == 0]) * log(p0[x == 0])) -
                   sum((1 - y[x == 1]) * log(1 - p1[x == 1]) +
                       (y[x == 1]) * log(p1[x == 1]))))) {
      stop("NaN values encountered in unpenalized.nllh. Please check the input matrices.")
    }
    unpenalized.nllh <- (-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) +
                                (y[x == 0]) * log(p0[x == 0])) -
                           sum((1 - y[x == 1]) * log(1 - p1[x == 1]) +
                                 (y[x == 1]) * log(p1[x == 1])))
    if (intercept == TRUE) {
      penalty <- lambda * (sum(abs(alpha[-1])) + sum(abs(beta[-1]))) # penalty does not apply to the intercept coefficient
    } else {
      penalty <- lambda * (sum(abs(alpha)) + sum(abs(beta)))
    }
    return(unpenalized.nllh + penalty)
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
    res_alpha <- proximal.gd.alpha.fista(y_alpha, step_size_alpha,
                                         lambda, t_alpha, last_alpha)
    alpha_new <- res_alpha$alpha_new
    t_alpha <- res_alpha$t_alpha
    y_alpha <- res_alpha$y_alpha
    
    # FISTA update for beta
    last_beta <- beta
    res_beta <- proximal.gd.beta.fista(y_beta, step_size_beta,
                                       lambda, t_beta, last_beta)
    beta_new <- res_beta$beta_new
    t_beta <- res_beta$t_beta
    y_beta <- res_beta$y_beta
    
    grad_alpha <- numDeriv::grad(nllh.alpha, y_alpha, method = "simple")
    grad_beta <- numDeriv::grad(nllh.beta, y_beta, method = "simple")
    if ((max(abs(grad_alpha)) < thres) & (max(abs(grad_beta)) < thres)) {
      break
    }
    
    # Update alpha and beta for the next iteration
    alpha <- alpha_new
    beta <- beta_new
  }
  opt <- list(
    point.est = c(alpha, beta), convergence = (step < max.step),
    value = penalized.neg.log.likelihood(c(alpha, beta)),
    step = step
  )
  
  return(opt)
}

