
#' @export
stop_crit <- function(eval_grad = T,
                      grad_thres = 1e-2,
                      grad_alpha = NULL,
                      grad_beta = NULL,
                      eval_rel_chang = T,
                      eval_rel_grad_thres = 1e-2,
                      alpha = NULL,
                      beta = NULL,
                      last_alpha = NULL,
                      last_beta = NULL,
                      message = F
                      ){
  
  # cli::cli_alert_success("grad_alpha: {grad_alpha} || grad_beta: {grad_beta}")
  
  if(eval_grad){
    
    if (any(is.nan(grad_alpha)) ||
        any(is.nan(grad_beta))) cli::cli_abort("Gradients has NaN(s).")
    
    if (is.null(grad_alpha) ||
        is.null(grad_beta)) cli::cli_abort("No enough information to compute relative change. grad_alpha: {grad_alpha}, grad_beta: {grad_beta}")
    
    grad_norm_alpha <- norm(grad_alpha, type="2")
    grad_norm_beta <- norm(grad_beta, type="2")
    
    grad_eval_return <- (grad_norm_alpha < grad_thres && 
                           grad_norm_beta < grad_thres)
    
    if(message) cli::cli_alert_success("grad_norm_alpha: {grad_norm_alpha} || grad_norm_beta: {grad_norm_beta} || Gradient eval: {grad_eval_return}")
    
  } else {
    grad_eval_return <- F
  }
  
  if(eval_rel_chang){
    
    if (is.null(alpha) ||
        is.null(beta) ||
        is.null(last_alpha) ||
        is.null(last_beta)) cli::cli_abort("No enough information to compute relative change.")
    
    rel_change_alpha <- norm(alpha - last_alpha, type="2") / pmax(1e-8, norm(last_alpha, type="2"))
 
    rel_change_beta <- norm(beta - last_beta, type="2") / pmax(1e-8, norm(last_beta, type="2"))
    
    rel_change_return <- (rel_change_alpha < eval_rel_grad_thres &&
                            rel_change_beta < eval_rel_grad_thres)
    
    if(message) cli::cli_alert_success("rel_change_alpha: {rel_change_alpha} || rel_change_beta: {rel_change_beta} || Relative change: {rel_change_return}")
    
  } else {
    rel_change_return <- F
  }
  return(grad_eval_return || rel_change_return)
  
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
step_fista <- function(alpha, beta,
                       value_old,
                       opt,
                       step_size, lambda, t_old,
                       intercept, va, vb, x, y,
                       prob_fun = getProbRR.org) {
  
  if (!(opt %in% c("alpha","beta"))) {
    cli::cli_abort("Option 'opt' must be either 'alpha' or 'beta'.")
  }
  
  if (opt == "alpha") value <- alpha
  if (opt == "beta") value <- beta
 
  # Momentum update
  t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
  
  a_new <- (t_old-1)/t_new
  
  y_value_new <- value + a_new * (value - value_old)
  
  if (opt == "alpha") {
    
    gradient <- numDeriv::grad(function(.x){nllh(.x, beta, va, vb, x, y,
                                                 prob_fun = prob_fun)}, 
                               y_value_new, method = "simple")
  }
  
  if (opt == "beta") {
    
    gradient <- numDeriv::grad(function(.x) {nllh(alpha, .x, va, vb, x, y,
                                                  prob_fun = prob_fun)},
                               y_value_new, method = "simple")
  }
  
  # Clean any NA gradients to prevent issues during computation
  gradient[is.na(gradient)] <- 0
  
  # Proximal gradient update with soft-thresholding
  input <- y_value_new - step_size * gradient
  
  value_new <- soft_thres(input, lambda * step_size)
  
  # Maintain intercept term if specified
  if (intercept) value_new[1] <- input[1]
  
  
  # Return updated values in a structured list
  return(list(value_new = value_new, t_value = t_new, y_value = y_value_new))
}

#' @export
fista_opt <- function(alpha.start, beta.start,
                      step_size_alpha, step_size_beta,
                      lambda, 
                      intercept, 
                      max.step, 
                      va, vb, x, y,
                      prob_fun = getProbRR.org,
                      opt_step = step_fista,
                      eval_grad = T,
                      eval_rel_chang = T){
  ## Optimization
  step <- 0
  alpha <- y_alpha <- last_alpha <- alpha.start
  beta <- y_beta <- last_beta <- beta.start
  t_alpha <- t_beta <- 1
  
  alphas <- matrix(0, max.step, ncol(v))
  betas <- matrix(0, max.step, ncol(v))
  g_alphas <- matrix(0, max.step, ncol(v))
  g_betas <- matrix(0, max.step, ncol(v))
  nllh_results <- vector("double", max.step)
  
  for (iter in 1:max.step) {
    step <- step + 1
    # FISTA update for alpha
    
    res_alpha <- opt_step(alpha = alpha, beta = beta, 
                          value_old = last_alpha,
                          opt = "alpha",
                          step_size = step_size_alpha, 
                          lambda = lambda, 
                          t_old = t_alpha,
                          intercept = intercept,
                          va = va, vb = vb, x = x, y = y,
                          prob_fun = prob_fun)
    last_alpha <- alpha
    step_size_alpha <- ifelse(is.null(res_alpha$step_size), 
                              step_size_alpha, res_alpha$step_size)
    alpha <- res_alpha$value_new
    t_alpha <- res_alpha$t_value
    y_alpha <- res_alpha$y_value
    
    
    # FISTA update for beta
    res_beta <- opt_step(alpha = alpha, beta = beta, 
                         value_old = last_beta,
                         opt = "beta",
                         step_size = step_size_beta, 
                         lambda = lambda, t_old = t_beta,
                         intercept = intercept,
                         va = va, vb = vb, x = x, y = y,
                         prob_fun = prob_fun)
    
    last_beta <- beta
    
    
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
    
    nllh_iter <- penalized_nllh(alpha, beta, va, vb, x, y, 
                                 lambda = lambda, intercept = intercept,
                                 prob_fun = prob_fun)

    alphas[step,] <- alpha
    betas[step,] <- beta
    g_alphas[step,] <- grad_alpha
    g_betas[step,] <- grad_beta
    nllh_results[step] <- nllh_iter
    
    stop_boolean <- stop_crit(eval_grad = eval_grad,
                              grad_alpha = grad_alpha,
                              grad_beta = grad_beta,
                              eval_rel_chang = eval_rel_chang,
                              alpha = alpha,
                              beta = beta,
                              last_alpha = last_alpha,
                              last_beta = last_beta)
    
    if (stop_boolean) {
      
      alphas <- alphas[1:step,] 
      betas <- betas[1:step,] 
      g_alphas <- g_alphas[1:step,] 
      g_betas <- g_betas[1:step,] 
      nllh_results <- nllh_results[1:step]
      
      break
    }
    
  }
  
  return(list(alpha = alpha,
              beta = beta,
              step = step,
              alphas = alphas,
              betas = betas,
              grad_alphas = g_alphas,
              grad_betas = g_betas,
              nllh_results = nllh_results))
}

fista_opt2 <- function(alpha.start, beta.start,
                      step_size_alpha, step_size_beta,
                      lambda, 
                      intercept, 
                      max.step, 
                      va, vb, x, y,
                      prob_fun = getProbRR.org,
                      opt_step = step_fista){
  ## Optimization
  step <- 0
  alpha <- y_alpha <- last_alpha <- alpha.start
  beta <- y_beta <- last_beta <- beta.start
  t_alpha <- t_beta <- 1
  
  alphas <- matrix(0, max.step, ncol(v))
  betas <- matrix(0, max.step, ncol(v))
  g_alphas <- matrix(0, max.step, ncol(v))
  g_betas <- matrix(0, max.step, ncol(v))
  nllh_results <- vector("double", max.step)
  
  for (iter in 1:max.step) {
    step <- step + 1
    # FISTA update for alpha
    
    res_alpha <- opt_step(alpha = alpha, beta = beta, 
                          value_old = last_alpha,
                          opt = "alpha",
                          step_size = step_size_alpha, 
                          lambda = lambda, 
                          t_old = t_alpha,
                          intercept = intercept,
                          va = va, vb = vb, x = x, y = y,
                          prob_fun = prob_fun)
    
    
    # FISTA update for beta
    res_beta <- opt_step(alpha = alpha, beta = beta, 
                         value_old = last_beta,
                         opt = "beta",
                         step_size = step_size_beta, 
                         lambda = lambda, t_old = t_beta,
                         intercept = intercept,
                         va = va, vb = vb, x = x, y = y,
                         prob_fun = prob_fun)
    
    # Update values 
    last_alpha <- alpha
    step_size_alpha <- ifelse(is.null(res_alpha$step_size), 
                              step_size_alpha, res_alpha$step_size)
    alpha <- res_alpha$value_new
    t_alpha <- res_alpha$t_value
    y_alpha <- res_alpha$y_value
    
    last_beta <- beta
    
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
    
    nllh_iter <- penalized_nllh(alpha, beta, va, vb, x, y, 
                                lambda = lambda, intercept = intercept,
                                prob_fun = prob_fun)
    
    alphas[step,] <- alpha
    betas[step,] <- beta
    g_alphas[step,] <- grad_alpha
    g_betas[step,] <- grad_beta
    nllh_results[step] <- nllh_iter
    
    stop_boolean <- stop_crit(grad_alpha = grad_alpha,
                              grad_beta = grad_beta,
                              alpha = alpha,
                              beta = beta,
                              last_alpha = last_alpha,
                              last_beta = last_beta)
    
    if (stop_boolean) {
      
      alphas <- alphas[1:step,] 
      betas <- betas[1:step,] 
      g_alphas <- g_alphas[1:step,] 
      g_betas <- g_betas[1:step,] 
      nllh_results <- nllh_results[1:step]
      
      break
    }
    
  }
  
  return(list(alpha = alpha,
              beta = beta,
              step = step,
              alphas = alphas,
              betas = betas,
              grad_alphas = g_alphas,
              grad_betas = g_betas,
              nllh_results = nllh_results))
}


#' @export
double_fista_opt <- function(alpha.start, beta.start,
                      step_size_alpha, step_size_beta,
                      lambda, 
                      intercept, 
                      max.step, 
                      va, vb, x, y,
                      prob_fun = getProbRR.org,
                      opt_step = step_fista,
                      cont_opt = 10){
  ## Optimization
  step <- 0
  alpha <- y_alpha <- last_alpha <- alpha.start
  beta <- y_beta <- last_beta <- beta.start
  t_alpha <- t_beta <- 1
  
  max.step <- ceiling(max.step/cont_opt)
  
  alphas <- matrix(0, max.step, ncol(va))
  betas <- matrix(0, max.step, ncol(va))
  g_alphas <- matrix(0, max.step, ncol(va))
  g_betas <- matrix(0, max.step, ncol(va))
  nllh_results <- vector("double", max.step)
  
  for (iter in 1:max.step) {
    step <- step + 1
    
    # FISTA update for beta
    for (i in 1:cont_opt) {
      res_beta <- opt_step(alpha = alpha, beta = beta, 
                           value_old = last_beta,
                           opt = "beta",
                           step_size = step_size_beta, 
                           lambda = lambda, t_old = t_beta,
                           intercept = intercept,
                           va = va, vb = vb, x = x, y = y,
                           prob_fun = prob_fun)
      
      last_beta <- beta
      
      step_size_beta <- ifelse(is.null(res_beta$step_size), 
                               step_size_beta, res_beta$step_size)
      beta <- res_beta$value_new
      t_beta <- res_beta$t_value
      y_beta <- res_beta$y_value
      
    }
    
    
    # FISTA update for alpha
    
    for (i in 1:cont_opt) {
    res_alpha <- opt_step(alpha = alpha, beta = beta, 
                          value_old = last_alpha,
                          opt = "alpha",
                          step_size = step_size_alpha, 
                          lambda = lambda, 
                          t_old = t_alpha,
                          intercept = intercept,
                          va = va, vb = vb, x = x, y = y,
                          prob_fun = prob_fun)
    last_alpha <- alpha
    step_size_alpha <- ifelse(is.null(res_alpha$step_size), 
                              step_size_alpha, res_alpha$step_size)
    alpha <- res_alpha$value_new
    t_alpha <- res_alpha$t_value
    y_alpha <- res_alpha$y_value
    }
    
    
    
    
    # Update alpha and beta for the next iteration
    
    grad_alpha <- numDeriv::grad(function(.x) {nllh(.x, beta, va, vb, x, y,
                                                    prob_fun = prob_fun)}, 
                                 alpha, method = "simple")
    
    grad_beta <- numDeriv::grad(function(.x) {nllh(alpha, .x, va, vb, x, y,
                                                   prob_fun = prob_fun)}, 
                                beta, method = "simple")
    
    
    nllh_iter <- penalized_nllh(alpha, beta, va, vb, x, y, 
                                lambda = lambda, intercept = intercept,
                                prob_fun = prob_fun)
    
    
    # Stopping criteria based on gradient norm, relative change in objective, and relative change in parameters
    # grad_norm_alpha <- max(abs(grad_alpha))
    # grad_norm_beta <- max(abs(grad_beta))
    
    # Clean any NA gradients to prevent issues during computation
    grad_alpha[is.na(grad_alpha)] <- 0
    grad_beta[is.na(grad_beta)] <- 0
    
    
    
    # pmax()
    # Break if all conditions are met
    alphas[step,] <- alpha
    betas[step,] <- beta
    g_alphas[step,] <- grad_alpha
    g_betas[step,] <- grad_beta
    nllh_results[step] <- nllh_iter
    
    
    stop_boolean <- stop_crit(grad_alpha = grad_alpha,
                              grad_beta = grad_beta,
                              alpha = alpha,
                              beta = beta,
                              last_alpha = last_alpha,
                              last_beta = last_beta)
    
    if (stop_boolean) {
      
      alphas <- alphas[1:step,] 
      betas <- betas[1:step,] 
      g_alphas <- g_alphas[1:step,] 
      g_betas <- g_betas[1:step,] 
      nllh_results <- nllh_results[1:step]
      
      break
    }
    
  }
  
  return(list(alpha = alpha,
              beta = beta,
              step = step*cont_opt,
              alphas = alphas,
              betas = betas,
              grad_alphas = g_alphas,
              grad_betas = g_betas,
              nllh_results = nllh_results))
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
                                backtrack_factor = 0.6,  # Step size reduction factor
                                beta_increase = 1.2,  # Factor to increase step size
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
fista <- function(alpha.start, beta.start,
                  step_size_alpha, step_size_beta,
                  lambda, 
                  intercept, 
                  max.step, 
                  va, vb, x, y,
                  prob_fun = getProbRR.org,
                  opt_step = proximal.gd.fista){
  ## Optimization
  step <- 0
  alpha <- y_alpha <- alpha.start
  beta <- y_beta <- beta.start
  t_alpha <- t_beta <- 1
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
    
    
    stop_boolean <- stop_crit(grad_alpha = grad_alpha,
                              grad_beta = grad_beta,
                              alpha = alpha,
                              beta = beta,
                              last_alpha = last_alpha,
                              last_beta = last_beta)
    
    if (stop_boolean) {
      break
    }
    
  }
  
  return(list(alpha = alpha,
              beta = beta,
              step = step))
}

#' ASFISTA
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
step_asfista <- function(alpha, beta,
                         value_old,
                         opt,
                         step_size, lambda, t_old,
                         intercept, va, vb, x, y,
                         prob_fun = getProbRR.org,
                         max_backtrack = 20,  # Max backtracking iterations
                         backtrack_factor = .95,  # Step size reduction factor
                         beta_increase = 1.5  # Factor to increase step size
) {
  
  if (!(opt %in% c("alpha", "beta"))) {
    cli::cli_abort("Option 'opt' must be either 'alpha' or 'beta'.")
  }
  
  if (opt == "alpha") value <- alpha
  if (opt == "beta") value <- beta
  
  # Momentum update
  t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
  a_new <- (t_old - 1) / t_new
  y_value_new <- value + a_new * (value - value_old)
  
  # Compute gradient based on the chosen parameter (alpha or beta)
  if (opt == "alpha") {
    gradient <- numDeriv::grad(function(.x) {
      nllh(.x, beta, va, vb, x, y, prob_fun = prob_fun)
    }, y_value_new, method = "simple")
  } else if (opt == "beta") {
    gradient <- numDeriv::grad(function(.x) {
      nllh(alpha, .x, va, vb, x, y, prob_fun = prob_fun)
    }, y_value_new, method = "simple")
  }
  
  # Clean any NA gradients to prevent computation issues
  gradient[is.na(gradient)] <- 0
  
  # Proximal gradient update with soft-thresholding
  step_size <- step_size * beta_increase
  input <- y_value_new - step_size * gradient
  value_new <- soft_thres(input, lambda * step_size)
  
  # Maintain intercept term if specified
  if (intercept) value_new[1] <- input[1]
  
  # Backtracking to ensure sufficient decrease
  for (bt_iter in 1:max_backtrack) {
    if (opt == "alpha") {
      loss_new <- nllh(value_new, beta, va, vb, x, y, prob_fun = prob_fun)
      loss_old <- nllh(value, beta, va, vb, x, y, prob_fun = prob_fun)
    } else if (opt == "beta") {
      loss_new <- nllh(alpha, value_new, va, vb, x, y, prob_fun = prob_fun)
      loss_old <- nllh(alpha, value, va, vb, x, y, prob_fun = prob_fun)
    }
    
    # Armijo-like condition
    sufficient_adj <- loss_new <= loss_old - 
      (sum((value_new - y_value_new)^2) / (2 * step_size))
    
    if (is.na(sufficient_adj)) {
    cli::cli_alert_success("loss_new {loss_new} || loss_old {loss_old} || step_size: {step_size}")
      sufficient_adj <- T
      }
    
    
    if (sufficient_adj) {

      break
      
      } else {
      
      step_size <- step_size * backtrack_factor
      if (step_size > 0.5) step_size <- 0.5
      
      input <- y_value_new - step_size * gradient
      
      value_new <- soft_thres(input, lambda * step_size)
    
    }
    
    
  }
  # cli::cli_alert_info("{step_size}")
  
  # Return updated values in a structured list
  return(list(
    value_new = value_new,
    t_value = t_new,
    y_value = y_value_new,
    step_size = step_size
  ))
}


#' @export
asfista <- function(alpha.start, beta.start,
                    step_size_alpha, step_size_beta,
                    lambda, 
                    intercept, 
                    max.step, 
                    va, vb, x, y,
                    prob_fun = getProbRR.org,
                    opt_step = step_asfista){
  ## Optimization
  step <- 0
  alpha <- y_alpha <- last_alpha <- alpha.start
  beta <- y_beta <- last_beta <- beta.start
  t_alpha <- t_beta <- 1
  
  alphas <- matrix(0, max.step, ncol(v))
  betas <- matrix(0, max.step, ncol(v))
  g_alphas <- matrix(0, max.step, ncol(v))
  g_betas <- matrix(0, max.step, ncol(v))
  nllh_results <- vector("double", max.step)
  
  step_size_alpha_loop <- step_size_alpha 
  step_size_beta_loop <- step_size_beta 
  
  
  for (iter in 1:max.step) {
    step <- step + 1
    # FISTA update for alpha
    
    # cli::cli_alert_success("step_size_alpha {step_size_alpha} || step_size_beta {step_size_beta}")
    
    res_alpha <- opt_step(alpha = alpha, beta = beta, 
                          value_old = last_alpha,
                          opt = "alpha",
                          step_size = step_size_alpha_loop, 
                          lambda = lambda, 
                          t_old = t_alpha,
                          intercept = intercept,
                          va = va, vb = vb, x = x, y = y,
                          prob_fun = prob_fun)
    last_alpha <- alpha
    step_size_alpha_loop <- ifelse(is.null(res_alpha$step_size),
                                   step_size_alpha_loop, res_alpha$step_size)
    alpha <- res_alpha$value_new
    t_alpha <- res_alpha$t_value
    y_alpha <- res_alpha$y_value
    
    
    # FISTA update for beta
    res_beta <- opt_step(alpha = alpha, beta = beta, 
                         value_old = last_beta,
                         opt = "beta",
                         step_size = step_size_beta_loop, 
                         lambda = lambda, t_old = t_beta,
                         intercept = intercept,
                         va = va, vb = vb, x = x, y = y,
                         prob_fun = prob_fun)
    
    last_beta <- beta
    
    
    step_size_beta_loop <- ifelse(is.null(res_beta$step_size),
                                  step_size_beta_loop, res_beta$step_size)
    
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
    
    nllh_iter <- penalized_nllh(alpha, beta, va, vb, x, y, 
                                lambda = lambda, intercept = intercept,
                                prob_fun = prob_fun)
    
    alphas[step,] <- alpha
    betas[step,] <- beta
    g_alphas[step,] <- grad_alpha
    g_betas[step,] <- grad_beta
    nllh_results[step] <- nllh_iter
    
    stop_boolean <- stop_crit(grad_alpha = grad_alpha,
                              grad_beta = grad_beta,
                              alpha = alpha,
                              beta = beta,
                              last_alpha = last_alpha,
                              last_beta = last_beta)
    
    # if (norm(grad_alpha , type="2") > 1e-2 &
    #     (step_size_alpha_loop < (step_size_alpha))) step_size_alpha_loop <- step_size_alpha #/ ceiling(step/100)
    # if (norm(grad_beta , type="2") > 1e-2 &
    #     (step_size_beta_loop < (step_size_beta))) step_size_beta_loop <- step_size_beta #/ ceiling(step/100)
  
      
  
    
    if (stop_boolean) {
      
      alphas <- alphas[1:step,] 
      betas <- betas[1:step,] 
      g_alphas <- g_alphas[1:step,] 
      g_betas <- g_betas[1:step,] 
      nllh_results <- nllh_results[1:step]
      
      break
    }
    
  }
  
  return(list(alpha = alpha,
              beta = beta,
              step = step,
              alphas = alphas,
              betas = betas,
              grad_alphas = g_alphas,
              grad_betas = g_betas,
              nllh_results = nllh_results))
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
greedy_fista <- function(alpha.start, beta.start,
    step_size_alpha, step_size_beta,
    lambda, 
    intercept, 
    max.step, thres,
    va, vb, x, y,
    prob_fun = getProbRR.org,
    opt_step = proximal.gd.fista,
    safeguard_factor = 10, gamma_shrink = 0.9
) {
  # Initialize variables
  step <- 0
  initial_diff <- NULL
  alpha <- y_alpha <- alpha.start
  beta <- y_beta <- beta.start
  t_alpha <- t_beta <- 1
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
                 max.step = 1000, lambda = 0,
                 lr.alpha = 0.01, lr.beta = 0.01,
                 intercept = TRUE,
                 prob_fun = getProbRR.org,
                 opt_fun = fista , save_opt = T) {
  
  # va <- v; vb <- v; alpha.start = NULL; beta.start = NULL;
  # max.step = 1000;  lambda = 0;
  # lr.alpha = 0.06; lr.beta = 0.02;
  # intercept = TRUE;  prob_fun = getProbRR.org;
  # opt_fun = fista
  tictoc::tic("Total time")
  
  if (is.null(vb)) {
    vb <- va
  }
  
  va <- as.matrix(va)
  vb <- as.matrix(vb)
  
  pa <- dim(va)[2]
  pb <- dim(vb)[2]
  
  # sanity check for the intercept term in va, vb
  
  if (all(va[, 1] == 1) & all(vb[, 1] == 1)) {
    intercept <- TRUE
  }
  ## starting values for parameter optimization
  if (is.null(alpha.start)) alpha.start <- c(rep(0, pa))
  if (length(alpha.start) < pa) alpha.start <- c(rep(alpha.start[1], pa))
  
  if (is.null(beta.start)) beta.start <- c(rep(0.01, pb))
  if (length(beta.start) < pa) beta.start <- c(rep(beta.start[1], pa))
  
  
  ## Optimization

  opt_result <- opt_fun(alpha.start, beta.start,
                        lr.alpha, lr.beta,
                        lambda, 
                        intercept, 
                        max.step,
                        va, vb, x, y,
                        prob_fun)
  
  step <- opt_result$step
  alpha <- opt_result$alpha
  beta <- opt_result$beta
  
  if(!save_opt) opt_result <- NULL
  
  time <- tictoc::toc(quiet = TRUE)
  
  opt <- list(
    point.est = c(alpha, beta), 
    optimization.info = opt_result,
    convergence = (step < max.step),
    value = penalized_nllh(alpha, beta, 
                           va, vb, x, y,
                           lambda, intercept),
    step = step,
    time = round(time$toc - time$tic, 4)
  )
  
  return(structure(opt, class = c("rbrm")))
}



