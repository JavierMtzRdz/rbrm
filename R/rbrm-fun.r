#' @export
stop_crit <- function(eval_grad = T,
                      grad_thres = 1e-7,
                      grad_alpha = NULL,
                      grad_beta = NULL,
                      eval_rel_chang = T,
                      eval_rel_grad_thres = 1e-5,
                      alpha = NULL,
                      beta = NULL,
                      last_alpha = NULL,
                      last_beta = NULL,
                      stop_nan = T,
                      message = F,
                      message_true = F
                      ){
  
  if (any(is.nan(grad_alpha)) ||
      any(is.nan(grad_beta))) {
    cli::cli_alert_danger("Gradients has NaN(s).")
    
    return(T)
  }
  
  if(eval_grad){
    
    if (is.null(grad_alpha) ||
        is.null(grad_beta)) cli::cli_abort("No enough information to compute relative change. grad_alpha: {grad_alpha}, grad_beta: {grad_beta}")
    
    grad_norm_alpha <- norm(grad_alpha, type="2")
    grad_norm_beta <- norm(grad_beta, type="2")
    
    grad_eval_return <- (grad_norm_alpha < grad_thres && 
                           grad_norm_beta < grad_thres)
    
    if(message | (message_true & grad_eval_return)) cli::cli_alert_success("grad_norm_alpha: {grad_norm_alpha} || grad_norm_beta: {grad_norm_beta} || Gradient eval: {grad_eval_return}")
    
  } else {
    grad_eval_return <- F
  }
  
  if(eval_rel_chang){
    
    if (is.null(alpha) ||
        is.null(beta) ||
        is.null(last_alpha) ||
        is.null(last_beta)) cli::cli_abort("No enough information to compute relative change.")
    
    # rel_change_alpha <- norm(alpha - last_alpha, type="2")^2 / pmax(1e-08, norm(alpha, type="2"))^2
    # 
    # rel_change_beta <- norm(beta - last_beta, type="2")^2 / pmax(1e-08, norm(beta, type="2"))^2
    
    # rel_change_return <- (rel_change_alpha < eval_rel_grad_thres &&
    #                         rel_change_beta < eval_rel_grad_thres)
    
    Diff = function(x,y) sum((x-y)^2)/sum(x^2+thres)
    
    rel_change_alpha <- sum((alpha-last_alpha)^2)/sum(alpha^2+eval_rel_grad_thres)
    
    rel_change_beta <- sum((beta-last_beta)^2)/sum(beta^2+eval_rel_grad_thres)
    
    diff <-  max(rel_change_alpha, rel_change_beta)
    
    rel_change_return <- (diff < eval_rel_grad_thres)
    
    if(message | (message_true & rel_change_return)) cli::cli_alert_success("rel_change_alpha: {round(rel_change_alpha, 5)} || rel_change_beta: {round(rel_change_beta, 5)}")
    
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
                       opt = c("alpha", "beta"),
                       step_size, lambda, t_old,
                       intercept, va, vb, x, y,
                       prob_fun = getProbRR.org) {
  
  opt <- rlang::arg_match(opt)
  
  if (opt == "alpha") value <- alpha
  if (opt == "beta") value <- beta
  
  # Momentum update
  t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
  
  a_new <- (t_old-1)/t_new
  
  y_value_new <- value + a_new * (value - value_old)
  
  if (opt == "alpha") gradient <- grad_nll(y_value_new, beta,
                                           x, y, va, vb,
                                           prob_fun, opt = "alpha")
  
  
  if (opt == "beta") gradient <- grad_nll(alpha, y_value_new, 
                                          x, y, va, vb,
                                          prob_fun, opt = "beta")
  
  # Clean any NA gradients to prevent issues during computation
  if (any(is.na(gradient))) {
    cli::cli_alert_danger("NaN in gradient, replacing with 0.")
    print(gradient)
    gradient[is.na(gradient)] <- 0
  }
  
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
                      max_step, 
                      va, vb, x, y,
                      prob_fun = getProbRR.org,
                      opt_step = step_fista,
                      eval_grad = T){
  ## Optimization
  step <- 0
  alpha <- y_alpha <- last_alpha <- alpha.start
  beta <- y_beta <- last_beta <- beta.start
  t_alpha <- t_beta <- 1
  
  alphas <- matrix(0, max_step, ncol(va))
  betas <- matrix(0, max_step, ncol(vb))
  g_alphas <- matrix(0, max_step, ncol(va))
  g_betas <- matrix(0, max_step, ncol(vb))
  nllh_results <- vector("double", max_step)
  
  for (iter in 1:max_step) {
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
    
    if(eval_grad){
      grad <- grad_nll(alpha, beta, x, y, va, vb,
                       prob_fun)
      grad_alpha <- grad$grad_alpha
      grad_beta <- grad$grad_beta
    } else {
      grad_alpha <- grad_beta <- 99
    }
    
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
                              # eval_rel_chang = eval_rel_chang,
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
L <- function(alpha, beta, va, vb, x, y, prob_fun, 
              opt) {
  
  if (opt == "alpha") {
    hessian_matrix <- numDeriv::hessian(function(.x) {
      nllh(.x, beta, va, vb, x, y, prob_fun = prob_fun)
    }, alpha)
  }
  
  if (opt == "beta") {
    hessian_matrix <- numDeriv::hessian(function(.x) {
      nllh(alpha, .x, va, vb, x, y, prob_fun = prob_fun)
    }, beta)
  }
  
  
  # Compute the largest eigenvalue (Lipschitz constant)
  L <- max(eigen(hessian_matrix, symmetric = TRUE, only.values = TRUE)$values)
  
  return(L)
}

#' @export
fista_opt2 <- function(alpha.start, beta.start,
                       step_size_alpha, step_size_beta,
                       lambda, 
                       intercept, 
                       max_step, 
                       va, vb, x, y,
                       prob_fun = getProbRR.org,
                       opt_step = step_fista,
                       eval_grad = T){
  ## Optimization
  step <- 0
  
  last_alpha <- alpha.start + 1
  last_beta <- beta.start + 1
  
  alpha <- y_alpha <- alpha.start
  beta <- y_beta <- beta.start
  t_alpha <- t_beta <- 1
  
  alphas <- matrix(0, max_step, ncol(va))
  betas <- matrix(0, max_step, ncol(vb))
  g_alphas <- matrix(0, max_step, ncol(va))
  g_betas <- matrix(0, max_step, ncol(vb))
  nllh_results <- vector("double", max_step)
  
  # L_alpha <- L(alpha, beta, va, vb, x, y, prob_fun,
  #              opt = "alpha")
  # 
  # L_beta <- L(alpha, beta, va, vb, x, y, prob_fun,
  #             opt = "beta")
  # 
  # step_size_alpha <- 1/L_alpha
  # step_size_beta <- 1/L_beta

     # cli::cli_alert("step_size_alpha: {step_size_alpha} | step_size_beta: {step_size_beta}")
  
  for (iter in 1:max_step) {
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
    
    
    # Update alpha and beta for the next iteration
    
    # Update values 
    last_alpha <- alpha
    step_size_alpha <- ifelse(is.null(res_alpha$step_size), 
                              step_size_alpha, res_alpha$step_size)
    alpha <- res_alpha$value_new
    t_alpha <- res_alpha$t_value
    last_y_alpha <- y_alpha
    y_alpha <- res_alpha$y_value
    
    last_beta <- beta
    
    step_size_beta <- ifelse(is.null(res_beta$step_size), 
                             step_size_beta, res_beta$step_size)
    beta <- res_beta$value_new
    t_beta <- res_beta$t_value
    last_y_beta <- y_beta
    y_beta <- res_beta$y_value
    
    if(eval_grad){
      grad <- grad_nll(alpha, beta, x, y, va, vb,
                             prob_fun)
      grad_alpha <- grad$grad_alpha
      grad_beta <- grad$grad_beta
      
    } else {
      grad_alpha <- grad_beta <- 99
    }
    
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


#' Coordinate Descent Optimization for Alpha and Beta
#'
#' Applies coordinate descent to optimize alpha and beta for the
#' L1-penalized negative log-likelihood.
#'
#' @param alpha.start Initial numeric vector for alpha coefficients.
#' @param beta.start Initial numeric vector for beta coefficients.
#' @param step_size_alpha Numeric value for the step size used in alpha updates.
#' @param step_size_beta Numeric value for the step size used in beta updates.
#' @param lambda Numeric value for the L1 regularization parameter.
#' @param intercept Logical. If TRUE, the first coefficient of alpha and beta is not penalized.
#' @param max_step Maximum number of iterations.
#' @param va Matrix of independent variables for alpha.
#' @param vb Matrix of independent variables for beta.
#' @param x Vector indicating group assignment (e.g., exposure).
#' @param y Vector of binary outcomes.
#' @param prob_fun Function to calculate probabilities (e.g., getProbRR.org).
#' @param nllh_fun Function to compute the negative log-likelihood (without penalty).
#'                   Its signature should be: nllh_fun(alpha, beta, va, vb, x, y, prob_fun).
#' @param penalized_nllh_fun Function to compute the penalized negative log-likelihood.
#'                   Its signature should be: penalized_nllh_fun(alpha, beta, va, vb, x, y, lambda, intercept, prob_fun).
#' @param tol Tolerance for convergence.
#' @param verbose Logical, if TRUE, prints progress.
#' @return A list with the optimized alpha, beta, number of iterations, and history.



#' @export
cd_opt <- function(alpha.start, beta.start,
                   step_size_alpha, step_size_beta,
                   lambda,
                   intercept,
                   max_step,
                   va, vb, x, y,
                   prob_fun = getProbRR.org,
                   nllh_fun = nllh,  
                   penalized_nllh_fun = penalized_nllh, 
                   tol = 1e-6,
                   verbose = FALSE)  {
  
  alpha <- alpha.start
  beta <- beta.start
  p_alpha <- length(alpha)
  p_beta <- length(beta)
  
  # History storage
  alphas_hist <- matrix(NA, nrow = max_step, ncol = p_alpha)
  betas_hist <- matrix(NA, nrow = max_step, ncol = p_beta)
  grad_alphas_hist <- matrix(NA, nrow = max_step, ncol = p_alpha)
  grad_betas_hist <- matrix(NA, nrow = max_step, ncol = p_beta)
  pen_nllh_values_hist <- numeric(max_step) # Stores penalized NLLH
  
  iter_count <- 0
  
  for (current_iter in 1:max_step) {
    iter_count <- current_iter
    alpha_old_iter <- alpha
    beta_old_iter <- beta
    
    # --- Update alpha coefficients ---
    for (k in 1:p_alpha) {
      
      grad_ak <- grad_nll_k(alpha, beta,
                          x, y, va, vb,
                          prob_fun, opt = "alpha",
                          k_index = k)
      grad_alphas_hist[iter_count, k] <- grad_ak
      alpha_k_unreg_update <- alpha[k] - step_size_alpha * grad_ak
      alpha[k] <- soft_thres(alpha_k_unreg_update, lambda * step_size_alpha)
      
      grad_bk <- grad_nll_k(alpha, beta, 
                          x, y, va, vb,
                          prob_fun, opt = "beta",
                          k_index = k)
      grad_betas_hist[iter_count, k] <- grad_bk
      beta_k_unreg_update <- beta[k] - step_size_beta * grad_bk
      beta[k] <- soft_thres(beta_k_unreg_update, lambda * step_size_beta)
    }
    
    # --- Store history for this iteration ---
    if(p_alpha > 0) alphas_hist[current_iter, ] <- alpha
    if(p_beta > 0) betas_hist[current_iter, ] <- beta
    
    
    current_pen_nllh <- penalized_nllh_fun(alpha, beta, va, vb, x, y, lambda, intercept, prob_fun)
    pen_nllh_values_hist[current_iter] <- current_pen_nllh
    
    if (verbose && (current_iter %% 10 == 0)) {
      cat("Iter: ", current_iter, ", Penalized NLLH: ", current_pen_nllh, "\n")
    }
    
    # --- Check convergence ---
    delta_alpha <- max(abs(alpha - alpha_old_iter))
    delta_beta <- max(abs(beta - beta_old_iter))
    
    # Ensure there's at least one parameter to check for change
    if (delta_alpha < tol && delta_beta < tol) {
      if (verbose) cat("Converged: Max coefficient change below tolerance at iteration ", current_iter, "\n")
      break
    }
    if (current_iter > 2) {
      rel_change_obj <- abs(pen_nllh_values_hist[current_iter] - pen_nllh_values_hist[current_iter - 1]) /
        (abs(pen_nllh_values_hist[current_iter - 1]) + 1e-8)
      if (rel_change_obj < tol) {
        if (verbose) cat("Converged: Relative change in objective below tolerance at iter ", current_iter, "\n")
        break
      }
    }
  }
  
  if (iter_count == max_step && verbose && max_step > 0) {
    cat("Reached max iterations (", max_step, ").\n")
  }
  
  # Truncate history to actual number of iterations
  actual_iters <- iter_count
  alphas_hist_out <- if(p_alpha > 0) alphas_hist[1:actual_iters, , drop = FALSE] else matrix(NA, nrow=actual_iters, ncol=0)
  betas_hist_out <- if(p_beta > 0) betas_hist[1:actual_iters, , drop = FALSE] else matrix(NA, nrow=actual_iters, ncol=0)
  grad_alphas_hist_out <- grad_alphas_hist[1:actual_iters, , drop = FALSE]
  grad_betas_hist_out <- grad_betas_hist[1:actual_iters, , drop = FALSE] 
  nllh_results_out <- pen_nllh_values_hist[1:actual_iters]
  
  
  
  return(list(alpha = alpha,
              beta = beta,
              step = iter_count, 
              alphas = alphas_hist_out,
              betas = betas_hist_out,
              grad_alphas = grad_alphas_hist_out,
              grad_betas = grad_betas_hist_out,
              nllh_results = nllh_results_out)) 
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
#' @param max_step Maximum number of optimization steps. Default is 3000.
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
rbrm.exp <- function(va, vb, x, y,
                              alpha.start = NULL, beta.start = NULL,
                              max_step = 1000, lambda = 0,
                              lr.alpha = 1, lr.beta = 1,
                              intercept = FALSE,
                              prob_fun = getProbRR.org,
                              opt_fun = fista , save_opt = T) {
  
  # va <- v; vb <- v; alpha.start = NULL; beta.start = NULL;
  # max_step = 1000;  lambda = 0;
  # lr.alpha = 0.06; lr.beta = 0.02;
  # intercept = TRUE;  prob_fun = getProbRR.org;
  # opt_fun = fista
  
  tictoc::tic("Total time")
  
  if (is.null(vb)) {
    vb <- va
  }
  
  va <- as.matrix(va)
  vb <- as.matrix(vb)
  
  # Add intercept column if not already present
  # if (intercept) {
  #   va <- cbind(1, va)
  #   vb <- cbind(1, vb)
  # }
  
  
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
                        max_step,
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
    convergence = (step < max_step),
    value = penalized_nllh(alpha, beta, 
                           va, vb, x, y,
                           lambda, intercept),
    step = step,
    time = round(time$toc - time$tic, 4)
  )
  
  return(structure(opt, class = c("rbrm")))
}


#' @export
step_bp <- function(alpha, beta,
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
    # browser()
    gradient <- grad_nll(y_value_new, beta,
                         x, y, va, vb,
                         prob_fun, opt = "alpha")
    
    lambda_step <- ((sum(abs(alpha))+ 1e-10)/(sum(abs(beta), 
                                                 abs(alpha))+2e-10))*lambda*2
  }
  
  if (opt == "beta") {
    # browser()
    gradient <- grad_nll(alpha, y_value_new, 
                         x, y, va, vb,
                         prob_fun, opt = "beta")
    
    lambda_step <- ((sum(abs(beta))+ 1e-10)/(sum(abs(beta), 
                                                   abs(alpha))+2e-10))*lambda*2
  }
  
  # Clean any NA gradients to prevent issues during computation
  if (any(is.na(gradient))) {
    cli::cli_alert_danger("NaN in gradient, replacing with 0.")
    print(gradient)
    gradient[is.na(gradient)] <- 0
  }
  
  # Proximal gradient update with soft-thresholding
  input <- y_value_new - step_size * gradient
  # value_new <- soft_thres(input, lambda_step * step_size)
  value_new <-  try(soft_thres(input, lambda_step * step_size))
  
  if(any(class(value_new) == "try-error")) browser()
  # if(any(abs(value_new) > 3)) browser()
  
  # Maintain intercept term if specified
  if (intercept) value_new[1] <- input[1]
  
  # Return updated values in a structured list
  return(list(value_new = value_new, t_value = t_new, y_value = y_value_new))
}


#' @export
double_fista_opt <- function(alpha.start, beta.start,
                      step_size_alpha, step_size_beta,
                      lambda, 
                      intercept, 
                      max_step, 
                      va, vb, x, y,
                      prob_fun = getProbRR.org,
                      opt_step = step_fista,
                      cont_opt = 10,
                      eval_grad = T){
  ## Optimization
  step <- 0
  alpha <- y_alpha <- last_alpha <- alpha.start
  beta <- y_beta <- last_beta <- beta.start
  t_alpha <- t_beta <- 1
  
  max_step <- ceiling(max_step/cont_opt)
  
  alphas <- matrix(0, max_step, ncol(va))
  betas <- matrix(0, max_step, ncol(va))
  g_alphas <- matrix(0, max_step, ncol(va))
  g_betas <- matrix(0, max_step, ncol(va))
  nllh_results <- vector("double", max_step)
  
  for (iter in 1:max_step) {
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
    if(eval_grad){
      grad <- grad_nll(alpha, beta, x, y, va, vb,
                       prob_fun)
      grad_alpha <- grad$grad_alpha
      grad_beta <- grad$grad_beta
      
    } else {
      grad_alpha <- grad_beta <- 99
    }
    
    
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
                  max_step, 
                  va, vb, x, y,
                  prob_fun = getProbRR.org,
                  opt_step = proximal.gd.fista,
                  grad_alpha = T){
  ## Optimization
  step <- 0
  alpha <- y_alpha <- alpha.start
  beta <- y_beta <- beta.start
  t_alpha <- t_beta <- 1
  for (iter in 1:max_step) {
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
    
    if(eval_grad){
      grad_alpha <- numDeriv::grad(function(.x) {nllh(.x, beta, va, vb, x, y,
                                                      prob_fun = prob_fun)}, 
                                   alpha, method = "simple")
      
      grad_beta <- numDeriv::grad(function(.x) {nllh(alpha, .x, va, vb, x, y,
                                                     prob_fun = prob_fun)}, 
                                  beta, method = "simple")
    } else {
      grad_alpha <- grad_beta <- 99
    }
    
    
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
  if (opt == "alpha") gradient <- grad_nll(y_value_new, beta,
                                           x, y, va, vb,
                                           prob_fun, opt = "alpha")
  
  
  if (opt == "beta") gradient <- grad_nll(alpha, y_value_new, 
                                          x, y, va, vb,
                                          prob_fun, opt = "beta")
  
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
                    max_step, 
                    va, vb, x, y,
                    prob_fun = getProbRR.org,
                    opt_step = step_asfista,
                    eval_grad = T){
  ## Optimization
  step <- 0
  alpha <- y_alpha <- last_alpha <- alpha.start
  beta <- y_beta <- last_beta <- beta.start
  t_alpha <- t_beta <- 1
  
  alphas <- matrix(0, max_step, ncol(v))
  betas <- matrix(0, max_step, ncol(v))
  g_alphas <- matrix(0, max_step, ncol(v))
  g_betas <- matrix(0, max_step, ncol(v))
  nllh_results <- vector("double", max_step)
  
  step_size_alpha_loop <- step_size_alpha 
  step_size_beta_loop <- step_size_beta 
  
  
  for (iter in 1:max_step) {
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
    
    if(eval_grad){
      grad_alpha <- numDeriv::grad(function(.x) {nllh(.x, beta, va, vb, x, y,
                                                      prob_fun = prob_fun)}, 
                                   alpha, method = "simple")
      
      grad_beta <- numDeriv::grad(function(.x) {nllh(alpha, .x, va, vb, x, y,
                                                     prob_fun = prob_fun)}, 
                                  beta, method = "simple")
    } else {
      grad_alpha <- grad_beta <- 99
    }
    
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
    max_step, thres,
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
  for (iter in 1:max_step) {
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


#' Initialize Coefficient Starting Vectors
#'
#' Handles NULL or incorrectly sized starting vectors, returning a vector
#' of the expected length.
#'
#' @param start_vec User-provided starting vector (or NULL).
#' @param expected_len The required length (number of columns).
#' @param default_val The default value to use if start_vec is NULL or invalid.
#' @param vec_name Character name of the vector for warning messages.
#' @return A numeric vector of length expected_len.
#' @keywords internal
.initialize_start_coeffs <- function(start_vec, expected_len, default_val = 0, vec_name = "coeffs") {
  # Handle zero-column case
  if (expected_len <= 0) {
    return(numeric(0))
  }
  
  if (is.null(start_vec)) {
    # Default initialization
    final_vec <- rep(default_val, expected_len)
  } else if (length(start_vec) == expected_len) {
    # Use user-provided if length is correct
    final_vec <- start_vec
  } else {
    # Provided vector has wrong length
    cli::cli_warn(
      "Length of {vec_name} ({length(start_vec)}) != expected ({expected_len}). Using default value {default_val} instead."
    )
    final_vec <- rep(default_val, expected_len)
  }
  return(as.numeric(final_vec)) # Ensure numeric type
}

#' Experimental RBRM Model Fitting Function
#'
#' Wraps an optimization function to fit the RBRM model with Lasso penalty.
#'
#' @param va Matrix of predictors for alpha.
#' @param vb Matrix of predictors for beta (defaults to `va` if `NULL`).
#' @param x Binary treatment vector (0/1).
#' @param y Binary outcome vector (0/1).
#' @param alpha_start Optional starting vector for alpha coefficients.
#' @param beta_start Optional starting vector for beta coefficients.
#' @param max_step Maximum iterations for the optimizer.
#' @param lambda Lasso penalty strength (non-negative scalar).
#' @param step_size_alpha_init Initial step size guess for alpha (used by some optimizers).
#' @param step_size_beta_init Initial step size guess for beta (used by some optimizers).
#' @param intercept Logical. Does the model include an intercept term (as the first
#'   element of alpha/beta)? Controls penalization. User should ensure `va`/`vb`
#'   matrices include/exclude a column of 1s accordingly.
#' @param prob_fun Function to calculate probabilities (e.g., `getProbRR.org`).
#' @param opt_fun The optimization function to use (e.g., `fista_opt2_ls_sc`).
#'   Must accept specific arguments (see code).
#' @param save_optimizer_details Logical. If `TRUE`, include the full raw output
#'   from `opt_fun` in the results.
#' @param ... Additional arguments passed directly to `opt_fun` (e.g., `tol`,
#'   `ls_max_iter`, `eval_grad`).
#'
#' @return An object of class "rbrm" (or similar), typically a list containing
#'   estimated coefficients, convergence status, objective value, etc.
#'
#' @importFrom utils modifyList
#' @importFrom cli cli_abort cli_warn cli_alert_info
#' @importFrom tictoc tic toc
#' @export
rbrm.exp2 <- function(va, vb = NULL, x, y,
                      alpha.start = NULL, beta.start = NULL,
                      max_step = 1000, lambda = 0,
                      lr.alpha = 1, lr.beta = 1,
                      intercept = F,
                      prob_fun = getProbRR.org,    
                      opt_fun = fista_opt, 
                      save_opt = F) {
  tictoc::tic("rbrm_experimental time")

  if (lambda < 0) { cli::cli_warn("lambda is negative ({lambda}), using 0 instead."); lambda <- 0 }

  if (is.null(vb)) {
    cli::cli_alert_info("vb is NULL, using va for beta predictors.")
    vb <- va
  }
  # Ensure matrix format
  va <- tryCatch(as.matrix(va), error = function(e) cli::cli_abort("Failed to coerce 'va' to matrix: {e$message}"))
  vb <- tryCatch(as.matrix(vb), error = function(e) cli::cli_abort("Failed to coerce 'vb' to matrix: {e$message}"))

  n <- length(y)
  pa <- ncol(va)
  pb <- ncol(vb)

  # Check intercept column based on user flag (guidance only)
  has_intercept_col <- isTRUE(all(va[, 1] == 1)) 
  if (intercept && !has_intercept_col) {
    cli::cli_warn("intercept=TRUE but a column of 1s was not detected as the first column of 'va'. Ensure data includes intercept if needed.")
  }
  if (!intercept && has_intercept_col) {
    cli::cli_warn("intercept=FALSE but a column of 1s was detected as the first column of 'va'. Ensure data excludes intercept if not desired.")
  }

  # --- 2. Initialize Starting Values ---
  alpha_start <- .initialize_start_coeffs(alpha.start, pa, default_val = 0, "alpha_start")
  beta_start  <- .initialize_start_coeffs(beta.start, pb, default_val = 0.01, "beta_start") # Default beta to 0


  # Prepare arguments list
  opt_args <- list(
    alpha.start = alpha_start,
    beta.start = beta_start,
    step_size_alpha = lr.alpha,
    step_size_beta = lr.beta,
    lambda = lambda,
    intercept = intercept, 
    max_step = max_step,
    va = va, vb = vb, x = x, y = y,
    prob_fun = prob_fun
  )
  
  # Call the optimizer
  opt_result <- tryCatch({
    do.call(opt_fun, opt_args)
  }, error = function(e){
    cli::cli_abort("Optimization failed: {e$message}") # Abort if optimizer itself errors
  })

  # Extract Results ---
  step  <- opt_result$step
  alpha <- opt_result$alpha
  beta  <- opt_result$beta
  
  # Check if optimizer returned expected results
  if(is.null(step) || is.null(alpha) || is.null(beta) || length(alpha) != pa || length(beta) != pb){
    cli::cli_abort("Optimizer returned NULL or coefficients/step of incorrect length! Check 'opt_fun'. alpha: {length(alpha)} (exp {pa}), beta: {length(beta)} (exp {pb})")
  }
  
  # Objective Value
  final_value <- tryCatch({penalized_nllh(alpha, beta, va, vb, x, y, lambda, intercept, prob_fun = prob_fun)}, error = function(e){
    cli::cli_warn("Calculation of final penalized NLLH failed: {e$message}")
    NA_real_
  })

  
  # Structure Output ---
  time_info <- tictoc::toc(quiet = TRUE)
  run_time <- round(time_info$toc - time_info$tic, 4)
  
  if(!save_opt)  opt_result <- NULL
  
  result <- list(
    point.est = c(alpha, beta),
    alpha = alpha,
    beta = beta,
    convergence = (!is.null(step) && is.numeric(step) && step < max_step), 
    step = step,
    optimizer_details = opt_result,
    lambda = lambda,
    intercept = intercept,
    dimensions = list(n = n, p_a = pa, p_b = pb),
    time = run_time
  )
  # cli::cli_alert_success("rbrm_experimental finished in {run_time} seconds.")
  return(structure(result, class = c("rbrm")))
}



opt_mle <- function(alpha.start, beta.start,
                     step_size_alpha, step_size_beta,
                     lambda, 
                     intercept, 
                     max_step, 
                     va, vb, x, y,
                     prob_fun = getProbRR.org,
                    thres = 1e-08) {
  
  
  alphas <- matrix(0, max_step, ncol(v))
  betas <- matrix(0, max_step, ncol(v))
  g_alphas <- matrix(0, max_step, ncol(v))
  g_betas <- matrix(0, max_step, ncol(v))
  nllh_results <- vector("double", max_step)
  
  Diff = function(x,y) sum((x-y)^2)/sum(x^2+thres)
  alpha = alpha.start; beta = beta.start
  diff = thres + 1; step = 0
  while(diff > thres & step < max_step){
    step = step + 1
    opt1 = stats::optim(alpha,
                        function(.x){penalized_nllh(.x, beta, va, vb, x, y, 
                                       lambda = lambda, intercept = intercept,
                                       prob_fun = prob_fun)},
                        control=list(maxit=max(100,max_step/10)))
    diff1 = Diff(opt1$par,alpha)
    alpha = opt1$par
    opt2 = stats::optim(beta,
                        function(.x){penalized_nllh(alpha, .x, va, vb, x, y, 
                                                    lambda = lambda, intercept = intercept,
                                                    prob_fun = prob_fun)},
                        ,control=list(maxit=max(100,max_step/10)))
    diff  = max(diff1,Diff(opt2$par,beta))
    beta = opt2$par
    nllh_iter <- penalized_nllh(alpha, beta, va, vb, x, y, 
                                lambda = lambda, intercept = intercept,
                                prob_fun = prob_fun)
    
    alphas[step,] <- alpha
    betas[step,] <- beta
    g_alphas[step,] <- 0
    g_betas[step,] <- 0
    nllh_results[step] <- nllh_iter
  }
  
  alphas <- alphas[1:step,] 
  betas <- betas[1:step,] 
  g_alphas <- g_alphas[1:step,] 
  g_betas <- g_betas[1:step,] 
  nllh_results <- nllh_results[1:step]
  
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
step_fista_ls <- function(current_param_val, # Current value of param being optimized (alpha_k or beta_k)
                          other_param_val,   # The other param (beta_k or updated alpha_{k+1})
                          param_val_old,     # Param value from iter k-1 (alpha_{k-1} or beta_{k-1})
                          opt_target,        # "alpha" or "beta"
                          initial_step_size, # Starting step size for line search
                          lambda,
                          t_old,
                          intercept,
                          va, vb, x, y,
                          prob_fun,
                          shrink_factor = 0.5,
                          max_ls_iter = 20) {
  
  # 1. FISTA Momentum update
  t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
  momentum_coeff <- (t_old - 1) / t_new
  # y_k = x_k + w_k * (x_k - x_{k-1})
  y_extrapolated <- current_param_val + momentum_coeff * (current_param_val - param_val_old)
  
  # 2. Gradient of the smooth part (NLL) at the extrapolated point y_extrapolated
  grad_f_at_y <- NULL
  if (opt_target == "alpha") {
    grad_f_at_y <- grad_nll(y_extrapolated, other_param_val, y, x, va, vb, prob_fun, opt = "alpha")
  } else { # opt_target == "beta"
    grad_f_at_y <- grad_nll(other_param_val, y_extrapolated, y, x, va, vb, prob_fun, opt = "beta")
  }
  
  if (any(is.na(grad_f_at_y))) { grad_f_at_y[is.na(grad_f_at_y)] <- 0 } # Simplified NA handling
  
  # 3. NLL value at y_extrapolated ( f(y_k) )
  f_val_at_y <- NULL
  if (opt_target == "alpha") {
    f_val_at_y <- nllh(y_extrapolated, other_param_val, va, vb, x, y, prob_fun)
  } else {
    f_val_at_y <- nllh(other_param_val, y_extrapolated, va, vb, x, y, prob_fun)
  }
  
  # 4. Line search loop
  current_s <- initial_step_size
  param_new_accepted <- NULL
  
  for (ls_iter in 1:max_ls_iter) {
    # Candidate update: x_{k+1} = prox_{s*g}(y_k - s * grad_f(y_k))
    prox_arg <- y_extrapolated - current_s * grad_f_at_y
    param_candidate <- soft_thres(prox_arg, lambda * current_s)
    if (intercept) { param_candidate[1] <- prox_arg[1] }
    
    # NLL value at param_candidate ( f(x_{k+1}) )
    f_val_at_candidate <- NULL
    if (opt_target == "alpha") {
      f_val_at_candidate <- nllh(param_candidate, other_param_val, va, vb, x, y, prob_fun)
    } else {
      f_val_at_candidate <- nllh(other_param_val, param_candidate, va, vb, x, y, prob_fun)
    }
    
    # Backtracking condition: f(x_new) <= f(y) + <grad_f(y), x_new-y> + (1/(2s))||x_new-y||^2
    # This is the standard condition for proximal gradient methods.
    rhs_condition <- f_val_at_y +
      sum(grad_f_at_y * (param_candidate - y_extrapolated)) +
      (1 / (2 * current_s)) * sum((param_candidate - y_extrapolated)^2)
    
    if (f_val_at_candidate <= rhs_condition + 1e-9) { # Added small tolerance
      param_new_accepted <- param_candidate
      break
    }
    current_s <- current_s * shrink_factor
  }
  
  if (is.null(param_new_accepted)) { # Line search failed to satisfy condition
    # Default to using the smallest step tried, or could issue a warning/error
    prox_arg <- y_extrapolated - current_s * grad_f_at_y 
    param_new_accepted <- soft_thres(prox_arg, lambda * current_s)
    if (intercept) { param_new_accepted[1] <- prox_arg[1] }
    # warning(paste("Line search for", opt_target, "may not have converged; using step_size =", current_s))
  }
  
  return(list(value_new = param_new_accepted, t_value = t_new, 
              y_value = y_extrapolated, step_size = current_s))
}

#' @export
fista_opt2_ls <- function(alpha.start, beta.start,
                          step_size_alpha, step_size_beta,
                                      lambda, intercept, max_step,
                                      va, vb, x, y,
                                      prob_fun,           # For nllh_fun, penalized_nllh_fun, grad_nll_fun
                                      eval_grad = TRUE,
                                      ls_shrink_factor = 0.5,
                                      ls_max_iter = 20) { 
  
  alpha <- alpha.start
  beta <- beta.start
  
  last_alpha <- alpha.start # alpha_{k-1} for the first iteration (k=0)
  last_beta <- beta.start   # beta_{k-1} for the first iteration (k=0)
  
  t_alpha <- 1.0
  t_beta <- 1.0
  
  current_s_alpha <- step_size_alpha
  current_s_beta <- step_size_beta
  
  # History storage
  p_alpha <- length(alpha.start)
  p_beta <- length(beta.start)
  alphas_hist <- matrix(NA_real_, nrow = max_step, ncol = p_alpha)
  betas_hist <- matrix(NA_real_, nrow = max_step, ncol = p_beta)
  grad_alphas_hist_sc <- matrix(NA_real_, nrow = max_step, ncol = p_alpha) # For stopping criteria
  grad_betas_hist_sc <- matrix(NA_real_, nrow = max_step, ncol = p_beta)   # For stopping criteria
  nllh_results_hist <- vector("double", max_step)
  
  final_iter <- 0
  
  for (iter in 1:max_step) {
    final_iter <- iter
    
    alpha_k_start <- alpha # Value of alpha at the beginning of iteration k
    beta_k_start <- beta   # Value of beta at the beginning of iteration k
    
    # FISTA update for alpha
    res_alpha <- step_fista_ls(
      current_param_val = alpha_k_start,    # x_k
      other_param_val = beta_k_start,     # beta_k (used to calculate grad NLL for alpha)
      param_val_old = last_alpha,         # x_{k-1}
      opt_target = "alpha",
      initial_step_size = current_s_alpha,
      lambda = lambda, t_old = t_alpha, intercept = intercept,
      va = va, vb = vb, x = x, y = y, prob_fun = prob_fun,
      shrink_factor = ls_shrink_factor, max_ls_iter = ls_max_iter
    )
    alpha_next <- res_alpha$value_new # This is x_{k+1} for alpha
    
    # FISTA update for beta
    res_beta <- step_fista_ls(
      current_param_val = beta_k_start,     # x_k for beta
      other_param_val = alpha_next,       # Use updated alpha for beta's gradient calc
      param_val_old = last_beta,          # x_{k-1} for beta
      opt_target = "beta",
      initial_step_size = current_s_beta,
      lambda = lambda, t_old = t_beta, intercept = intercept,
      va = va, vb = vb, x = x, y = y, prob_fun = prob_fun,
      shrink_factor = ls_shrink_factor, max_ls_iter = ls_max_iter
    )
    beta_next <- res_beta$value_new # This is x_{k+1} for beta
    
    # Update values for the next iteration (k becomes k-1, k+1 becomes k)
    last_alpha <- alpha_k_start
    last_beta <- beta_k_start
    
    alpha <- alpha_next
    beta <- beta_next
    
    t_alpha <- res_alpha$t_value
    t_beta <- res_beta$t_value
    
    current_s_alpha <- res_alpha$step_size
    current_s_beta <- res_beta$step_size
    
    # Store results
    if(p_alpha > 0) alphas_hist[iter, ] <- alpha
    if(p_beta > 0) betas_hist[iter, ] <- beta
    nllh_results_hist[iter] <- penalized_nllh(alpha, beta, va, vb, x, y,
                                                  lambda = lambda, intercept = intercept,
                                                  prob_fun = prob_fun)
    
    # Gradient for stopping criterion
    
    if (eval_grad) {
      # If not precomputed, calculate them now based on *updated* alpha and beta
      grads_sc <- grad_nll(alpha, beta, x, y, va, vb, prob_fun)
      g_alpha_sc <- grads_sc$grad_alpha
      g_beta_sc <- grads_sc$grad_beta
    }
    
    if(eval_grad){
      grad_alphas_hist_sc[iter, ] <- g_alpha_sc
      grad_betas_hist_sc[iter, ] <- g_beta_sc
    }
    
    # Check stopping criterion
    # Pass current alpha/beta and alpha/beta from start of this iteration
    # (alpha_k_start, beta_k_start act as "last_alpha", "last_beta" for change calculation)
    if (iter > 0) { # iter > 1 if stop_crit needs change from previous iter
      stop_boolean <- stop_crit(
        grad_alpha = if(eval_grad && p_alpha > 0) g_alpha_sc else NULL,
        grad_beta = if(eval_grad && p_beta > 0) g_beta_sc else NULL,
        alpha = alpha, beta = beta,
        last_alpha = alpha_k_start, last_beta = beta_k_start,
        eval_grad = eval_grad
        # ... ensure all necessary args for stop_crit_fun are passed
      )
      if (stop_boolean) {
        break
      }
    }
  }
  
  # Truncate history matrices
  alphas_hist <- alphas_hist[1:final_iter, , drop = FALSE]
  betas_hist <- betas_hist[1:final_iter, , drop = FALSE]
  nllh_results_hist <- nllh_results_hist[1:final_iter]
  grad_alphas_hist_sc <- grad_alphas_hist_sc[1:final_iter, , drop = FALSE]
  grad_betas_hist_sc <- grad_betas_hist_sc[1:final_iter, , drop = FALSE]
  
  return(list(
    alpha = alpha, beta = beta,
    step = final_iter,
    alphas = alphas_hist, betas = betas_hist,
    grad_alphas = grad_alphas_hist_sc, grad_betas = grad_betas_hist_sc,
    nllh_results = nllh_results_hist,
    final_step_size_alpha = current_s_alpha,
    final_step_size_beta = current_s_beta
  ))
}
