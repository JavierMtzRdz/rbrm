#' @export
stop_crit <- function(eval_grad = T,
                      grad_thres = 1e-4,
                      grad_alpha = NULL,
                      grad_beta = NULL,
                      eval_rel_chang = T,
                      eval_rel_grad_thres = 1e-4,
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
    
    rel_change_alpha <- norm(alpha - last_alpha, type="2")^2 / pmax(1e-08, norm(alpha, type="2"))^2
 
    rel_change_beta <- norm(beta - last_beta, type="2")^2 / pmax(1e-08, norm(beta, type="2"))^2
    
    rel_change_return <- (rel_change_alpha < eval_rel_grad_thres &&
                            rel_change_beta < eval_rel_grad_thres)
    
    if(message | (message_true & rel_change_return)) cli::cli_alert_success("rel_change_alpha: {round(rel_change_alpha, 5)} || rel_change_beta: {round(rel_change_beta, 5)}")
    
  } else {
    rel_change_return <- F
  }
  return(grad_eval_return || rel_change_return)
  
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
    # browser()
    gradient <- numDeriv::grad(function(.x){nllh(.x, beta, va, vb, x, y,
                                                 prob_fun = prob_fun)},
                               y_value_new, method = "simple")
    # gradient <- -grad_nll_alpha(y_value_new, beta, x, y, va, vb, prob_fun)
    
    
  }
  
  if (opt == "beta") {
    # browser()
    gradient <- numDeriv::grad(function(.x) {nllh(alpha, .x, va, vb, x, y,
                                                  prob_fun = prob_fun)},
                               y_value_new, method = "simple")
    # gradient <- -grad_nll_beta(alpha, y_value_new, x, y, va, vb, prob_fun)
  
  }

  
  
  
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
step_fista_an <- function(alpha, beta,
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
    # gradient <- numDeriv::grad(function(.x){nllh(.x, beta, va, vb, x, y,
    #                                              prob_fun = prob_fun)},
    #                            y_value_new, method = "simple")
    
    
    gradient <- grad_nll(alpha, y_value_new, x, y, va, vb,
                         prob_fun, opt = "alpha")

    # if(any(abs(gradient - gradient2) > 3)) {browser()} else {gradient <- gradient2}
    
    
  }
  
  if (opt == "beta") {
    
    # gradient <- numDeriv::grad(function(.x) {nllh(alpha, .x, va, vb, x, y,
    #                                               prob_fun = prob_fun)},
    #                            y_value_new, method = "simple")

    gradient <- grad_nll(alpha, y_value_new, x, y, va, vb,
                                      prob_fun, opt = "beta")

    # if(any(abs(gradient - gradient2) > 1.5)) {browser()} else {gradient <- gradient2}

  }
  
  
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
                      max.step, 
                      va, vb, x, y,
                      prob_fun = getProbRR.org,
                      opt_step = step_fista,
                      eval_grad = T){
  ## Optimization
  step <- 0
  alpha <- y_alpha <- last_alpha <- alpha.start
  beta <- y_beta <- last_beta <- beta.start
  t_alpha <- t_beta <- 1
  
  alphas <- matrix(0, max.step, ncol(va))
  betas <- matrix(0, max.step, ncol(vb))
  g_alphas <- matrix(0, max.step, ncol(va))
  g_betas <- matrix(0, max.step, ncol(vb))
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
fista_opt_an <- function(alpha.start, beta.start,
                      step_size_alpha, step_size_beta,
                      lambda, 
                      intercept, 
                      max.step, 
                      va, vb, x, y,
                      prob_fun = getProbRR.org,
                      opt_step = step_fista,
                      eval_grad = T){
  ## Optimization
  step <- 0
  alpha <- y_alpha <- last_alpha <- alpha.start
  beta <- y_beta <- last_beta <- beta.start
  t_alpha <- t_beta <- 1
  
  alphas <- matrix(0, max.step, ncol(va))
  betas <- matrix(0, max.step, ncol(vb))
  g_alphas <- matrix(0, max.step, ncol(va))
  g_betas <- matrix(0, max.step, ncol(vb))
  nllh_results <- vector("double", max.step)
  

  
  # L_alpha <- L(alpha, beta, va, vb, x, y, prob_fun,
  #              opt = "alpha")
  # 
  # L_beta <- L(alpha, beta, va, vb, x, y, prob_fun,
  #             opt = "beta")
  # 
  # step_size_alpha <- 1/L_alpha
  # step_size_beta <- 1/L_beta
  # 
  #    cli::cli_alert("step_size_alpha: {step_size_alpha} | step_size_beta: {step_size_beta}")
    

  
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
    

    grad_alpha <- grad_nll(alpha, beta, x, y, va, vb,
                         prob_fun, opt = "alpha")
    
    grad_beta <- grad_nll(alpha, beta, x, y, va, vb,
                            prob_fun, opt = "beta")
    
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
                              eval_rel_grad_thres = 1e-4,
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
fista_opt2 <- function(alpha.start, beta.start,
                      step_size_alpha, step_size_beta,
                      lambda, 
                      intercept, 
                      max.step, 
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
  
  alphas <- matrix(0, max.step, ncol(va))
  betas <- matrix(0, max.step, ncol(vb))
  g_alphas <- matrix(0, max.step, ncol(va))
  g_betas <- matrix(0, max.step, ncol(vb))
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
    
    if (stop_boolean) {
      
      alphas <- alphas[1:step,] 
      betas <- betas[1:step,] 
      g_alphas <- g_alphas[1:step,] 
      g_betas <- g_betas[1:step,] 
      nllh_results <- nllh_results[1:step]
      
      break
    } else {
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
                      cont_opt = 10,
                      eval_grad = T){
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
                  opt_step = proximal.gd.fista,
                  grad_alpha = T){
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
                    opt_step = step_asfista,
                    eval_grad = T){
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
                 intercept = FALSE,
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
rbrm.experimental2 <- function(va, vb = NULL, x, y,
                               alpha.start = NULL, beta.start = NULL,
                              max_step = 1000, lambda = 0,
                              lr.alpha = 0.01, lr.beta = 0.01,
                              intercept = F,
                              prob_fun = getProbRR.org,    
                              opt_fun = fista_opt, 
                              save_opt = F) {
  # browser()
  tictoc::tic("rbrm_experimental time")

  # --- 1. Input Validation and Preparation ---
  if (!is.function(opt_fun)) cli::cli_abort("'opt_fun' must be a function.")
  if (!is.null(prob_fun) && !is.function(prob_fun)) cli::cli_abort("If provided, 'prob_fun' must be a function.")
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
  if (nrow(va) != n || nrow(vb) != n || length(x) != n) {
    cli::cli_abort("Input dimension mismatch (va, vb, x, y rows/lengths).")
  }

  # Check intercept column based on user flag (guidance only)
  has_intercept_col <- pa > 0 && isTRUE(all(va[, 1] == 1)) # Check only va
  if (intercept && !has_intercept_col) {
    cli::cli_warn("intercept=TRUE but a column of 1s was not detected as the first column of 'va'. Ensure data includes intercept if needed.")
  }
  if (!intercept && has_intercept_col) {
    cli::cli_warn("intercept=FALSE but a column of 1s was detected as the first column of 'va'. Ensure data excludes intercept if not desired.")
  }

  # --- 2. Initialize Starting Values ---
  alpha_start_final <- .initialize_start_coeffs(alpha.start, pa, default_val = 0, "alpha_start")
  beta_start_final  <- .initialize_start_coeffs(beta.start, pb, default_val = 0.01, "beta_start") # Default beta to 0


  # Prepare arguments list
  opt_args <- list(
    alpha.start = alpha_start_final,
    beta.start = beta_start_final,
    # Pass step sizes using names expected
    step_size_alpha = lr.alpha,
    step_size_beta = lr.beta,
    lambda = lambda,
    intercept = intercept, # Pass the user's intent
    max.step = max_step,
    va = va, vb = vb, x = x, y = y,
    prob_fun = prob_fun
  )
  
  # Call the optimizer
  opt_result <- tryCatch({
    do.call(opt_fun, opt_args)
  }, error = function(e){
    cli::cli_abort("Optimization failed: {e$message}") # Abort if optimizer itself errors
  })

# browser()
  # --- 4. Validate and Extract Results ---
  step  <- opt_result$step
  alpha <- opt_result$alpha
  beta  <- opt_result$beta
  
  # Check if optimizer returned expected results
  if(is.null(step) || is.null(alpha) || is.null(beta) || length(alpha) != pa || length(beta) != pb){
    cli::cli_abort("Optimizer returned NULL or coefficients/step of incorrect length! Check 'opt_fun'. alpha: {length(alpha)} (exp {pa}), beta: {length(beta)} (exp {pb})")
  }
  
  # --- 5. Calculate Final Objective Value ---
  final_value <- tryCatch({penalized_nllh(alpha, beta, va, vb, x, y, lambda, intercept, prob_fun = prob_fun)}, error = function(e){
    cli::cli_warn("Calculation of final penalized NLLH failed: {e$message}")
    NA_real_
  })
  if (!is.finite(final_value)) {
    cli::cli_warn("Final penalized NLLH is non-finite (NA/Inf).")
  }
  
  # --- 6. Structure Output ---
  time_info <- tictoc::toc(quiet = TRUE)
  run_time <- round(time_info$toc - time_info$tic, 4)
  
  if(!save_opt)  opt_result <- NULL
  
  result <- list(
    point.est = c(alpha, beta),
    alpha = alpha,
    beta = beta,
    convergence = (!is.null(step) && is.numeric(step) && step < max_step), # Check step validity
    value = final_value, # Penalized NLLH
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
                     max.step, 
                     va, vb, x, y,
                     prob_fun = getProbRR.org) {
  
  thres = 1e-08
  alphas <- matrix(0, max.step, ncol(v))
  betas <- matrix(0, max.step, ncol(v))
  g_alphas <- matrix(0, max.step, ncol(v))
  g_betas <- matrix(0, max.step, ncol(v))
  nllh_results <- vector("double", max.step)
  
  Diff = function(x,y) sum((x-y)^2)/sum(x^2+thres)
  alpha = alpha.start; beta = beta.start
  diff = thres + 1; step = 0
  while(diff > thres & step < max.step){
    step = step + 1
    opt1 = stats::optim(alpha,
                        function(.x){penalized_nllh(.x, beta, va, vb, x, y, 
                                       lambda = lambda, intercept = intercept,
                                       prob_fun = prob_fun)},
                        control=list(maxit=max(100,max.step/10)))
    diff1 = Diff(opt1$par,alpha)
    alpha = opt1$par
    opt2 = stats::optim(beta,
                        function(.x){penalized_nllh(alpha, .x, va, vb, x, y, 
                                                    lambda = lambda, intercept = intercept,
                                                    prob_fun = prob_fun)},
                        ,control=list(maxit=max(100,max.step/10)))
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


fista_opt2_ls <-  function(alpha.start, beta.start,
                           # Initial step size guesses (will be adapted)
                           step_size_alpha_init = 1.0,
                           step_size_beta_init = 1.0,
                           lambda,
                           intercept, # Boolean: Is the first element an intercept?
                           max.step,
                           va, vb, x, y,
                           prob_fun = getProbRR.org,
                           # Stopping criterion function (MUST handle NULL gradients if eval_grad=FALSE)
                           # Line Search parameters
                           ls_beta_shrink = 0.5, # Step size reduction factor
                           ls_max_iter = 20,     # Max iterations for line search step
                           # Gradient function (replace with analytic if possible!)
                           grad_fun = numDeriv::grad,
                           grad_method = "simple", # Method for numDeriv::grad
                           # Objective function components
                           nllh_fun = nllh,
                           penalized_nllh_fun = penalized_nllh,
                           soft_thres_fun = soft_thres,
                           # Other controls
                           eval_grad = TRUE, # Evaluate gradients for stop_crit?
                           tol = 1e-6        # Tolerance (can be used within stop_crit)
){
  
  ## Input dimensions
  pa <- ncol(va)
  pb <- ncol(vb)
  if(length(alpha.start) != pa) { cli::cli_abort("alpha.start length mismatch: %d vs %d", length(alpha.start), pa); }
  if(length(beta.start) != pb) { cli::cli_abort("beta.start length mismatch: %d vs %d", length(beta.start), pb); }
  
  ## Initialization
  step <- 0
  alpha <- alpha.start
  beta <- beta.start
  last_alpha <- alpha.start # Store previous iteration's value for momentum
  last_beta <- beta.start
  t_alpha <- 1
  t_beta <- 1
  current_step_size_alpha <- step_size_alpha_init
  current_step_size_beta <- step_size_beta_init
  
  # Preallocate storage
  alphas <- matrix(NA_real_, nrow = max.step + 1, ncol = pa)
  betas <- matrix(NA_real_, nrow = max.step + 1, ncol = pb)
  # Store gradients computed at momentum points 'y' during update/line search
  g_alphas <- if(pa > 0) matrix(NA_real_, nrow = max.step + 1, ncol = pa) else matrix(NA_real_, 0, 0)
  g_betas <- if(pb > 0) matrix(NA_real_, nrow = max.step + 1, ncol = pb) else matrix(NA_real_, 0, 0)
  nllh_results <- vector("double", max.step + 1)
  
  # Store initial values
  if(pa > 0) alphas[1,] <- alpha.start
  if(pb > 0) betas[1,] <- beta.start
  nllh_results[1] <- tryCatch(penalized_nllh_fun(alpha, beta, va, vb, x, y, lambda, intercept, prob_fun, nllh_fun), error = function(e) NA_real_)
  
  ## Optimization Loop
  for (iter in 1:max.step) {
    step <- iter
    storage_idx <- iter + 1
    
    alpha_iter_start <- alpha # Value at the start of this iteration
    beta_iter_start <- beta
    
    # --- Alpha Update with Line Search ---
    grad_alpha_y <- NULL # Gradient computed at y_alpha
    if (pa > 0) {
      t_new_alpha <- (1 + sqrt(1 + 4 * t_alpha^2)) / 2
      momentum_coeff_alpha <- (t_alpha - 1) / t_new_alpha
      if (!is.finite(momentum_coeff_alpha)) momentum_coeff_alpha <- 0
      y_alpha <- alpha + momentum_coeff_alpha * (alpha - last_alpha)
      
      grad_alpha_y <- tryCatch({ # grad at y_alpha
        grad_fun(function(.a) nllh_fun(.a, beta, va, vb, x, y, prob_fun),
                 y_alpha, method = grad_method)
      }, error = function(e) { NULL })
      
      if (is.null(grad_alpha_y) || length(grad_alpha_y) != pa || any(!is.finite(grad_alpha_y))) {
        cli::cli_warn("Invalid gradient grad_alpha_y at step %d. Skipping alpha update.", step)
        grad_alpha_y <- rep(NA_real_, pa) # Store NA gradient
        t_new_alpha <- t_alpha # Keep old momentum
      } else {
        f_y_alpha <- nllh_fun(y_alpha, beta, va, vb, x, y, prob_fun)
        if(!is.finite(f_y_alpha)){
          cli::cli_warn("Non-finite NLLH at y_alpha (step %d). Skipping alpha update.", step)
          t_new_alpha <- t_alpha
          grad_alpha_y <- rep(NA_real_, pa) # Store NA gradient
        } else {
          ss_alpha <- current_step_size_alpha
          alpha_updated <- FALSE
          for (ls_iter in 1:ls_max_iter) {
            prox_arg_alpha <- y_alpha - ss_alpha * grad_alpha_y
            trial_alpha <- soft_thres_fun(prox_arg_alpha, lambda * ss_alpha)
            if (intercept && pa > 0) trial_alpha[1] <- prox_arg_alpha[1]
            f_trial_alpha <- nllh_fun(trial_alpha, beta, va, vb, x, y, prob_fun)
            if (!is.finite(f_trial_alpha)) {
              ss_alpha <- ss_alpha * ls_beta_shrink; next
            }
            quad_approx <- f_y_alpha + sum(grad_alpha_y * (trial_alpha - y_alpha)) +
              (0.5 / ss_alpha) * sum((trial_alpha - y_alpha)^2)
            if (f_trial_alpha <= quad_approx + 1e-8) {
              alpha <- trial_alpha; current_step_size_alpha <- ss_alpha
              t_alpha <- t_new_alpha; alpha_updated <- TRUE; break
            } else {
              ss_alpha <- ss_alpha * ls_beta_shrink
            }
          }
          if (!alpha_updated) {
            cli::cli_warn("Line search for alpha failed at iteration %d. Skipping alpha update.", step)
            t_new_alpha <- t_alpha; grad_alpha_y <- rep(NA_real_, pa) # Store NA
          }
        }
      }
      # Store the gradient computed at y_alpha (even if NA)
      g_alphas[storage_idx, ] <- grad_alpha_y
    } else {
      if (pa > 0) g_alphas[storage_idx, ] <- NA_real_ # No update, store NA
    }
    
    
    # --- Beta Update with Line Search (using updated alpha) ---
    grad_beta_y <- NULL # Gradient computed at y_beta
    if (pb > 0) {
      t_new_beta <- (1 + sqrt(1 + 4 * t_beta^2)) / 2
      momentum_coeff_beta <- (t_beta - 1) / t_new_beta
      if (!is.finite(momentum_coeff_beta)) momentum_coeff_beta <- 0
      y_beta <- beta + momentum_coeff_beta * (beta - last_beta)
      
      grad_beta_y <- tryCatch({ # grad at y_beta
        grad_fun(function(.b) nllh_fun(alpha, .b, va, vb, x, y, prob_fun),
                 y_beta, method = grad_method)
      }, error = function(e) { NULL })
      
      if (is.null(grad_beta_y) || length(grad_beta_y) != pb || any(!is.finite(grad_beta_y))) {
        cli::cli_warn("Invalid gradient grad_beta_y at step %d. Skipping beta update.", step)
        grad_beta_y <- rep(NA_real_, pb) # Store NA gradient
        t_new_beta <- t_beta # Keep old momentum
      } else {
        f_y_beta <- nllh_fun(alpha, y_beta, va, vb, x, y, prob_fun)
        if(!is.finite(f_y_beta)){
          cli::cli_warn("Non-finite NLLH at y_beta (step %d). Skipping beta update.", step)
          t_new_beta <- t_beta
          grad_beta_y <- rep(NA_real_, pb) # Store NA gradient
        } else {
          ss_beta <- current_step_size_beta
          beta_updated <- FALSE
          for (ls_iter in 1:ls_max_iter) {
            prox_arg_beta <- y_beta - ss_beta * grad_beta_y
            trial_beta <- soft_thres_fun(prox_arg_beta, lambda * ss_beta)
            if (intercept && pb > 0) trial_beta[1] <- prox_arg_beta[1]
            f_trial_beta <- nllh_fun(alpha, trial_beta, va, vb, x, y, prob_fun)
            if (!is.finite(f_trial_beta)) {
              ss_beta <- ss_beta * ls_beta_shrink; next
            }
            quad_approx <- f_y_beta + sum(grad_beta_y * (trial_beta - y_beta)) +
              (0.5 / ss_beta) * sum((trial_beta - y_beta)^2)
            if (f_trial_beta <= quad_approx + 1e-8) {
              beta <- trial_beta; current_step_size_beta <- ss_beta
              t_beta <- t_new_beta; beta_updated <- TRUE; break
            } else {
              ss_beta <- ss_beta * ls_beta_shrink
            }
          }
          if (!beta_updated) {
            cli::cli_warn("Line search for beta failed at iteration %d. Skipping beta update.", step)
            t_new_beta <- t_beta; grad_beta_y <- rep(NA_real_, pb) # Store NA
          }
        }
      }
      # Store the gradient computed at y_beta (even if NA)
      g_betas[storage_idx, ] <- grad_beta_y
    } else {
      if (pb > 0) g_betas[storage_idx, ] <- NA_real_ # No update, store NA
    }
    
    # --- Store results ---
    if (pa > 0) alphas[storage_idx,] <- alpha
    if (pb > 0) betas[storage_idx,] <- beta
    nllh_iter <- tryCatch({
      penalized_nllh_fun(alpha, beta, va, vb, x, y,
                         lambda = lambda, intercept = intercept, prob_fun = prob_fun, nllh_fun = nllh_fun)
    }, error = function(e) { NA_real_ })
    nllh_results[storage_idx] <- nllh_iter
    
    # --- Check Stopping Criterion ---
    grad_alpha_for_stop <- NULL
    grad_beta_for_stop <- NULL
    
    if (eval_grad) {
      # Calculate gradients at the updated alpha/beta *if* needed by stop_crit
      if (pa > 0) {
        grad_alpha_for_stop <- tryCatch({
          grad_fun(function(.a) nllh_fun(.a, beta, va, vb, x, y, prob_fun),
                   alpha, method = grad_method)
        }, error = function(e) { cli::cli_warn("Stop_crit grad alpha failed: %s", e$message); NULL })
        if (is.null(grad_alpha_for_stop) || length(grad_alpha_for_stop) != pa || any(!is.finite(grad_alpha_for_stop))) {
          grad_alpha_for_stop <- NULL # Ensure it's NULL if invalid
        }
      }
      if (pb > 0) {
        grad_beta_for_stop <- tryCatch({
          grad_fun(function(.b) nllh_fun(alpha, .b, va, vb, x, y, prob_fun),
                   beta, method = grad_method)
        }, error = function(e) { cli::cli_warn("Stop_crit grad beta failed: %s", e$message); NULL })
        if (is.null(grad_beta_for_stop) || length(grad_beta_for_stop) != pb || any(!is.finite(grad_beta_for_stop))) {
          grad_beta_for_stop <- NULL # Ensure it's NULL if invalid
        }
      }
    }
    
    # Call the user-provided stop_crit function
    # It MUST handle NULL for grad_alpha_for_stop/grad_beta_for_stop if eval_grad=FALSE
    # stop_boolean <- stop_crit(
    #   grad_alpha = grad_alpha_for_stop,
    #   grad_beta = grad_beta_for_stop,
    #   alpha = alpha,
    #   beta = beta,
    #   last_alpha = alpha_iter_start, # Value from start of this iteration
    #   last_beta = beta_iter_start,   # Value from start of this iteration
    #   iter = iter,
    #   max_step = max.step,
    #   tol = tol
    #   # Add any other arguments your specific stop_crit function needs
    # )
    stop_boolean <- stop_crit(eval_grad = eval_grad,
                              grad_alpha = grad_alpha_for_stop,
                              grad_beta = grad_beta_for_stop,
                              # eval_rel_chang = eval_rel_chang,
                              alpha = alpha,
                              beta = beta,
                              last_alpha = last_alpha,
                              last_beta = last_beta)
    
    if (stop_boolean) {
      # cli::cli_alert_info("Convergence reached at step %d by stop_crit.", step)
      break # Exit main optimization loop
    }
    
    # --- Update 'last' values for next iteration's momentum ---
    # Use values from the *start* of the current iteration
    last_alpha <- alpha_iter_start
    last_beta <- beta_iter_start
    
    # Note: t_alpha and t_beta were updated within the line search success blocks
    
  } # End main optimization loop (for iter ...)
  
  # --- Prepare results ---
  # (Same as previous version)
  final_step <- step
  final_idx <- final_step + 1
  
  alphas <- alphas[1:final_idx, , drop = FALSE]
  betas <- betas[1:final_idx, , drop = FALSE]
  if(pa > 0) g_alphas <- g_alphas[1:final_idx, , drop = FALSE] # Gradients at y_alpha
  if(pb > 0) g_betas <- g_betas[1:final_idx, , drop = FALSE]   # Gradients at y_beta
  nllh_results <- nllh_results[1:final_idx]
  
  final_alpha <- if (pa > 0) alphas[final_idx, ] else numeric(0)
  final_beta <- if (pb > 0) betas[final_idx, ] else numeric(0)
  
  if (final_step == 0 && max.step > 0) {
    cli::cli_warn("Optimization loop did not complete any steps.")
    return(list(alpha = alpha.start, beta = beta.start, step = 0,
                alphas = matrix(alpha.start, nrow=1), betas = matrix(beta.start, nrow=1),
                grad_alphas = if(pa > 0) matrix(NA_real_, nrow=1, ncol=pa) else matrix(NA_real_,0,0),
                grad_betas = if(pb > 0) matrix(NA_real_, nrow=1, ncol=pb) else matrix(NA_real_,0,0),
                nllh_results = nllh_results[1]))
  }
  
  return(list(alpha = final_alpha,
              beta = final_beta,
              step = final_step,
              alphas = alphas,
              betas = betas,
              grad_alphas = g_alphas, # Note: Gradients are at momentum point 'y'
              grad_betas = g_betas,   # Note: Gradients are at momentum point 'y'
              nllh_results = nllh_results))
}



# L_alpha <- L(alpha, beta, va, vb, x, y, prob_fun, 
#                      opt = "alpha")
# 
# L_beta <- L(alpha, beta, va, vb, x, y, prob_fun, 
#                     opt = "beta")
# 
# step_size_alpha <- 1/L_alpha
# step_size_beta <- 1/L_beta
# 
# cli::cli_alert("step_size_alpha: {step_size_alpha} | step_size_beta: {step_size_beta}")