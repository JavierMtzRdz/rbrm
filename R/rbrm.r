#' Soft Thresholding Operator
#' 
#' Applies the soft thresholding operator to the input.
#' 
#' @param x A numeric vector to apply the soft threshold.
#' @param lambda A positive numeric value representing the threshold parameter.
#' @return A numeric vector where the soft thresholding has been applied.
#' @examples
#' #soft_thres(c(3, -1.5, 0.2), 0.5)
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
nllh.alpha <- function(alpha, beta, va, vb, x, y) {
  logrr <- (va %*% alpha)
  logop <- (vb %*% beta)
  p0 <- brm::getProbRR(logrr, logop)[, 1]
  p1 <- brm::getProbRR(logrr, logop)[, 2]
  
  # Clipping probabilities
  p0 <- pmin(pmax(p0, 1e-15), 1 - 1e-15)
  p1 <- pmin(pmax(p1, 1e-15), 1 - 1e-15)
  
  return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) + (y[x == 0]) * log(p0[x == 0])) -
           sum((1 - y[x == 1]) * log(1 - p1[x == 1]) + (y[x == 1]) * log(p1[x == 1])))
}





#' Negative Log-Likelihood for Beta
#' 
#' Computes the negative log-likelihood for the beta coefficients.
#' 
#' @param beta Coefficients.
#' @return The negative log-likelihood for the given beta.
#' @examples
#' #beta <- c(1, 2, 3)
#' #nllh.beta(beta)
nllh.beta <- function(beta, alpha, va, vb, x, y) {
  logrr <- (va %*% alpha)
  logop <- (vb %*% beta)
  p0 <- brm::getProbRR(logrr, logop)[, 1]
  p1 <- brm::getProbRR(logrr, logop)[, 2]
  
  # Clipping probabilities
  p0 <- pmin(pmax(p0, 1e-15), 1 - 1e-15)
  p1 <- pmin(pmax(p1, 1e-15), 1 - 1e-15)
  
  return(-sum((1 - y[x == 0]) * log1p(-p0[x == 0]) + (y[x == 0]) * log(p0[x == 0])) -
           sum((1 - y[x == 1]) * log1p(-p1[x == 1]) + (y[x == 1]) * log(p1[x == 1])))
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
proximal.gd.alpha.fista <- function(alpha, step_size, lambda, t_old, last_alpha,
                                    intercept,
                                    beta, va, vb, x, y) {
  gradient <- numDeriv::grad(function(.x){nllh.alpha(.x, beta, va, vb, x, y)}, alpha, method = "simple")
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
proximal.gd.beta.fista <- function(beta, step_size, lambda, t_old, last_beta,
                                   intercept,
                                   alpha, va, vb, x, y) {
  gradient <- numDeriv::grad(function(.x) {nllh.beta(.x, alpha, va, vb, x, y)}, beta, method = "simple")
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
penalized.neg.log.likelihood <- function(pars, pa, pb, intercept, lambda,
                                         alpha, beta, va, vb, x, y) {
  alpha <- pars[1:pa]
  beta <- pars[(pa + 1):(pa + pb)]
  logrr <- (va %*% alpha)
  logop <- (vb %*% beta)
  p0 <- brm::getProbRR(logrr, logop)[, 1]
  p1 <- brm::getProbRR(logrr, logop)[, 2]
  
  # Clipping probabilities
  p0 <- pmin(pmax(p0, 1e-15), 1 - 1e-15)
  p1 <- pmin(pmax(p1, 1e-15), 1 - 1e-15)
  
  unpenalized.nllh <- (-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) + (y[x == 0]) * log(p0[x == 0])) -
                         sum((1 - y[x == 1]) * log(1 - p1[x == 1]) + (y[x == 1]) * log(p1[x == 1])))
  
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
rbrm <- function(va, vb, y, x, alpha.start = NULL, beta.start = NULL,
                 max.step = 3000, thres = 1e-04, lambda = 0,
                 lr.alpha = 0.06, lr.beta = 0.02,
                 intercept = TRUE, early_stopping_rounds = 10) {
  # 
  # va <- v; vb <- v; alpha.start = NULL; beta.start = NULL;
  # max.step = 3000; thres = 1e-04; lambda = 0;
  # lr.alpha = 0.06; lr.beta = 0.02;
  # intercept = TRUE; early_stopping_rounds = 10

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

        grad_alpha <- numDeriv::grad(function(.x) {nllh.alpha(.x, beta, va, vb, x, y)}, y_alpha, method = "simple")
        grad_beta <- numDeriv::grad(function(.x) {nllh.beta(.x, alpha, va, vb, x, y)}, y_beta, method = "simple")
        
        if ((max(abs(grad_alpha)) < thres) & (max(abs(grad_beta)) < thres)) {
            break
        }

        # Update alpha and beta for the next iteration
        alpha <- alpha_new
        beta <- beta_new
    }
    opt <- list(
        point.est = c(alpha, beta), convergence = (step < max.step),
        value = penalized.neg.log.likelihood(c(alpha, beta), 
                                             pa, pb, intercept, lambda,
                                             alpha, beta, va, vb, x, y),
        step = step
    )

    return(opt)
}
