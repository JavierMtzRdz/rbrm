# Assuming prob_fun is defined elsewhere and takes logrr and logop as input
# and returns a matrix with two columns: p0 and p1.

# Function to calculate the gradient of the log-likelihood with respect to alpha
# grad_nll_alpha <- function(alpha, beta, X, Y, va, vb, prob_fun) {
#   logrr <- va %*% alpha
#   logop <- vb %*% beta
#   ps <- prob_fun(logrr, logop)
#   p0 <- ps$p0
#   p1 <- ps$p1
#   n <- length(X)
#   grad_alpha <- rep(0, length(alpha))
#   for (i in 1:n) {
#     Ai <- exp(logrr[i, 1]) # Assuming logrr is an n x 1 matrix
#     Wi <- va[i, ]           # Assuming va is an n x length(alpha) matrix
#     term_numerator <- (Y[i] - p0[i]) * (X[i] == 0) - (Y[i] - p1[i]) * (X[i] == 1)
#     term_denominator <- p0[i] * (1 + Ai) - 2
#     grad_alpha <- grad_alpha + (term_numerator / term_denominator) * Wi
#   }
#   return(-grad_alpha)
# }

grad_nll_alpha <- function(alpha, beta, X, Y, va, vb, prob_fun) {
  # Compute the linear predictors for the two models:
  logrr <- as.vector(va %*% alpha)  # theta
  logop <- as.vector(vb %*% beta)     # phi
  
  # Get p0 and p1 from the provided probability function
  ps <- prob_fun(logrr, logop)
  p0 <- ps$p0
  p1 <- ps$p1
  
  # Define A = exp(theta) and B = exp(phi)
  A <- exp(logrr)
  B <- exp(logop)
  
  # Compute the derivative of p0 with respect to theta.
  # This is derived from the closed-form expression for p0.
  dp0_dtheta <- - A * (B * p0 * (1 - p0) + p0^2) / 
    (2 * A * p0 + B * (1 - A * p0) + A * B * (1 - p0))
  
  # Initialize the score vector for each observation
  score <- numeric(length(Y))
  
  # Indices for untreated (X == 0) and treated (X == 1)
  idx0 <- which(X == 0)
  idx1 <- which(X == 1)
  
  # For observations with X == 0:
  # The derivative contribution is (Y/p0 - (1-Y)/(1-p0)) * dp0_dtheta
  score[idx0] <- (Y[idx0] / p0[idx0] - (1 - Y[idx0]) / (1 - p0[idx0])) * 
    dp0_dtheta[idx0]
  
  # For observations with X == 1:
  # Here, note that p1 = p0 * exp(theta) so the derivative is:
  # e^(theta) * (dp0_dtheta + p0)
  score[idx1] <- (Y[idx1] / p1[idx1] - (1 - Y[idx1]) / (1 - p1[idx1])) * 
    ((dp0_dtheta[idx1] + p0[idx1]) * A[idx1])
  
  # The gradient with respect to alpha is the negative sum over observations of:
  # score * the covariate vector (each row of va)
  grad_alpha <- - t(va) %*% score
  return(as.vector(grad_alpha))
}


# Function to calculate the gradient of the log-likelihood with respect to beta
grad_nll_beta <- function(alpha, beta, X, Y, va, vb, prob_fun){
  logrr <- va %*% alpha
  logop <- vb %*% beta
  ps <- prob_fun(logrr, logop)
  p0 <- ps$p0
  p1 <- ps$p1
  
  exp_logrr <- exp(logrr)
  exp_logop <- exp(logop)
  
  # Calculate dp0/dlogop (derivative of p0 with respect to logop)
  term1 <- -(exp_logrr + 1) * exp_logop
  term2 <- sqrt(exp_logop^2 * (exp_logrr + 1)^2 + 4 * exp_logrr * exp_logop * (1 - exp_logop))
  dp0_dlogop <- (
    -(exp_logrr + 1) * exp_logop +
      (exp_logop * (exp_logrr + 1)^2 + 2 * exp_logrr * (1 - 2 * exp_logop)) /
      term2 * exp_logop
  ) / (2 * exp_logrr * (1 - exp_logop)) +
    (term1 + term2) / (2 * exp_logrr * (1 - exp_logop)^2) * exp_logop
  
  
  # Calculate dp1/dlogop
  dp1_dlogop <- dp0_dlogop * exp_logrr
  
  # Calculate the gradient
  dL_dp0 <- (Y * (1 - X)) / p0 - ((1 - Y) * (1 - X)) / (1 - p0)
  dL_dp1 <- (Y * X) / p1 - ((1 - Y) * X) / (1 - p1)
  
  grad <- t(vb) %*% (dL_dp0 * dp0_dlogop + dL_dp1 * dp1_dlogop)
  
  return(-as.vector(grad))
}


# alpha_init <- rep(.2, ncol(va))
# beta_init <- rep(-1, ncol(va))
# 
# (gradient_result <- grad_nll_alpha(alpha_init, beta_init, x, y, va, vb, getProbRR.org))
# 
# (gradient <- numDeriv::grad(function(.x){nllh(.x, beta_init, va, vb, x, y,
#                                               prob_fun = getProbRR.org)},
#                             alpha_init, method = "Richardson"))
# 
# 
# (gradient_result <- grad_nll_beta(alpha_init, beta_init, x, y, va, vb, getProbRR.org))
# 
# (gradient <- numDeriv::grad(function(.x){nllh(alpha_init, .x, va, vb, x, y,
#                                               prob_fun = getProbRR.org)},
#                             beta_init, method = "Richardson"))




##
## Gradient of the NEGATIVE log-likelihood for the *alternative* Pozza et al. specification
## p0 = (M - D) / (2 * C * A),  p1 = p0 * A,
## where M = 1 + C(1 + A),  D = sqrt(M^2 - 4 C^2 A),
##       A = exp(theta), C = exp(phiA).
##
ll_gradient_pozza_alt <- function(alpha, beta, X, Y, va, vb, prob_fun) {
  
  ## 1) Compute linear predictors:
  ##    theta = alpha^T * W,   phiA = beta^T * Z
  theta <- as.vector(va %*% alpha)   # "logRR" part
  phiA  <- as.vector(vb %*% beta)    # "alternative nuisance" part
  
  ## 2) Convert to exponentiated terms for convenience
  A <- exp(theta)    # A = e^theta
  C <- exp(phiA)     # C = e^phiA
  
  ## 3) Define M and D for the closed-form expression of p0
  M <- 1 + C * (1 + A)                # M = 1 + C(1 + A)
  D <- sqrt(M^2 - 4 * C^2 * A)        # D = sqrt(M^2 - 4 C^2 A)
  
  ## 4) Probabilities:
  ##    p0 = (M - D)/(2 C A),  p1 = p0 * A
  # p0 <- (M - D) / (2 * C * A)
  # p1 <- p0 * A
  # logrr <- va %*% alpha
  # logop <- vb %*% beta
  ps <- prob_fun(theta, phiA)
  p0 <- ps$p0
  p1 <- ps$p1
  
  
  ## 5) Derivatives of p0 wrt theta (alpha) and phiA (beta).
  ##
  ##    p0(theta, phiA) = (M - D) / (2 C A)
  ##    with M = 1 + C(1 + A),  D = sqrt(M^2 - 4 C^2 A).
  ##
  ##    Let us define partials carefully:
  ##    - partial(1/(C A))/partial theta = - (1/(C A))
  ##    - partial(M)/partial theta = C*A
  ##    - partial(D)/partial theta = [C A (M - 2 C)] / D
  ##
  ##    Final result for dp0/dtheta is:
  ##        1/2 * [ - (M - D)/(C A) + (C A - C A(M - 2C)/D)/(C A ) ]
  ##     => 1/2 * [ 1 - (M - D)/(C A) - (M - 2C)/D ]
  ##
  dp0_dtheta <- 0.5 * (
    1 -
      (M - D) / (C * A) -    # subtract 2 p0 if you like, since (M-D)/(C A) = 2 p0
      (M - 2 * C) / D
  )
  
  ##
  ##    For dp0/dphiA:
  ##    - partial(1/(C A))/partial phiA = - 1/(C A)
  ##    - partial(M)/partial phiA = (1 + A)*C
  ##    - partial(D)/partial phiA = [ C ( M(1 + A) - 4 C A ) ] / D
  ##
  ##    => dp0/dphiA = 1/2 * [
  ##         - (M - D)/(C A)
  ##         + (1/(C A)) * C { (1 + A) - [ M(1 + A) - 4 C A ]/D }
  ##       ]
  ##    => 1/2 * [
  ##         - (M - D)/(C A)
  ##         + (1/A) { (1 + A) - [M(1 + A) - 4 C A]/D }
  ##       ]
  ##
  term1 <- -(M - D)/(C * A)
  term2 <- (1 + A)/A - (1/A) * ( (M * (1 + A) - 4 * C * A) / D )
  dp0_dphi <- 0.5 * (term1 + term2)
  
  ## 6) p1 = p0 * A => derivatives by product rule:
  ##    dp1/dtheta = A * dp0/dtheta + p0 * dA/dtheta = A * dp0/dtheta + p0 * A
  ##                = A [ dp0/dtheta + p0 ]
  ##    dp1/dphi   = A * dp0_dphi
  ##
  dp1_dtheta <- A * (dp0_dtheta + p0)
  dp1_dphi   <- A * dp0_dphi
  
  ## 7) For a Bernoulli log-likelihood contribution:
  ##    log f_i =  X_i [Y_i log p1_i + (1 - Y_i) log(1 - p1_i)]
  ##             + (1 - X_i)[Y_i log p0_i + (1 - Y_i) log(1 - p0_i)]
  ##
  ##    derivative wrt p0 or p1 => factor = Y/p - (1-Y)/(1-p).
  ##    We'll compute factor0 for p0 and factor1 for p1, then do chain rule.
  ##
  factor0 <- Y / p0 - (1 - Y) / (1 - p0)   # used if X=0
  factor1 <- Y / p1 - (1 - Y) / (1 - p1)   # used if X=1
  
  dL_dtheta <- (1 - X) * dp0_dtheta * factor0  +  X * dp1_dtheta * factor1
  dL_dphi   <- (1 - X) * dp0_dphi   * factor0  +  X * dp1_dphi   * factor1
  
  ## 8) Gradient of the *negative* log-likelihood is the negative of the sum of dL
  grad_alpha <- - t(dL_dtheta) %*% (va)
  grad_beta  <- - t(dL_dphi) %*% (vb)
  
  ## Return as list
  list(grad_alpha = as.vector(grad_alpha), grad_beta = as.vector(grad_beta))
}


# ll_gradient_pzz <- function(gamma, beta, x, y, va, vb, prob_fun = getProbRR.alt) {
#   logrr <- va %*% gamma
#   logop <- vb %*% beta
#   ps <- prob_fun(logrr, logop)
#   p0 <- ps$p0
#   p1 <- ps$p1
#   
#   
#   derivative_p0_dlogrr <- function(eta1, eta2)
#   {
#     dp0.eta1 <-  ((exp(-eta2-eta1)*(exp(eta2+eta1)-(2*((exp(eta1)+1)*exp(eta2)+1)*exp(eta2+eta1)-4*exp(2*eta2+eta1))/
#                                       (2*sqrt(((exp(eta1)+1)*exp(eta2)+1)^2-4*exp(2*eta2+eta1)))))/2-
#                     (exp(-eta2-eta1)*(-sqrt(((exp(eta1)+1)*exp(eta2)+1)^2
#                                             -4*exp(2*eta2+eta1))+(exp(eta1)+1)*exp(eta2)+1))/2)
#     return(dp0.eta1)
#   }
#   
#   derivative_p0_dlogop <- function(eta1,eta2)
#   {
#     dp0.eta2 <- ((exp(-eta2-eta1)*((exp(eta1)+1)*exp(eta2)-(2*(exp(eta1)+1)*exp(eta2)*((exp(eta1)+1)*exp(eta2)+1)-8*exp(2*eta2+eta1))/
#                                      (2*sqrt(((exp(eta1)+1)*exp(eta2)+1)^2-4*exp(2*eta2+eta1)))))/2-
#                    (exp(-eta2-eta1)*(-sqrt(((exp(eta1)+1)*exp(eta2)+1)^2-4*exp(2*eta2+eta1))+(exp(eta1)+1)*exp(eta2)+1))/2)
#     return(dp0.eta2)
#   }
#   
#   dp0_dlogrr_val <- derivative_p0_dlogrr(logrr, logop)
#   dp0_dlogop_val <- derivative_p0_dlogop(logrr, logop)
#   
#   d.1.logrr <- (dp0_dlogrr_val * exp(logrr * x) + pi * x)
#   d.1.logop <- (dp0_dlogop_val * exp(logrr * x))
#   
#   #complete the formula 
#   ...
#   
#   grad_alpha <- d.1.logrr
#   grad_beta <- 
#   
#   return(list(grad_alpha = as.vector(-grad_alpha), grad_beta = as.vector(-grad_beta)))
# }



ll_gradient_pzz <- function(alpha, beta, x, y, va, vb, prob_fun = getProbRR.alt) {
  logrr <- va %*% alpha
  logop <- vb %*% beta
  ps <- prob_fun(logrr, logop)
  p0 <- ps$p0
  p1 <- ps$p1
  
  
  derivative_p0_dlogrr <- function(eta1, eta2)
  {
    dp0.eta1 <-  ((exp(-eta2-eta1)*(exp(eta2+eta1)-(2*((exp(eta1)+1)*exp(eta2)+1)*exp(eta2+eta1)-4*exp(2*eta2+eta1))/
                                      (2*sqrt(((exp(eta1)+1)*exp(eta2)+1)^2-4*exp(2*eta2+eta1)))))/2-
                    (exp(-eta2-eta1)*(-sqrt(((exp(eta1)+1)*exp(eta2)+1)^2
                                            -4*exp(2*eta2+eta1))+(exp(eta1)+1)*exp(eta2)+1))/2)
    return(dp0.eta1)
  }
  
  derivative_p0_dlogop <- function(eta1,eta2)
  {
    dp0.eta2 <- ((exp(-eta2-eta1)*((exp(eta1)+1)*exp(eta2)-(2*(exp(eta1)+1)*exp(eta2)*((exp(eta1)+1)*exp(eta2)+1)-8*exp(2*eta2+eta1))/
                                     (2*sqrt(((exp(eta1)+1)*exp(eta2)+1)^2-4*exp(2*eta2+eta1)))))/2-
                   (exp(-eta2-eta1)*(-sqrt(((exp(eta1)+1)*exp(eta2)+1)^2-4*exp(2*eta2+eta1))+(exp(eta1)+1)*exp(eta2)+1))/2)
    return(dp0.eta2)
  }
  
  dp0_dlogrr_val <- derivative_p0_dlogrr(logrr, logop)
  dp0_dlogop_val <- derivative_p0_dlogop(logrr, logop)
  
  grad_alpha <- numeric(length(alpha))
  for (j in 1:length(alpha)) {
    score_alpha_j <- 0
    for (i in 1:length(y)) {
      term1 <- 0
      if (x[i] == 1) {
        term1 <- (y[i] - p1[i]) / (p1[i] * (1 - p1[i])) * (dp0_dlogrr_val[i] * exp(logrr[i]) + p1[i]) * va[i, j]
      }
      term2 <- 0
      if (x[i] == 0) {
        term2 <- (y[i] - p0[i]) / (p0[i] * (1 - p0[i])) * dp0_dlogrr_val[i] * va[i, j]
      }
      score_alpha_j <- score_alpha_j + term1 + term2
    }
    grad_alpha[j] <- score_alpha_j
  }
  
  grad_beta <- numeric(length(beta))
  for (k in 1:length(beta)) {
    score_beta_k <- 0
    for (i in 1:length(y)) {
      term1 <- 0
      if (x[i] == 1) {
        term1 <- (y[i] - p1[i]) / (p1[i] * (1 - p1[i])) * (dp0_dlogop_val[i] * exp(logrr[i])) * vb[i, k]
      }
      term2 <- 0
      if (x[i] == 0) {
        term2 <- (y[i] - p0[i]) / (p0[i] * (1 - p0[i])) * dp0_dlogop_val[i] * vb[i, k]
      }
      score_beta_k <- score_beta_k + term1 + term2
    }
    grad_beta[k] <- score_beta_k
  }
  
  return(list(grad_alpha = as.vector(-grad_alpha), grad_beta = as.vector(-grad_beta)))
}





# alpha_init <- rep(1, ncol(va))
# beta_init <- rep(.4, ncol(va))
# 
# (gradient_result <- ll_gradient_pzz(alpha_init, beta_init, x, y, va, vb, getProbRR.alt)$grad_alpha)
# 
# (gradient_result <- ll_gradient_pozza_alt(alpha_init, beta_init, x, y, va, vb, getProbRR.alt)$grad_alpha)
# 
# (gradient <- numDeriv::grad(function(.x){nllh(.x, beta_init, va, vb, x, y,
#                                               prob_fun = getProbRR.alt)},
#                             alpha_init, method = "simple"))
# 
# 
# (gradient_result <- ll_gradient_pzz(alpha_init, beta_init, x, y, va, vb, getProbRR.alt)$grad_beta)
# (gradient_result <- ll_gradient_pozza_alt(alpha_init, beta_init, x, y, va, vb, getProbRR.alt)$grad_beta)
# 
# (gradient <- numDeriv::grad(function(.x){nllh(alpha_init, .x, va, vb, x, y,
#                                               prob_fun = getProbRR.alt)},
#                             beta_init, method = "Richardson"))





#' @title Gradient of the Negative Log-Likelihood for Risk Ratio Models
#'
#' @description
#' Computes the gradient of the negative log-likelihood function for risk ratio models
#' with respect to the parameters alpha and beta.  This is an internal function
#' and not meant to be called directly by users. It supports both the Richardson
#' et al. (2017) and the alternative (Pozza et al., 2023) parameterizations.
#'
#' @param alpha A numeric vector of coefficients for the treatment model.
#' @param beta A numeric vector of coefficients for the control model.
#' @param X A numeric vector indicating treatment status (1 = treated, 0 = control).
#' @param Y A numeric vector of binary outcomes (1 = event, 0 = no event).
#' @param va A numeric matrix of covariates for the treatment model (design matrix for alpha).
#' @param vb A numeric matrix of covariates for the control model (design matrix for beta).
#' @param prob_fun A function that computes probabilities p0 and p1 based on logrr and logop.  This
#'   function determines the specific parameterization (e.g., Richardson or alternative).
#'   It should return a list with elements `p0`, `p1`, and `class`.
#' @param opt A character string indicating which parameters to compute the gradient for.
#'   Must be either "alpha" or "beta".
#'
#' @return A numeric vector representing the gradient of the negative log-likelihood.
#'
#' @keywords internal
#'
.grad_nll_alpha <- function(alpha, beta, X, Y, va, vb, prob_fun) {
  # Compute the linear predictors for the two models:
  logrr <- as.vector(va %*% alpha)  # theta
  logop <- as.vector(vb %*% beta)    # phi or phi_A
  
  # Get p0 and p1 from the provided probability function
  ps <- prob_fun(logrr, logop)
  p0 <- ps$p0
  p1 <- ps$p1
  
  prob_fun_class <- ps$class
  is_richardson <- grepl("ProbRR\\.org", prob_fun_class)
  is_alternative <- grepl("ProbRR\\.alt", prob_fun_class)
  
  if (is_richardson) {
    # Richardson et al. (2017) specification
    A <- exp(logrr)
    B <- exp(logop)
    
    # Compute the derivative of p0 with respect to theta.
    # This is derived from the closed-form expression for p0.
    dp0_dtheta <- - A * (B * p0 * (1 - p0) + p0^2) /
      (2 * A * p0 + B * (1 - A * p0) + A * B * (1 - p0))
    
    # Initialize the score vector for each observation
    score <- numeric(length(Y))
    
    # Indices for untreated (X == 0) and treated (X == 1)
    idx0 <- which(X == 0)
    idx1 <- which(X == 1)
    
    # For observations with X == 0:
    # The derivative contribution is (Y/p0 - (1-Y)/(1-p0)) * dp0_dtheta
    score[idx0] <- (Y[idx0] / p0[idx0] - (1 - Y[idx0]) / (1 - p0[idx0])) *
      dp0_dtheta[idx0]
    
    # For observations with X == 1:
    # Here, note that p1 = p0 * exp(theta) so the derivative is:
    # e^(theta) * (dp0_dtheta + p0)
    score[idx1] <- (Y[idx1] / p1[idx1] - (1 - Y[idx1]) / (1 - p1[idx1])) *
      ((dp0_dtheta[idx1] + p0[idx1]) * A[idx1])
    
    # The gradient with respect to alpha is the negative sum over observations of:
    # score * the covariate vector (each row of va)
    grad_alpha <- - t(va) %*% score
    return(as.vector(grad_alpha))
    
  } else if (is_alternative) {
    # Pozza et al. (2023) alternative specification
    
    derivative_p0_dlogrr <- function(eta1, eta2)
    {
      dp0.eta1 <-  ((exp(-eta2-eta1)*(exp(eta2+eta1)-(2*((exp(eta1)+1)*exp(eta2)+1)*exp(eta2+eta1)-4*exp(2*eta2+eta1))/
                                        (2*sqrt(((exp(eta1)+1)*exp(eta2)+1)^2-4*exp(2*eta2+eta1)))))/2-
                      (exp(-eta2-eta1)*(-sqrt(((exp(eta1)+1)*exp(eta2)+1)^2
                                              -4*exp(2*eta2+eta1))+(exp(eta1)+1)*exp(eta2)+1))/2)
      return(dp0.eta1)
    }
    
    dp0_dlogrr_val <- derivative_p0_dlogrr(logrr, logop)
    
    grad_alpha <- numeric(length(alpha))
    for (j in 1:length(alpha)) {
      score_alpha_j <- 0
      for (i in 1:length(Y)) {
        term1 <- 0
        if (X[i] == 1) {
          term1 <- (Y[i] - p1[i]) / (p1[i] * (1 - p1[i])) * (dp0_dlogrr_val[i] * exp(logrr[i]) + p1[i]) * va[i, j]
        }
        term2 <- 0
        if (X[i] == 0) {
          term2 <- (Y[i] - p0[i]) / (p0[i] * (1 - p0[i])) * dp0_dlogrr_val[i] * va[i, j]
        }
        score_alpha_j <- score_alpha_j + term1 + term2
      }
      grad_alpha[j] <- score_alpha_j
    }
    return(as.vector(-grad_alpha))
  } else {
    stop("prob_fun not recognized.")
  }
}


#' @rdname dot-grad_nll_alpha
#' @keywords internal
.grad_nll_beta <- function(alpha, beta, X, Y, va, vb, prob_fun) {
  
  # Compute the linear predictors for the two models:
  logrr <- as.vector(va %*% alpha)  # theta
  logop <- as.vector(vb %*% beta)    # phi or phi_A
  
  # Get p0 and p1 from the provided probability function
  ps <- prob_fun(logrr, logop)
  p0 <- ps$p0
  p1 <- ps$p1
  
  prob_fun_class <- ps$class
  is_richardson <- grepl("ProbRR\\.org", prob_fun_class)
  is_alternative <- grepl("ProbRR\\.alt", prob_fun_class)
  if (is_richardson) {
    exp_logrr <- exp(logrr)
    exp_logop <- exp(logop)
    
    # Calculate dp0/dlogop (derivative of p0 with respect to logop)
    term1 <- -(exp_logrr + 1) * exp_logop
    term2 <- sqrt(exp_logop^2 * (exp_logrr + 1)^2 + 4 * exp_logrr * exp_logop * (1 - exp_logop))
    dp0_dlogop <- (
      -(exp_logrr + 1) * exp_logop +
        (exp_logop * (exp_logrr + 1)^2 + 2 * exp_logrr * (1 - 2 * exp_logop)) /
        term2 * exp_logop
    ) / (2 * exp_logrr * (1 - exp_logop)) +
      (term1 + term2) / (2 * exp_logrr * (1 - exp_logop)^2) * exp_logop
    
    
    # Calculate dp1/dlogop
    dp1_dlogop <- dp0_dlogop * exp_logrr
    
    # Calculate the gradient
    dL_dp0 <- (Y * (1 - X)) / p0 - ((1 - Y) * (1 - X)) / (1 - p0)
    dL_dp1 <- (Y * X) / p1 - ((1 - Y) * X) / (1 - p1)
    
    grad <- t(vb) %*% (dL_dp0 * dp0_dlogop + dL_dp1 * dp1_dlogop)
    
    return(-as.vector(grad))
  } else if (is_alternative){
    derivative_p0_dlogop <- function(eta1,eta2)
    {
      dp0.eta2 <- ((exp(-eta2-eta1)*((exp(eta1)+1)*exp(eta2)-(2*(exp(eta1)+1)*exp(eta2)*((exp(eta1)+1)*exp(eta2)+1)-8*exp(2*eta2+eta1))/
                                       (2*sqrt(((exp(eta1)+1)*exp(eta2)+1)^2-4*exp(2*eta2+eta1)))))/2-
                     (exp(-eta2-eta1)*(-sqrt(((exp(eta1)+1)*exp(eta2)+1)^2-4*exp(2*eta2+eta1))+(exp(eta1)+1)*exp(eta2)+1))/2)
      return(dp0.eta2)
    }
    
    dp0_dlogop_val <- derivative_p0_dlogop(logrr, logop)
    
    grad_beta <- numeric(length(beta))
    for (k in 1:length(beta)) {
      score_beta_k <- 0
      for (i in 1:length(Y)) {
        term1 <- 0
        if (X[i] == 1) {
          term1 <- (Y[i] - p1[i]) / (p1[i] * (1 - p1[i])) * (dp0_dlogop_val[i] * exp(logrr[i])) * vb[i, k]
        }
        term2 <- 0
        if (X[i] == 0) {
          term2 <- (Y[i] - p0[i]) / (p0[i] * (1 - p0[i])) * dp0_dlogop_val[i] * vb[i, k]
        }
        score_beta_k <- score_beta_k + term1 + term2
      }
      grad_beta[k] <- score_beta_k
    }
    return(as.vector(-grad_beta))
  } else {
    stop("prob_fun not recognized.")
  }
}



#' @title Gradient of the Negative Log-Likelihood (Unified Interface)
#'
#' @description This function provides a unified interface to compute the gradient,
#'   dispatching to either `.grad_nll_alpha` or `.grad_nll_beta` based on the `opt` argument.
#'   It's the primary exported function for gradient calculation.
#'
#' @inheritParams .grad_nll_alpha
#' @export
grad_nll <- function(alpha, beta, X, Y, va, vb, prob_fun, opt = c("alpha", "beta")) {
  opt <- match.arg(opt)
  
  if (opt == "alpha") {
    return(.grad_nll_alpha(alpha, beta, X, Y, va, vb, prob_fun))
  } else if (opt == "beta") {
    return(.grad_nll_beta(alpha, beta, X, Y, va, vb, prob_fun))
  }
}




# (gradient_result <- grad_nll(alpha_init, beta_init, x, y, va, vb, getProbRR.org, opt = 'alpha'))
# 
# (gradient_result <- grad_nll_alpha(alpha_init, beta_init, x, y, va, vb, getProbRR.org))
# 
# (gradient <- numDeriv::grad(function(.x){nllh(.x, beta_init, va, vb, x, y,
#                                               prob_fun = getProbRR.org)},
#                             alpha_init, method = "Richardson"))
# 
# (gradient_result <- grad_nll(alpha_init, beta_init, x, y, va, vb, getProbRR.org, opt = 'beta'))
# 
# (gradient_result <- grad_nll_beta(alpha_init, beta_init, x, y, va, vb, getProbRR.org))
# 
# 
# (gradient <- numDeriv::grad(function(.x){nllh(alpha_init, .x, va, vb, x, y,
#                                               prob_fun = getProbRR.org)},
#                             beta_init, method = "Richardson"))




# (gradient_result <- grad_nll(alpha_init, beta_init, x, y, va, vb,
#                              purrr::partial(getProbRR.org,
#                                             clipping = T), opt = 'alpha'))
# 
# (gradient_result <- ll_gradient_pzz(alpha_init, beta_init, x, y, va, vb, getProbRR.alt)$grad_alpha)
# 
# (gradient <- numDeriv::grad(function(.x){nllh(.x, beta_init, va, vb, x, y,
#                                               prob_fun = getProbRR.alt)},
#                             alpha_init, method = "Richardson"))
# 
# (gradient_result <- grad_nll(alpha_init, beta_init, x, y, va, vb, getProbRR.alt, opt = 'beta'))
# (gradient_result <- ll_gradient_pzz(alpha_init, beta_init, x, y, va, vb, getProbRR.alt)$grad_beta)
# 
# (gradient <- numDeriv::grad(function(.x){nllh(alpha_init, .x, va, vb, x, y,
#                                               prob_fun = getProbRR.alt)},
#                             beta_init, method = "Richardson"))