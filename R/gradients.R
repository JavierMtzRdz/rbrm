
# derivative of pi0 with respect to theta
dp0_theta <- function(theta, phi, 
                      ps_spe = c("Richardson", "Pozza"), ep = 1e-8){
  
  ps_spe <- rlang::arg_match(ps_spe)
  
  if(ps_spe == "Richardson"){
  dp0.theta <- ( -exp(phi - theta) / (2 * (exp(phi) - 1)) + exp(phi) /
                   ((exp(phi) - 1) * sqrt(4 * exp(phi + theta) + 
                   (exp(theta) - 1)^2 *  exp(2 * phi))) +
                   (exp( - theta) * (1 - exp(theta)) *
                    exp(2 * phi)) / (2 * (exp(phi) - 1)*
                    sqrt(4 * exp(phi + theta) + (exp(theta) - 1)^2 * 
                           exp(2*phi))) )
  
  # extension for continuity
  dp0.theta[abs(phi)< ep] <- (- exp(theta[abs(phi)< ep]) /
                                (exp(theta[abs(phi)< ep]) + 1) ^ 2)
  return(dp0.theta)}
  
  if(ps_spe == "Pozza"){
    dp0.theta <-  ((exp(-phi-theta)*(exp(phi+theta)-
                 (2*((exp(theta)+1)*exp(phi)+1)*exp(phi+theta)-
                 4*exp(2*phi+theta))/
                 (2*sqrt(((exp(theta)+1)*exp(phi)+1)^2-
                 4*exp(2*phi+theta)))))/2-(exp(-phi-theta)*
                  (-sqrt(((exp(theta)+1)*exp(phi)+1)^2-
                  4*exp(2*phi+theta))+(exp(theta)+1)*exp(phi)+1))/2)
    return(dp0.theta)}
}
# derivative of pi0 with respect to phi
# written in a different way compared to thesis
dp0_phi <- function(theta, phi, 
                    ps_spe = c("Richardson", "Pozza"), ep = 1e-8){
  
  ps_spe <- rlang::arg_match(ps_spe)
  
  if(ps_spe == "Richardson"){
  dp0.phi <- ( - ((exp(theta) + 1) * exp(phi)) / 
              (2 * exp(theta) * (exp(phi) - 1) ^ 2) + exp(phi) / 
              ((exp(phi) - 1) ^ 2 * sqrt(4 * exp(phi + theta) + 
              (exp(theta) - 1) ^ 2 * exp(2 * phi))) + 
              (exp( - theta) * (exp(2 * theta) + 1) * exp(2 * phi)) /
              (2 * (exp(phi) - 1) ^ 2 * sqrt(4 * exp(phi + theta) +
              (exp(theta) - 1) ^ 2 * exp(2 * phi))))
  
  # extension for continuity
  dp0.phi[abs(phi)< ep] <- (exp(theta[abs(phi)< ep]) /
                                (exp(theta[abs(phi)< ep]) + 1) ^ 3)
  return(dp0.phi)}
  
  if(ps_spe == "Pozza"){
    dp0.phi <- ((exp(-phi-theta)*((exp(theta)+1)*exp(phi)-(2*(exp(theta)+1)*exp(phi)*((exp(theta)+1)*exp(phi)+1)-8*exp(2*phi+theta))/
                                     (2*sqrt(((exp(theta)+1)*exp(phi)+1)^2-4*exp(2*phi+theta)))))/2-
                   (exp(-phi-theta)*(-sqrt(((exp(theta)+1)*exp(phi)+1)^2-4*exp(2*phi+theta))+(exp(theta)+1)*exp(phi)+1))/2)
    return(dp0.phi)}
}




#' Calculate Analytical Gradient of User's nllh Function (RR Model) - Refined Check
#' Calls the stable derivative helper function.
#'
#' @param alpha Vector of parameters for theta.
#' @param beta Vector of parameters for phi.
#' @param va Design matrix for alpha parameters (W in paper). Rows are observations.
#' @param vb Design matrix for beta parameters (Z in paper). Rows are observations.
#' @param x Vector of binary exposures/treatments (A in paper).
#' @param y Vector of binary outcomes (Y in paper).
#' @param prob_fun Function that takes theta (theta), phi (phi) and returns list(p0, p1).
#' @param opt Character string specifying which gradient to return ("alpha", "beta", or "both").
#'
#' @return Analytical gradient vector(s) of the negative log-likelihood defined by user's nllh.
#' @export
grad_nll <- function(alpha, beta, y, x, va, vb, 
                     prob_fun, opt = c("both", "alpha", "beta"),
                     method = c("analytical", "numerical")) {
  # Argument matching
  method <- rlang::arg_match(method)
  opt <- rlang::arg_match(opt)
  
  n <- length(y)
  pa <- length(alpha)
  pb <- length(beta)
  
  # Calculate theta and phi
  theta <- as.vector(va %*% alpha)
  phi <- as.vector(vb %*% beta)
  
  ps <- prob_fun(theta, phi)
  
  p0 <- ps$p0
  p1 <- ps$p1
  
  # Derivatives of log-likelihood_i w.r.t p1_i and p0_i
  dllh_dp1 <- (y * x) / p1 - ((1 - y) * x) / (1 - p1)
  dllh_dp0 <- (y * (1 - x)) / p0 - ((1 - y) * (1 - x)) / (1 - p0)
  
  if (opt != "beta"){
    if(method == "analytical") dp0_dtheta <- dp0_theta(theta, phi,
    ps_spe = class(ps))
    # if(method == "analytical") dp0_dtheta <- -(1 - p0)/(1 - p0 + 1 - p1)
    if(method == "numerical") dp0_dtheta <- numDeriv::grad(
      func = (function(theta_val)return(prob_fun(theta_val, phi)$p0)),
      x = theta)
    dp1_dtheta <- (p0 + dp0_dtheta)*exp(theta)
    grad_alpha_sum <- numeric(pa)
    inner_alpha <- (dllh_dp1*dp1_dtheta + dllh_dp0*dp0_dtheta)
    grad_alpha <- -(t(va)%*%inner_alpha)/n
  }
  if (opt != "alpha"){
    
    if(method == "analytical") dp0_dphi <- dp0_phi(theta, phi,
                                                   ps_spe = class(ps))
    # if(method == "analytical") dp0_dphi <- (1 - p0) * (1 - p1)/(1 - p0 + 1 - p1)
    if(method == "numerical") dp0_dphi <- numDeriv::grad(
      func = (function(phi_val) return(prob_fun(theta, phi_val)$p0)),
      x = phi)

    dp1_dphi <- dp0_dphi * exp(theta)
    grad_beta_sum <- numeric(pb)
    inner_beta <- (dllh_dp1*dp1_dphi + dllh_dp0*dp0_dphi)
    grad_beta <- -t(inner_beta%*%vb)/n

    # neg.log.likelihood.beta = function(beta){
    #   p0p1 = brm::getProbRR(va %*% alpha, vb %*% beta)
    #   p0 = p0p1[, 1];   p1 = p0p1[, 2]
    #   weights = rep(1, length(y))
    #   
    #   return((-sum((1-y[x==0])*log(1-p0[x==0])*weights[x==0] +
    #                  (y[x==0])*log(p0[x==0])*weights[x==0]) -
    #             sum((1-y[x==1])*log(1-p1[x==1])*weights[x==1] +
    #                   (y[x==1])*log(p1[x==1])*weights[x==1]))/length(y))  
    # }
    # grad_beta2 <- numDeriv::grad(function(.x) {neg.log.likelihood.beta(.x)},
    #                             beta)
  }
  
  if (opt == "alpha") return(grad_alpha)
  if (opt == "beta") return(grad_beta)
  if (opt == "both") return(list(grad_alpha = grad_alpha, 
                                 grad_beta = grad_beta))
  
}

#' @export
grad_nll_k <- function(alpha, beta, y, x, va, vb,
                       prob_fun,
                       opt = c("alpha", "beta"), 
                       k_index,              
                       method = c("numerical", "analytical"))  {
  # Argument matching
  opt <- rlang::arg_match(opt)
  method <- rlang::arg_match(method)
  
  n <- length(y)
  pa <- length(alpha)
  pb <- length(beta)
  
  # Calculate theta and phi
  theta <- as.vector(va %*% alpha)
  phi <- as.vector(vb %*% beta)
  
  ps <- prob_fun(theta, phi)
  
  p0 <- ps$p0
  p1 <- ps$p1
  
  
  # Derivatives of log-likelihood_i w.r.t p1_i and p0_i
  dllh_dp1 <- (y * x) / p1 - ((1 - y) * x) / (1 - p1)
  dllh_dp0 <- (y * (1 - x)) / p0 - ((1 - y) * (1 - x)) / (1 - p0)
  
  # Handle cases where terms are not applicable to avoid NaN 
  # dllh_dp1[x == 0] <- 0
  # dllh_dp0[x == 1] <- 0
  
  
  # Initialize the required gradient component
  grad_kth_component <- NA
  
  if (opt == "alpha") {
    dp0_dtheta_val <- NULL
    if (method == "analytical") {
      dp0_dtheta_val <- dp0_theta(theta, phi, ps_spe = class(ps))
    } else if (method == "numerical") {
      get_p0_for_theta <- function(theta_val) return(prob_fun(theta_val, phi)$p0)
      jacobian_p0_theta <- numDeriv::jacobian(func = get_p0_for_theta, x = theta)
      dp0_dtheta_val <- diag(jacobian_p0_theta)
    }
    
    dp1_dtheta_val <- (p0 + dp0_dtheta_val) * exp(theta)
    inner_alpha_terms <- (dllh_dp1 * dp1_dtheta_val + dllh_dp0 * dp0_dtheta_val)
    grad_kth_component <- -sum(va[, k_index] * inner_alpha_terms) / n
    
  } else if (opt == "beta") {
    dp0_dphi_val <- NULL
    if (method == "analytical") {
      dp0_dphi_val <- dp0_phi(theta, phi, ps_spe = class(ps))
    } else if (method == "numerical") {
      get_p0_for_phi <- function(phi_val) return(prob_fun(theta, phi_val)$p0)
      jacobian_p0_phi <- numDeriv::jacobian(func = get_p0_for_phi, x = phi)
      dp0_dphi_val <- diag(jacobian_p0_phi)
    }
    
    dp1_dphi_val <- dp0_dphi_val * exp(theta) 
    inner_beta_terms <- (dllh_dp1 * dp1_dphi_val + dllh_dp0 * dp0_dphi_val)
    grad_kth_component <- -sum(vb[, k_index] * inner_beta_terms) / n
  }
  
  return(grad_kth_component)
}


#' @export
hessian_or = function(y, x, va, vb, alpha.ml, beta.ml, weights) {
  # calculating the Hessian using the second derivative have to do so
  # because under mis-specification of models Hessian no longer equals the
  # square of the first order derivatives
  
  p0p1 = brm::getProbRR(va %*% alpha.ml, vb %*% beta.ml)
  # p0p1 = cbind(p0, p1): n * 2 matrix
  p0 = p0p1[, 1]
  p1 = p0p1[, 2]
  n = nrow(va)
  pA = p0
  pA[x == 1] = p1[x == 1]
  
  
  ### Building blocks
  
  dpsi0.by.dtheta = -(1 - p0)/(1 - p0 + 1 - p1)
  dpsi0.by.dphi = (1 - p0) * (1 - p1)/(1 - p0 + 1 - p1)
  
  dtheta.by.dalpha = va
  dphi.by.dbeta = vb
  
  dl.by.dpsi0 = (y - pA)/(1 - pA)
  d2l.by.dpsi0.2 = (y - 1) * pA/((1 - pA)^2)
  
  
  
  ###### d2l.by.dalpha.2
  
  d2psi0.by.dtheta.2 = ((p0 - p1) * dpsi0.by.dtheta - (1 - p0) * p1)/((1 - 
                                                                         p0 + 1 - p1)^2)
  
  d2l.by.dtheta.2 = d2l.by.dpsi0.2 * (dpsi0.by.dtheta + x)^2 + dl.by.dpsi0 * 
    d2psi0.by.dtheta.2
  
  d2l.by.dalpha.2 = t(dtheta.by.dalpha * d2l.by.dtheta.2 * weights) %*% 
    dtheta.by.dalpha
  
  
  ###### d2l.by.dalpha.dbeta
  
  d2psi0.by.dtheta.dphi = (1 - p0) * (1 - p1) * (p0 - p1)/(1 - p0 + 1 - 
                                                             p1)^3
  
  d2l.by.dtheta.dphi = d2l.by.dpsi0.2 * (dpsi0.by.dtheta + x) * dpsi0.by.dphi + 
    dl.by.dpsi0 * d2psi0.by.dtheta.dphi
  
  d2l.by.dalpha.dbeta = t(dtheta.by.dalpha * d2l.by.dtheta.dphi * weights) %*% 
    dphi.by.dbeta
  d2l.by.dbeta.dalpha = t(d2l.by.dalpha.dbeta)
  # d2l.by.dalpha.dbeta is symmetric itself if (because) va=vb
  
  
  #### d2l.by.dbeta2
  
  d2psi0.by.dphi.2 = (-(p0 * (1 - p1)^2 + p1 * (1 - p0)^2)/(1 - p0 + 1 - 
                                                              p1)^2) * dpsi0.by.dphi
  
  d2l.by.dphi.2 = d2l.by.dpsi0.2 * (dpsi0.by.dphi)^2 + dl.by.dpsi0 * d2psi0.by.dphi.2
  
  d2l.by.dbeta.2 = t(dphi.by.dbeta * d2l.by.dphi.2 * weights) %*% dphi.by.dbeta
  
  
  
  hessian = -rbind(cbind(d2l.by.dalpha.2, d2l.by.dalpha.dbeta), cbind(d2l.by.dbeta.dalpha, 
                                                                      d2l.by.dbeta.2))
  ### NB Note the extra minus sign here
  
  return(list(hessian = hessian, p0 = p0, p1 = p1, pA = pA, dpsi0.by.dtheta = dpsi0.by.dtheta, 
              dpsi0.by.dphi = dpsi0.by.dphi, dtheta.by.dalpha = dtheta.by.dalpha, 
              dphi.by.dbeta = dphi.by.dbeta, dl.by.dpsi0 = dl.by.dpsi0,
              hess_alpha = -d2l.by.dalpha.2, hess_beta = -d2l.by.dbeta.2))
  
} 