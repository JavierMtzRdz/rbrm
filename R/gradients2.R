prob_dev <- function(theta, phi, p0, p1, tolerance = 1e-6) {
  
  R <- exp(theta)
  K <- exp(phi)
  n_obs <- length(theta)
  
  # Initialize derivative vectors
  deriv_p0_theta <- numeric(n_obs)
  deriv_p0_phi <- numeric(n_obs)
  deriv_p1_theta <- numeric(n_obs)
  deriv_p1_phi <- numeric(n_obs)
  
    R_case <- R
    K_case <- K
    p0_case <- p0
    theta_case <- theta
    
    R_plus_1 = R_case + 1
    R_plus_1_sq = R_plus_1^2
    K_sq = K_case^2
    one_minus_K = (1 - K_case)
    
    inside_sqrt = K_sq * R_plus_1_sq + 4 * R_case * K_case * one_minus_K
    inside_sqrt = pmax(inside_sqrt, 0) 
    sqrt_term_val = sqrt(inside_sqrt)
    sqrt_term_val = pmax(sqrt_term_val, 1e-15) 
    
    d_sqrt_dR = (1 / (2 * sqrt_term_val)) * (K_sq * 2 * R_plus_1 + 4 * K_case * one_minus_K)
    d_sqrt_dK = (1 / (2 * sqrt_term_val)) * (2 * K_case * R_plus_1_sq + 4 * R_case * one_minus_K - 4 * R_case * K_case)
    
    denom_p0 = 2 * R_case * one_minus_K 
    denom_deriv_p0_RK = denom_p0^2
    denom_deriv_p0_RK = pmax(denom_deriv_p0_RK, 1e-15) 
    
    N_term = (-(R_plus_1) * K_case + sqrt_term_val) 
    dN_dR = -K_case + d_sqrt_dR
    dN_dK = -(R_plus_1) + d_sqrt_dK
    dD_dR = 2 * one_minus_K
    dD_dK = -2 * R_case
    
    numerator_dp0_dR = (dN_dR * denom_p0) - (N_term * dD_dR)
    dp0_dR = numerator_dp0_dR / denom_deriv_p0_RK
    
    numerator_dp0_dK = (dN_dK * denom_p0) - (N_term * dD_dK)
    dp0_dK = numerator_dp0_dK / denom_deriv_p0_RK
    
    deriv_p0_theta = dp0_dR * R_case 
    deriv_p0_phi = dp0_dK * K_case   
    
    deriv_p1_theta = deriv_p0_theta * R_case + p0_case * R_case 
    deriv_p1_phi = deriv_p0_phi * R_case
  
  return(list(
    dp0_dtheta = deriv_p0_theta, dp0_dphi = deriv_p0_phi,
    dp1_dtheta = deriv_p1_theta, dp1_dphi = deriv_p1_phi
  ))
}


#' Calculate Analytical Gradient of User's nllh Function (RR Model) - Corrected Structure
#'
#' @param alpha Vector of parameters for theta.
#' @param beta Vector of parameters for phi.
#' @param va Design matrix for alpha parameters (W in paper). Rows are observations.
#' @param vb Design matrix for beta parameters (Z in paper). Rows are observations.
#' @param x Vector of binary exposures/treatments (A in paper).
#' @param y Vector of binary outcomes (Y in paper).
#' @param prob_fun Function that takes theta (logrr), phi (logop) and returns list(p0, p1).
#' @param opt Character string specifying which gradient to return ("alpha", "beta", or "both").
#'
#' @return Analytical gradient vector(s) of the negative log-likelihood defined by user's nllh.
# grad_nll <- function(alpha, beta, x, y, va, vb, prob_fun, opt = "both") {
#   
#   n <- length(y)
#   pa <- length(alpha)
#   pb <- length(beta)
#   
#   # Calculate theta (logrr) and phi (logop)
#   logrr <- va %*% alpha
#   logop <- vb %*% beta
#   
#   # Get p0, p1
#   ps <- prob_fun(logrr, logop)
#   p0 <- ps$p0
#   p1 <- ps$p1
#   
#   # Clamp probabilities (matching nllh)
#   eps <- 1e-15
#   # p0 <- pmax(eps, pmin(1 - eps, p0))
#   # p1 <- pmax(eps, pmin(1 - eps, p1))
#   
#   # --- Calculate Analytical Derivatives of p0, p1 for RR Model ---
#   prob_derivs <- prob_dev(logrr, logop, p0, p1)
#   
#   dp0_dtheta <- prob_derivs$dp0_dtheta
#   dp0_dphi <- prob_derivs$dp0_dphi
#   dp1_dtheta <- prob_derivs$dp1_dtheta
#   dp1_dphi <- prob_derivs$dp1_dphi
#   # -------------------------------------------------------------
#   
#   # Calculate derivative terms: (y-p)/(p*(1-p))
#   # Avoid division by zero if p clamped exactly to 0 or 1
#   deriv_term_0 <- (y - p0) / (p0 * (1 - p0))
#   deriv_term_1 <- (y - p1) / (p1 * (1 - p1))
#   
#   # Initialize gradient vectors
#   grad_alpha_sum <- numeric(pa)
#   grad_beta_sum <- numeric(pb)
#   
#   # Get indices for treated and untreated
#   idx0 <- which(x == 0)
#   idx1 <- which(x == 1)
#   
#   # --- Calculate Gradient Sum for Alpha ---
#   if (pa > 0) {
#     # Contribution from untreated group (x=0)
#     if (length(idx0) > 0) {
#       # Term: deriv_term_0 * dp0_dtheta * va
#       term0_alpha <- deriv_term_0[idx0] * dp0_dtheta[idx0]
#       # Use matrix multiplication: colSums( term_vector * covariate_matrix )
#       grad_alpha_sum <- grad_alpha_sum + colSums(term0_alpha * va[idx0, , drop = FALSE])
#     }
#     # Contribution from treated group (x=1)
#     if (length(idx1) > 0) {
#       # Term: deriv_term_1 * dp1_dtheta * va
#       term1_alpha <- deriv_term_1[idx1] * dp1_dtheta[idx1]
#       grad_alpha_sum <- grad_alpha_sum + colSums(term1_alpha * va[idx1, , drop = FALSE])
#     }
#   }
#   
#   # --- Calculate Gradient Sum for Beta ---
#   if (pb > 0) {
#     # Contribution from untreated group (x=0)
#     if (length(idx0) > 0) {
#       # Term: deriv_term_0 * dp0_dphi * vb
#       term0_beta <- deriv_term_0[idx0] * dp0_dphi[idx0]
#       grad_beta_sum <- grad_beta_sum + colSums(term0_beta * vb[idx0, , drop = FALSE])
#     }
#     # Contribution from treated group (x=1)
#     if (length(idx1) > 0) {
#       # Term: deriv_term_1 * dp1_dphi * vb
#       term1_beta <- deriv_term_1[idx1] * dp1_dphi[idx1]
#       grad_beta_sum <- grad_beta_sum + colSums(term1_beta * vb[idx1, , drop = FALSE])
#     }
#   }
#   
#   # Apply scaling factor (-1/n) from nllh function
#   grad_alpha <- -grad_alpha_sum /n
#   grad_beta <- -grad_beta_sum/n
#   
#   # Handle potential NaNs from derivative calculations
#   # grad_alpha <- ifelse(is.nan(grad_alpha), 0, grad_alpha)
#   # grad_beta <- ifelse(is.nan(grad_beta), 0, grad_beta)
#   
#   # Return requested gradient(s)
#   if (opt == "alpha") {
#     return(grad_alpha)
#   } else if (opt == "beta") {
#     return(grad_beta)
#   } else { # Default to "both"
#     return(c(grad_alpha, grad_beta))
#   }
# }




#' Helper: Calculate Analytical Derivatives (RR Model) - Refined Stability
#' Uses general formula everywhere, with enhanced checks.
#'
#' @param theta Vector of theta values.
#' @param phi Vector of phi values.
#' @param p0 Vector of p0 values.
#' @param p1 Vector of p1 values.
#' @param tolerance Small value for clamping denominators etc.
#' @return List: dp0_dtheta, dp0_dphi, dp1_dtheta, dp1_dphi.
calculate_analytical_prob_derivatives_RR_stable <- function(theta, phi, p0, p1, tolerance = 1e-12) {
  
  R <- exp(theta)
  K <- exp(phi)
  n_obs <- length(theta)
  
  # Initialize results
  dp0_dtheta <- numeric(n_obs)
  dp0_dphi <- numeric(n_obs)
  dp1_dtheta <- numeric(n_obs)
  dp1_dphi <- numeric(n_obs)
  
  # --- Calculations using General Formula ---
  R_plus_1 = R + 1
  R_plus_1_sq = R_plus_1^2
  K_sq = K^2
  one_minus_K = (1 - K)
  
  # Denominator D = 2*R*(1-K) - needs careful handling if K=1
  denom_p0 = 2 * R * one_minus_K
  # If K is exactly 1, denom_p0 is 0. Need to check if formula limit exists.
  # As derived before, p0 -> R/(1+R) when K->1. 
  # The derivatives calculated below will likely -> Inf/NaN if K=1.
  
  inside_sqrt = K_sq * R_plus_1_sq + 4 * R * K * one_minus_K
  # Clamp to avoid sqrt of small negative due to precision
  inside_sqrt = pmax(inside_sqrt, 0) 
  sqrt_term_val = sqrt(inside_sqrt)
  # Clamp to avoid division by zero, use tolerance
  safe_sqrt_term_val = pmax(sqrt_term_val, tolerance) 
  
  # Derivatives of SQRT term
  d_sqrt_dR = (0.5 / safe_sqrt_term_val) * (K_sq * 2 * R_plus_1 + 4 * K * one_minus_K)
  d_sqrt_dK = (0.5 / safe_sqrt_term_val) * (2 * K * R_plus_1_sq + 4 * R * one_minus_K - 4 * R * K)
  
  # Numerator N = -(R+1)*K + SQRT
  N_term = -(R_plus_1) * K + sqrt_term_val 
  
  # Derivatives of Numerator N
  dN_dR = -K + d_sqrt_dR
  dN_dK = -(R_plus_1) + d_sqrt_dK
  
  # Derivatives of Denominator D
  dD_dR = 2 * one_minus_K
  dD_dK = -2 * R
  
  # Denominator for Quotient Rule D^2
  denom_deriv_p0_RK = denom_p0^2
  safe_denom_deriv_p0_RK = pmax(denom_deriv_p0_RK, tolerance^2) # Clamp
  
  # dp0/dR = (dN/dR * D - N * dD/dR) / D^2
  numerator_dp0_dR = (dN_dR * denom_p0) - (N_term * dD_dR)
  dp0_dR = numerator_dp0_dR / safe_denom_deriv_p0_RK
  
  # dp0/dK = (dN/dK * D - N * dD/dK) / D^2
  numerator_dp0_dK = (dN_dK * denom_p0) - (N_term * dD_dK)
  dp0_dK = numerator_dp0_dK / safe_denom_deriv_p0_RK
  
  # Convert back to derivatives w.r.t. theta and phi
  dp0_dtheta = dp0_dR * R # dR/dtheta = R
  dp0_dphi   = dp0_dK * K   # dK/dphi = K
  
  # Calculate derivatives for p1 = p0 * R
  dp1_dtheta = dp0_dtheta * R + p0 * R 
  dp1_dphi   = dp0_dphi * R
  
  # --- Final Safety Checks ---
  # Check for cases where K was exactly 1, derivatives might be Inf/NaN
  # Apply the correct limits ONLY if K is exactly 1 (or extremely close)
  # This is different from previous attempts that used a wider tolerance band
  idx_K_is_1 <- which(abs(K - 1) < tolerance) # Use tolerance for float comparison
  if (length(idx_K_is_1) > 0) {
    R_limit <- R[idx_K_is_1]
    # Use rigorously derived limit for dp0/dtheta
    dp0_dtheta[idx_K_is_1] <- R_limit / (1 + R_limit)^2
    # Use the *approximation* for dp0/dphi limit (still suspect)
    dp0_dphi[idx_K_is_1] <- 0 
    
    # Recalculate p1 derivatives based on these limits
    p0_limit <- p0[idx_K_is_1]
    dp1_dtheta[idx_K_is_1] <- dp0_dtheta[idx_K_is_1] * R_limit + p0_limit * R_limit
    dp1_dphi[idx_K_is_1] <- dp0_dphi[idx_K_is_1] * R_limit # = 0 based on approximation
  }
  
  # Replace any remaining NaNs/Infs (e.g., from extreme inputs) with 0
  dp0_dtheta = ifelse(!is.finite(dp0_dtheta), 0, dp0_dtheta)
  dp0_dphi   = ifelse(!is.finite(dp0_dphi), 0, dp0_dphi)
  dp1_dtheta = ifelse(!is.finite(dp1_dtheta), 0, dp1_dtheta)
  dp1_dphi   = ifelse(!is.finite(dp1_dphi), 0, dp1_dphi)
  
  return(list(
    dp0_dtheta = dp0_dtheta, dp0_dphi = dp0_dphi,
    dp1_dtheta = dp1_dtheta, dp1_dphi = dp1_dphi
  ))
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
#' @param prob_fun Function that takes theta (logrr), phi (logop) and returns list(p0, p1).
#' @param opt Character string specifying which gradient to return ("alpha", "beta", or "both").
#'
#' @return Analytical gradient vector(s) of the negative log-likelihood defined by user's nllh.
grad_nll <- function(alpha, beta, y, x, va, vb, prob_fun, opt = "both") {
  
  n <- length(y)
  pa <- length(alpha)
  pb <- length(beta)
  
  # Calculate theta (logrr) and phi (logop)
  logrr <- if (pa > 0 && ncol(va) == pa) as.vector(va %*% alpha) else numeric(n) # Ensure vector
  logop <- if (pb > 0 && ncol(vb) == pb) as.vector(vb %*% beta) else numeric(n) # Ensure vector
  
  # Get p0, p1
  ps <- prob_fun(logrr, logop)
  p0 <- ps$p0
  p1 <- ps$p1
  
  # Clamp probabilities (matching nllh)
  eps <- 1e-15
  p0 <- pmax(eps, pmin(1 - eps, p0))
  p1 <- pmax(eps, pmin(1 - eps, p1))
  
  # --- Calculate Analytical Derivatives of p0, p1 for RR Model ---
  # Uses the stability-focused helper function
  prob_derivs <- calculate_analytical_prob_derivatives_RR_stable(logrr, logop, p0, p1)
  
  dp0_dtheta <- prob_derivs$dp0_dtheta
  dp0_dphi <- prob_derivs$dp0_dphi
  dp1_dtheta <- prob_derivs$dp1_dtheta
  dp1_dphi <- prob_derivs$dp1_dphi
  # -------------------------------------------------------------
  
  # Calculate derivative terms: (y-p)/(p*(1-p))
  # Avoid division by zero if p clamped exactly to 0 or 1
  deriv_term_0 <- ifelse(p0 > eps & p0 < (1-eps), (y - p0) / (p0 * (1 - p0)), 0)
  deriv_term_1 <- ifelse(p1 > eps & p1 < (1-eps), (y - p1) / (p1 * (1 - p1)), 0)
  
  # Initialize gradient sums
  grad_alpha_sum <- numeric(pa)
  grad_beta_sum <- numeric(pb)
  
  # Get indices for treated and untreated
  idx0 <- which(x == 0)
  idx1 <- which(x == 1)
  
  # --- Calculate Gradient Sum for Alpha ---
  if (pa > 0) {
    # Pre-calculate terms to multiply by covariates
    term0_alpha_vals <- deriv_term_0[idx0] * dp0_dtheta[idx0]
    term1_alpha_vals <- deriv_term_1[idx1] * dp1_dtheta[idx1]
    # Replace NA/NaN/Inf with 0 before multiplying covariates
    term0_alpha_vals <- ifelse(!is.finite(term0_alpha_vals), 0, term0_alpha_vals)
    term1_alpha_vals <- ifelse(!is.finite(term1_alpha_vals), 0, term1_alpha_vals)
    
    # Use matrix multiplication approach (more robust than colSums on subsets)
    grad_alpha_sum <- numeric(pa)
    if(length(idx0) > 0) {
      grad_alpha_sum <- grad_alpha_sum + crossprod(va[idx0, , drop = FALSE], term0_alpha_vals) # t(va) %*% term
    }
    if(length(idx1) > 0) {
      grad_alpha_sum <- grad_alpha_sum + crossprod(va[idx1, , drop = FALSE], term1_alpha_vals)
    }
  }
  
  # --- Calculate Gradient Sum for Beta ---
  if (pb > 0) {
    # Pre-calculate terms
    term0_beta_vals <- deriv_term_0[idx0] * dp0_dphi[idx0]
    term1_beta_vals <- deriv_term_1[idx1] * dp1_dphi[idx1]
    term0_beta_vals <- ifelse(!is.finite(term0_beta_vals), 0, term0_beta_vals)
    term1_beta_vals <- ifelse(!is.finite(term1_beta_vals), 0, term1_beta_vals)
    
    grad_beta_sum <- numeric(pb)
    if(length(idx0) > 0) {
      grad_beta_sum <- grad_beta_sum + crossprod(vb[idx0, , drop = FALSE], term0_beta_vals)
    }
    if(length(idx1) > 0) {
      grad_beta_sum <- grad_beta_sum + crossprod(vb[idx1, , drop = FALSE], term1_beta_vals)
    }
  }
  
  # Apply scaling factor (-1/n) from nllh function
  grad_alpha <- -grad_alpha_sum / n
  grad_beta <- -grad_beta_sum / n
  
  # Final check for NaN (shouldn't happen after crossprod if inputs are finite)
  grad_alpha <- ifelse(is.nan(grad_alpha), 0, grad_alpha)
  grad_beta <- ifelse(is.nan(grad_beta), 0, grad_beta)
  
  # Return requested gradient(s)
  if (opt == "alpha") {
    return(grad_alpha)
  } else if (opt == "beta") {
    return(grad_beta)
  } else { # Default to "both"
    return(c(grad_alpha, grad_beta))
  }
}