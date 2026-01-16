# derivative of pi0 with respect to theta
dp0_theta <- function(theta, phi,
                      ps_spe = c("Richardson", "Pozza"), ep = 1e-8) {
  ps_spe <- rlang::arg_match(ps_spe)

  x <- phi - theta
  expm1_theta <- expm1(theta)
  expm1_phi <- expm1(phi)

  if (ps_spe == "Richardson") {
    # Find indices for each condition.
    is_boundary <- (phi < -12) | (phi > 12) | (theta < -12) | (theta > 12)
    is_south_edge <- (theta < -12) | ((phi < -12) & (theta < 0))
    is_west_edge <- (theta > 12) | ((phi < -12) & (theta > 0))
    is_phi_zero <- same(phi, 0)

    x <- phi - theta

    dp0.theta <- rep(NA, length(theta))

    # Apply the conditions to calculate p0.

    # Case on the boundary
    idx_south_not_ext <- (is_boundary & is_south_edge) & ((x < 17) & (x > (-500)))
    idx_south_ext <- (is_boundary & is_south_edge) & !((x < 17) & (x > (-500)))

    y_val <- 4 * exp(-x[idx_south_not_ext])
    sqrt_term <- sqrt(y_val + 1)
    term1 <- -y_val / (1 + sqrt_term)
    dp0.theta[idx_south_not_ext] <-
      (0.5 * term1 * sqrt_term * exp(x[idx_south_not_ext]) + 1.0) / sqrt_term

    dp0.theta[idx_south_ext] <- 0
    dp0.theta[is_boundary & is_west_edge] <- 0
    dp0.theta[is_boundary & !(is_south_edge | is_west_edge)] <- ifelse(-exp(-theta[is_boundary & !(is_south_edge | is_west_edge)]) < -1, 0,
      -exp(-theta[is_boundary & !(is_south_edge | is_west_edge)])
    )


    # Case not on the boundary
    not_boundary_indices <- !is_boundary

    dp0.theta[not_boundary_indices & is_phi_zero] <- -exp(theta[not_boundary_indices & is_phi_zero]) / (1 + exp(theta[not_boundary_indices & is_phi_zero]))^2

    # Calculate for the quadratic equation case
    quadratic_indices <- not_boundary_indices & !is_phi_zero

    dp0.theta[quadratic_indices] <- (-exp(phi - theta) / (2 * expm1_phi) + exp(phi) /
      (expm1_phi * sqrt(4 * exp(phi + theta) + expm1_theta^2 * exp(2 * phi))) + (exp(-theta) * (-expm1_theta) *
      exp(2 * phi)) / (2 * expm1_phi * sqrt(4 * exp(phi +
      theta) + expm1_theta^2 * exp(2 * phi))))[quadratic_indices]

    return(dp0.theta)
  }

  if (ps_spe == "Pozza") {
    dp0.theta <- ((exp(-phi - theta) * (exp(phi + theta) -
      (2 * ((exp(theta) + 1) * exp(phi) + 1) * exp(phi + theta) -
        4 * exp(2 * phi + theta)) /
        (2 * sqrt(((exp(theta) + 1) * exp(phi) + 1)^2 -
          4 * exp(2 * phi + theta))))) / 2 - (exp(-phi - theta) *
      (-sqrt(((exp(theta) + 1) * exp(phi) + 1)^2 -
        4 * exp(2 * phi + theta)) + (exp(theta) + 1) * exp(phi) + 1)) / 2)
    return(dp0.theta)
  }
}
# derivative of pi0 with respect to phi

dp0_phi <- function(theta, phi,
                    ps_spe = c("Richardson", "Pozza"), ep = 1e-8) {
  ps_spe <- rlang::arg_match(ps_spe)

  x <- phi - theta
  expm1_theta <- expm1(theta)
  expm1_phi <- expm1(phi)

  if (ps_spe == "Richardson") {
    # Find indices for each condition.
    is_boundary <- (phi < -12) | (phi > 12) | (theta < -12) | (theta > 12)
    is_south_edge <- (theta < -12) | ((phi < -12) & (theta < 0))
    is_west_edge <- (theta > 12) | ((phi < -12) & (theta > 0))
    is_phi_zero <- same(phi, 0)


    dp0.phi <- rep(NA, length(phi))

    # Apply the conditions to calculate p0.

    # Case on the boundary
    idx_south_not_ext <- (is_boundary & is_south_edge) & ((x < 17) & (x > (-500)))
    idx_south_ext <- (is_boundary & is_south_edge) & !((x < 17) & (x > (-500)))

    y_val <- 4 * exp(-x[idx_south_not_ext])
    sqrt_term <- sqrt(y_val + 1)
    term1 <- -y_val / (1 + sqrt_term) # (1 - sqrt(1+y))
    dp0.phi[idx_south_not_ext] <- (0.5 * exp(x[idx_south_not_ext]) * term1 + 1) / sqrt_term

    dp0.phi[idx_south_ext] <- 0
    dp0.phi[is_boundary & is_west_edge] <- 0
    dp0.phi[is_boundary & !(is_south_edge | is_west_edge)] <- 0

    # Case not on the boundary
    not_boundary_indices <- !is_boundary

    dp0.phi[not_boundary_indices & is_phi_zero] <- (exp(theta[not_boundary_indices & is_phi_zero]) /
      (exp(theta[not_boundary_indices & is_phi_zero]) + 1)^3)


    # Calculate for the quadratic equation case
    quadratic_indices <- not_boundary_indices & !is_phi_zero

    dp0.phi[quadratic_indices] <- (-((exp(theta) + 1) * exp(phi)) / (2 * exp(theta) *
      expm1_phi^2) + exp(phi) / (expm1_phi^2 *
      sqrt(4 * exp(phi + theta) + expm1_theta^2 * exp(2 * phi))
    ) + (exp(-theta) * (exp(2 * theta) + 1) * exp(2 * phi)) /
      (2 * expm1_phi^2 * sqrt(4 * exp(phi + theta) + expm1_theta^2 * exp(2 * phi))))[quadratic_indices]

    return(dp0.phi)
  }

  if (ps_spe == "Pozza") {
    dp0.phi <- ((exp(-phi - theta) * ((exp(theta) + 1) * exp(phi) - (2 * (exp(theta) + 1) * exp(phi) * ((exp(theta) + 1) * exp(phi) + 1) - 8 * exp(2 * phi + theta)) /
      (2 * sqrt(((exp(theta) + 1) * exp(phi) + 1)^2 - 4 * exp(2 * phi + theta))))) / 2 -
      (exp(-phi - theta) * (-sqrt(((exp(theta) + 1) * exp(phi) + 1)^2 - 4 * exp(2 * phi + theta)) + (exp(theta) + 1) * exp(phi) + 1)) / 2)
    return(dp0.phi)
  }
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
#' @param clipping Small epsilon for probability clipping (defaults to 1e-10).
#'
#' @return Analytical gradient vector(s) of the negative log-likelihood defined by user's nllh.
#' @export
grad_nll <- function(alpha, beta, y, x, va, vb,
                     prob_fun, opt = c("both", "alpha", "beta"),
                     method = c("analytical", "numerical"),
                     clipping = 1e-10) {
  # Argument matching
  method <- rlang::arg_match(method)
  opt <- rlang::arg_match(opt)

  n <- length(y)
  pa <- length(alpha)
  pb <- length(beta)

  # Calculate theta, phi
  theta <- as.vector(va %*% alpha)
  phi <- as.vector(vb %*% beta)

  # Calculate probabilities with explicit clipping
  # We check if prob_fun accepts 'clipping' argument to be safe, though internal ones do.
  if ("clipping" %in% names(formals(prob_fun))) {
    ps <- prob_fun(theta, phi, clipping = clipping)
  } else {
    ps <- prob_fun(theta, phi)
  }

  p0 <- ps$p0
  p1 <- ps$p1

  # Derivatives of log-likelihood_i w.r.t p1_i and p0_i
  # Note: if p1 is clipped, this value might be large, but will be zeroed out later.
  # We use the clipped probabilities for this calculation to match the NLL evaluation point.
  dllh_dp1 <- (y * x) / p1 - ((1 - y) * x) / (1 - p1)
  dllh_dp0 <- (y * (1 - x)) / p0 - ((1 - y) * (1 - x)) / (1 - p0)

  # Identify observations where the relevant probability was clipped
  # pi = p1 if x=1, p0 if x=0
  pi_val <- p1 * x + p0 * (1 - x)
  is_clipped <- (pi_val <= clipping) | (pi_val >= (1 - clipping))

  if (opt != "beta") {
    if (method == "analytical") {
      dp0_dtheta <- dp0_theta(theta, phi,
        ps_spe = class(ps)
      )
    }
    if (method == "numerical") {
      dp0_dtheta <- numDeriv::grad(
        func = function(theta_val) prob_fun(theta_val, phi, clipping = clipping)$p0,
        x = theta,
        method = "simple"
      )
    }

    dp1_dtheta <- (p0 + dp0_dtheta) * exp(theta)
    # Check for boundary consistency explicitly for p1 derivative if needed,
    # but the generic zeroing below covers the NLL gradient.
    dp1_dtheta[which(same(p0 + dp0_dtheta, 0))] <- 0

    grad_alpha_sum <- numeric(pa)
    inner_alpha <- (dllh_dp1 * dp1_dtheta + dllh_dp0 * dp0_dtheta)

    # Zero out gradient for clipped observations
    inner_alpha[is_clipped] <- 0

    grad_alpha <- -(t(va) %*% inner_alpha) / n
  }
  if (opt != "alpha") {
    if (method == "analytical") {
      dp0_dphi <- dp0_phi(theta, phi,
        ps_spe = class(ps)
      )
    }
    if (method == "numerical") {
      dp0_dphi <- numDeriv::grad(
        func = (function(phi_val) prob_fun(theta, phi_val, clipping = clipping)$p0),
        x = phi,
        method = "simple"
      )
    }

    dp1_dphi <- dp0_dphi * exp(theta)
    dp1_dphi[which(same(dp0_dphi, 0))] <- 0

    grad_beta_sum <- numeric(pb)
    inner_beta <- (dllh_dp1 * dp1_dphi + dllh_dp0 * dp0_dphi)

    # Zero out gradient for clipped observations
    inner_beta[is_clipped] <- 0

    grad_beta <- -(t(vb) %*% inner_beta) / n
  }
  if (opt == "alpha") {
    return(grad_alpha)
  }
  if (opt == "beta") {
    return(grad_beta)
  }
  if (opt == "both") {
    return(list(
      grad_alpha = grad_alpha,
      grad_beta = grad_beta
    ))
  }
}

#' @export
grad_nll_k <- function(alpha, beta, y, x, va, vb,
                       prob_fun,
                       opt = c("alpha", "beta"),
                       k_index,
                       method = c("numerical", "analytical")) {
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


  # Initialize the required gradient component
  grad_kth_component <- NA

  if (opt == "alpha") {
    dp0_dtheta_val <- NULL
    if (method == "analytical") {
      dp0_dtheta_val <- dp0_theta(theta, phi, ps_spe = class(ps))
    } else if (method == "numerical") {
      get_p0_for_theta <- function(theta_val) {
        return(prob_fun(theta_val, phi)$p0)
      }
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
      get_p0_for_phi <- function(phi_val) {
        return(prob_fun(theta, phi_val)$p0)
      }
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
hessian_or <- function(y, x, va, vb, alpha.ml, beta.ml, weights) {
  # calculating the Hessian using the second derivative have to do so
  # because under mis-specification of models Hessian no longer equals the
  # square of the first order derivatives

  p0p1 <- brm::getProbRR(va %*% alpha.ml, vb %*% beta.ml)
  # p0p1 = cbind(p0, p1): n * 2 matrix
  p0 <- p0p1[, 1]
  p1 <- p0p1[, 2]
  n <- nrow(va)
  pA <- p0
  pA[x == 1] <- p1[x == 1]


  ### Building blocks

  dpsi0.by.dtheta <- -(1 - p0) / (1 - p0 + 1 - p1)
  dpsi0.by.dphi <- (1 - p0) * (1 - p1) / (1 - p0 + 1 - p1)

  dtheta.by.dalpha <- va
  dphi.by.dbeta <- vb

  dl.by.dpsi0 <- (y - pA) / (1 - pA)
  d2l.by.dpsi0.2 <- (y - 1) * pA / ((1 - pA)^2)


  ###### d2l.by.dalpha.2

  d2psi0.by.dtheta.2 <- ((p0 - p1) * dpsi0.by.dtheta - (1 - p0) * p1) / ((1 -
    p0 + 1 - p1)^2)

  d2l.by.dtheta.2 <- d2l.by.dpsi0.2 * (dpsi0.by.dtheta + x)^2 + dl.by.dpsi0 *
    d2psi0.by.dtheta.2

  d2l.by.dalpha.2 <- t(dtheta.by.dalpha * d2l.by.dtheta.2 * weights) %*%
    dtheta.by.dalpha


  ###### d2l.by.dalpha.dbeta

  d2psi0.by.dtheta.dphi <- (1 - p0) * (1 - p1) * (p0 - p1) / (1 - p0 + 1 -
    p1)^3

  d2l.by.dtheta.dphi <- d2l.by.dpsi0.2 * (dpsi0.by.dtheta + x) * dpsi0.by.dphi +
    dl.by.dpsi0 * d2psi0.by.dtheta.dphi

  d2l.by.dalpha.dbeta <- t(dtheta.by.dalpha * d2l.by.dtheta.dphi * weights) %*%
    dphi.by.dbeta
  d2l.by.dbeta.dalpha <- t(d2l.by.dalpha.dbeta)
  # d2l.by.dalpha.dbeta is symmetric itself if (because) va=vb


  #### d2l.by.dbeta2

  d2psi0.by.dphi.2 <- (-(p0 * (1 - p1)^2 + p1 * (1 - p0)^2) / (1 - p0 + 1 -
    p1)^2) * dpsi0.by.dphi

  d2l.by.dphi.2 <- d2l.by.dpsi0.2 * (dpsi0.by.dphi)^2 + dl.by.dpsi0 * d2psi0.by.dphi.2

  d2l.by.dbeta.2 <- t(dphi.by.dbeta * d2l.by.dphi.2 * weights) %*% dphi.by.dbeta


  hessian <- -rbind(cbind(d2l.by.dalpha.2, d2l.by.dalpha.dbeta), cbind(
    d2l.by.dbeta.dalpha,
    d2l.by.dbeta.2
  ))
  ### NB Note the extra minus sign here

  return(list(
    hessian = hessian, p0 = p0, p1 = p1, pA = pA, dpsi0.by.dtheta = dpsi0.by.dtheta,
    dpsi0.by.dphi = dpsi0.by.dphi, dtheta.by.dalpha = dtheta.by.dalpha,
    dphi.by.dbeta = dphi.by.dbeta, dl.by.dpsi0 = dl.by.dpsi0,
    hess_alpha = -d2l.by.dalpha.2, hess_beta = -d2l.by.dbeta.2
  ))
}
