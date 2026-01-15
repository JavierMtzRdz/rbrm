#' Generative model for rbrm
#'
#' This function generates simulated data for testing and studying the RBRM model.
#' It allows control over sample size, number of covariates, true coefficients,
#' and treatment prevalence.
#'
#' @param pa Number of alpha coefficients (predictors for Relative Risk) without the intercept.
#' @param pb Number of beta coefficients (predictors for Baseline Risk) without the intercept.
#' @param n Sample size for the training set.
#' @param n_test Sample size for the test set.
#' @param alpha True alpha coefficients (including intercept).
#' @param beta True beta coefficients (including intercept).
#' @param gamma True gamma coefficients for propensity score model (EXCLUDING intercept).
#'              The intercept is estimated automatically to match `treatment_prob`.
#' @param treatment_prob Target probability of treatment (default 0.5).
#'
#' @return A list containing training and test data: v (covariates), x (treatment), y (outcome),
#'         and true probabilities.
#' @export
generate_data <- function(pa, pb, n, n_test, alpha = NULL, beta = NULL, gamma = NULL, treatment_prob = 0.5,
                          target_beta_prob = 0.1) {
  if (any(is.null(alpha), is.null(beta), is.null(gamma))) {
    gen_t_vals <- true_vals(pa, pb)
  }
  if (is.null(alpha)) alpha <- gen_t_vals$true_alphas
  if (is.null(beta)) beta <- gen_t_vals$true_betas
  if (is.null(gamma)) gamma <- gen_t_vals$true_gammas

  if (length(alpha) != pa + 1) cli::cli_abort("Length of `alpha` ({length(alpha)}) must match `pa + 1` ({pa + 1}).")
  if (length(beta) != pb + 1) cli::cli_abort("Length of `beta` ({length(beta)}) must match `pb + 1` ({pb + 1}).")
  if (length(gamma) != pa) cli::cli_warn("Length of `gamma` ({length(gamma)}) typically matches `pa` ({pa}) for simulation, using provided length.")

  # Simulate Covariates (Uniform -1 to 1)
  v.train <- matrix(stats::runif(n * pa, min = -1, max = 1), nrow = n, ncol = pa)

  # Propensity Score Model
  # Compute intercept for propensity score to match target prevalence

  # Helper to find intercept for Sigmoid link
  find_int_sigmoid <- function(lp, target) {
    f <- function(b) mean(sigmoid(b + lp)) - target
    tryCatch(uniroot(f, c(-100, 100))$root, error = function(e) 0)
  }

  gamma_linear_pred <- v.train %*% gamma
  gamma_intercept <- find_int_sigmoid(gamma_linear_pred, treatment_prob)
  gamma_true <- c(gamma_intercept, gamma)

  # Calculate PS and Treatment
  pscore.true <- sigmoid(cbind(1, v.train) %*% gamma_true)
  x.train <- stats::rbinom(n, 1, pscore.true)

  # Outcome Model
  alpha_true <- alpha
  beta_true <- beta

  # Linear Predictors
  theta_train <- cbind(1, v.train) %*% alpha_true
  phi_train <- cbind(1, v.train) %*% beta_true

  # True Probabilities using brm or internal implementation
  p0p1.true <- getProbRR.org(theta_train, phi_train)
  pA.true <- p0p1.true$p0 # P(Y=1|X=0)
  pA.true[x.train == 1] <- p0p1.true$p1[x.train == 1] # P(Y=1|X=1)

  y.train <- stats::rbinom(n, 1, pA.true) # Observed outcome

  # Simulate Test Data
  v.test <- matrix(stats::runif(n_test * pa, min = -1, max = 1), nrow = n_test, ncol = pa)

  pscore.true.test <- sigmoid(cbind(1, v.test) %*% gamma_true)
  x.test <- stats::rbinom(n_test, 1, pscore.true.test)

  theta_test <- cbind(1, v.test) %*% alpha_true
  phi_test <- cbind(1, v.test) %*% beta_true

  p0p1.true.test <- getProbRR.org(theta_test, phi_test)
  pA.true.test <- p0p1.true.test$p0
  pA.true.test[x.test == 1] <- p0p1.true.test$p1[x.test == 1]

  y.test <- stats::rbinom(n_test, 1, pA.true.test)

  # Return structured list
  return(list(
    # Training
    v.train = v.train,
    x.train = x.train,
    y.train = y.train,
    p0.train = p0p1.true$p0,
    p1.train = p0p1.true$p1,
    pscore.train = pscore.true,

    # Testing
    v.test = v.test,
    x.test = x.test,
    y.test = y.test,
    p0.test = p0p1.true.test$p0,
    p1.test = p0p1.true.test$p1,
    pscore.test = pscore.true.test,

    # Truth
    true_alpha = alpha_true,
    true_beta = beta_true,
    true_gamma = gamma_true
  ))
}

#' @export
map_with_interpolation <- function(value) {
  mapping <- c(
    `5` = 0.26,
    `50` = 0.15,
    `150` = 0.07,
    `500` = .05
  )
  mapping <- mapping[order(as.numeric(names(mapping)))]
  keys <- as.numeric(names(mapping))
  values <- as.numeric(mapping)

  if (value %in% keys) {
    return(mapping[as.character(value)])
  }
  if (value < min(keys)) {
    return(values[1])
  }
  if (value > max(keys)) {
    return(values[length(values)])
  }

  lower_index <- max(which(keys < value))
  upper_index <- min(which(keys > value))
  x0 <- keys[lower_index]
  x1 <- keys[upper_index]
  y0 <- values[lower_index]
  y1 <- values[upper_index]

  return(y0 + (y1 - y0) * (value - x0) / (x1 - x0))
}

#' @export
true_vals <- function(pa, pb = pa) {
  # Helper to find intercept for generic link function via simulation
  # link_fun(intercept, linear_pred) -> probabilities
  compute_intercept_sim <- function(coefs, target_prob, link_fun, linear_pred_other = NULL) {
    # Simulate covariates if not provided
    n_vars <- length(coefs)
    n_sim <- 5000
    v_sim <- matrix(stats::runif(n_sim * n_vars, min = -1, max = 1), nrow = n_sim)
    lp <- v_sim %*% coefs

    # If there's another linear predictor (alpha) needed for the link
    if (is.null(linear_pred_other) && !is.null(formals(link_fun)$lp_other)) {
      stop("Need other linear predictor for this link.")
    }

    # Objective: find b such that mean(link(b, lp, lp_other)) = target
    f <- function(b) {
      probs <- link_fun(b, lp, linear_pred_other)
      mean(probs) - target_prob
    }

    tryCatch(uniroot(f, c(-20, 20))$root, error = function(e) 0)
  }

  # Generate Alpha
  p_a <- map_with_interpolation(pa)
  n_eff_a <- round(pa * p_a)
  n_null_a <- pa - n_eff_a * 2

  alpha_slopes <- c(rep(1, n_eff_a), rep(-1, n_eff_a), rep(0, n_null_a))
  true_alphas <- c(0, alpha_slopes) # Intercept = 0

  # Generate Gamma
  gamma_slopes <- c(rep(0.1, n_eff_a), rep(-0.5, n_eff_a), rep(0, n_null_a)) * stats::rnorm(pa)
  true_gammas <- gamma_slopes

  # Generate Beta
  p_b <- map_with_interpolation(pb)
  n_eff_b <- round(pb * p_b)
  n_null_b <- pb - n_eff_b * 2

  beta_slopes <- c(rep(-0.5, n_eff_b), rep(1, n_eff_b), rep(0, n_null_b))

  # Beta Intercept (Target P0 = 0.1)
  # P0 = f(theta, phi). phi = b + beta'V. theta = 0 + alpha'V.
  # We simulate V using PA dimensions for theta.
  n_sim <- 5000
  v_sim_a <- matrix(stats::runif(n_sim * pa, min = -1, max = 1), nrow = n_sim)
  theta_sim <- v_sim_a %*% alpha_slopes # Intercept is 0

  # For beta, if pb != pa, use specific V simulation if needed
  if (pb == pa) {
    lp_beta_sim <- v_sim_a %*% beta_slopes
  } else {
    v_sim_b <- matrix(stats::runif(n_sim * pb, min = -1, max = 1), nrow = n_sim)
    lp_beta_sim <- v_sim_b %*% beta_slopes
  }

  # Define RB-Link for calculating P0
  # b is intercept for phi
  link_rb_p0 <- function(b, lp_b, theta) {
    phi <- b + lp_b
    getProbRR.org(as.vector(theta), as.vector(phi))$p0
  }

  f_beta <- function(b) mean(link_rb_p0(b, lp_beta_sim, theta_sim)) - 0.1
  beta_int <- tryCatch(uniroot(f_beta, c(-20, 20))$root, error = function(e) -2.3)

  true_betas <- c(beta_int, beta_slopes)

  return(list(
    true_alphas = true_alphas,
    true_betas = true_betas,
    true_gammas = true_gammas
  ))
}
