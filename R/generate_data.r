#' Generative model for rbrm
#' 
#' This file/function generates simulated data for studying the rbrm model
#'
#' @param pa the number of alpha coefficients without the intercept
#' @param pb the number of beta coefficients without the intercept
#' @param n number of samples in training set
#' @param n_test the number of samples in the test set
#' @param alpha the true alpha coefficients that are used to generate data.
#' This list should include the intercept coefficient.
#' @param beta the true beta coefficients that are used to generate data.
#' This list should include the intercept coefficient.
#' @param gamma the true gamma coefficients that are used to generate
#' the propensity for exposure This list should not include the intercept.
#' The intercept value will be estimated numerically to achieve an
#' even proportion of exposure in the training data.
#'
#' @return data
#' @export
#' @examples
#' #generate_data(
#' #    pa = 10, pb = 10, n = 40, n_test = 100,
#' #    alpha = c(0.5, rep(5, 5), rep(0, 5)), beta = c(0.5, rep(5, 5), rep(0, 5)),
#' #    gamma = c(rep(5, 5), rep(0, 5)))
generate_data <- function(pa, pb, n, n_test, alpha, beta, gamma, misspec_nuisance = FALSE, misspec_propensity = FALSE) {
  
  # Simulate training data
  v.train <- matrix(stats::runif(n * pa, min = -1, max = 1), nrow = n, ncol = pa)
  
  # Generate propensity scores for training data
  if (misspec_propensity) {
    # Misspecified propensity model: use only the first half of covariates
    gamma_true <- c(0, gamma[1:(pa/2)], rep(0, pa/2))
  } else {
    # Correctly specified propensity model: use all covariates
    gamma_true <- c(0, gamma)
  }
  pscore.true <- sigmoid(cbind(rep(1, n), v.train) %*% gamma_true)

  x.train <- stats::rbinom(n, 1, pscore.true)
  
  # Generate outcome probabilities for training data
  if (misspec_nuisance) {
    # Misspecified nuisance model: use only the first half of covariates
    alpha_true <- c(alpha[1:(pa/2)], rep(0, pa/2))
    beta_true <- c(beta[1:(pa/2)], rep(0, pa/2))
  } else {
    # Correctly specified nuisance model: use all covariates
    alpha_true <- alpha
    beta_true <- beta
  }
  p0p1.true <- brm::getProbRR(v.train %*% alpha_true, v.train %*% beta_true)
  pA.true <- p0p1.true[, 1] # P(Y=1|X=0)
  pA.true[x.train == 1] <- p0p1.true[x.train == 1, 2] # P(Y=1|X=1)
  y.train <- stats::rbinom(n, 1, pA.true) # P(Y=1|X)
  
  # Simulate test data
  v.test <- matrix(stats::runif(n_test * pa, min = -1, max = 1), nrow = n_test, ncol = pa)
  
  # Generate propensity scores for test data
  pscore.true.test <- sigmoid(cbind(rep(1, n_test), v.test) %*% gamma_true)
  x.test <- stats::rbinom(n_test, 1, pscore.true.test)
  
  # Generate outcome probabilities for test data
  p0p1.true.test <- brm::getProbRR(v.test %*% alpha_true, v.test %*% beta_true)
  pA.true.test <- p0p1.true.test[, 1] # P(Y=1|X=0)
  pA.true.test[x.test == 1] <- p0p1.true.test[x.test == 1, 2] # P(Y=1|X=1)
  y.test <- stats::rbinom(n_test, 1, pA.true.test) # P(Y=1|X)
  
  # Return training and test data
  return(list(v.train = v.train, x.train = x.train, y.train = y.train, 
              v.test = v.test, x.test = x.test, y.test = y.test))
}


#' @export
generate_data2 <- function(pa, pb, n, n_test, alpha, beta, gamma, treatment_prob = 0.5, 
                          misspec_nuisance = FALSE, misspec_propensity = FALSE) {

  # Function to compute intercept for desired treatment prevalence
  compute_intercept <- function(gamma, target_prevalence) {
    # Simulate a large dataset to estimate the intercept
    n_large <- 10000
    v_large <- matrix(stats::runif(n_large * length(gamma)), nrow = n_large)
    linear_predictor <- v_large %*% gamma
    intercept <- uniroot(function(b) {
      mean(sigmoid(b + linear_predictor)) - target_prevalence
    }, interval = c(-20, 20))$root
    return(intercept)
  }
  

  # Simulate training data
  v.train <- matrix(stats::runif(n * pa, min = -1, max = 1), nrow = n, ncol = pa)
  
  # Compute intercept for propensity score model to achieve desired treatment probability
  gamma_intercept <- compute_intercept(gamma, target_prevalence = treatment_prob)
  gamma_true <- c(gamma_intercept, gamma)
  
  # Generate propensity scores for training data
  if (misspec_propensity) {
    # Misspecified propensity model: use only the first half of covariates
    gamma_true <- c(gamma_intercept, gamma[1:(pa/2)], rep(0, pa/2))
  }
  pscore.true <- sigmoid(cbind(rep(1, n), v.train) %*% gamma_true)
  x.train <- stats::rbinom(n, 1, pscore.true)
  
  # Generate outcome probabilities for training data
  if (misspec_nuisance) {
    # Misspecified nuisance model: use only the first half of covariates
    alpha_true <- c(alpha[1:(pa/2)], rep(0, pa/2))
    beta_true <- c(beta[1:(pa/2)], rep(0, pa/2))
  } else {
    # Correctly specified nuisance model: use all covariates
    alpha_true <- alpha
    beta_true <- beta
  }
  p0p1.true <- brm::getProbRR(v.train %*% alpha_true, v.train %*% beta_true)
  pA.true <- p0p1.true[, 1] # P(Y=1|X=0)
  pA.true[x.train == 1] <- p0p1.true[x.train == 1, 2] # P(Y=1|X=1)
  y.train <- stats::rbinom(n, 1, pA.true) # P(Y=1|X)
  
  # Simulate test data
  v.test <- matrix(stats::runif(n_test * pa, min = -1, max = 1), nrow = n_test, ncol = pa)
  
  # Generate propensity scores for test data
  pscore.true.test <- sigmoid(cbind(rep(1, n_test), v.test) %*% gamma_true)
  x.test <- stats::rbinom(n_test, 1, pscore.true.test)
  
  # Generate outcome probabilities for test data
  p0p1.true.test <- brm::getProbRR(v.test %*% alpha_true, v.test %*% beta_true)
  pA.true.test <- p0p1.true.test[, 1] # P(Y=1|X=0)
  pA.true.test[x.test == 1] <- p0p1.true.test[x.test == 1, 2] # P(Y=1|X=1)
  y.test <- stats::rbinom(n_test, 1, pA.true.test) # P(Y=1|X)
  
  # Return training and test data, including predicted p0 and p1
  return(list(v.train = v.train, x.train = x.train, y.train = y.train, 
              p0.train = p0p1.true[, 1], p1.train = p0p1.true[, 2], 
              v.test = v.test, x.test = x.test, y.test = y.test, 
              p0.test = p0p1.true.test[, 1], p1.test = p0p1.true.test[, 2]))
}
