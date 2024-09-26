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
generate_data <- function(pa, pb, n, n_test, alpha, beta, gamma) {
    # simulate training data
    v.train <-  matrix(stats::runif(n * pa, min = -1, max = 1), nrow = n, ncol = pa)
    # mvnfast::rmvn(n, mu = rep(0, pa), sigma = diag(rep(1, pa)))
    # center the training data v.train
    # v.train.mean <- colMeans(v.train)
    # v.train <- scale(v.train)
    
    # generate propensity scores for training data and create binary exposure x
    gamma_intercept <- offset.compute(10000, gamma, diag(rep(1, pa)), 0.5)
    gamma_true <- c(gamma_intercept, gamma)
    pscore.true <- sigmoid(cbind(rep(1, n), v.train) %*% gamma_true)
    x.train <- stats::rbinom(n, 1, pscore.true)

    # generate binary outcome y for training data
    p0p1.true <- brm::getProbRR(v.train %*% alpha, v.train %*% beta)
    pA.true <- p0p1.true[, 1] # P(Y=1|X=0)
    pA.true[x.train == 1] <- p0p1.true[x.train == 1, 2] # replace with # P(Y=1|X=1)
    y.train <- stats::rbinom(n, 1, pA.true) # P(Y=1|X)

    ## generate test data
    # set.seed(2)
    v.test <-  matrix(stats::runif(n * pa, min = -1, max = 1), nrow = n, ncol = pa)
    # mvnfast::rmvn(n_test, mu = rep(0, pa), sigma = diag(rep(1, pa)))

    # Center the test data using the mean of the training data
    # v.test <- scale(v.test, center = v.train.mean, scale = FALSE)

    # generate propensity scores for test data and create binary exposure x
    pscore.true.test <- sigmoid(cbind(rep(1, n_test), v.test) %*% gamma_true)
    x.test <- stats::rbinom(n_test, 1, pscore.true.test)

    # generate binary outcome y for test data
    p0p1.true.test <- brm::getProbRR(v.test %*% alpha, v.test %*% beta)
    pA.true.test <- p0p1.true.test[, 1] # P(Y=1|X=0)
    pA.true.test[x.test == 1] <- p0p1.true.test[x.test == 1, 2] # replace with # P(Y=1|X=1)
    y.test <- stats::rbinom(n_test, 1, pA.true.test) # P(Y=1|X)

    return(list(v.train, x.train, y.train, v.test, x.test, y.test))
}
