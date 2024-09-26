#' Sigmoid (expit) function
#'
#' Computes the sigmoid (or expit) function, which maps real-valued inputs to the (0, 1) range.
#'
#' @param t A numeric value or vector of values.
#' 
#' @return The sigmoid of the input value(s).
#' 
#' @examples
#' sigmoid(0) # returns 0.5
#' sigmoid(c(-1, 0, 1)) # returns c(0.2689414, 0.5, 0.7310586)
#' 
#' @export
sigmoid <- function(t) {
  return(exp(t) / (1 + exp(t)))
}

#' Compute Offset for a Desired Proportion Using the Sigmoid Function
#'
#' This function computes an offset to apply to a linear combination of predictors such that 
#' the expected proportion of positive outcomes matches the specified proportion.
#'
#' @param M Integer. The number of random samples to generate.
#' @param gamma A numeric vector of coefficients.
#' @param Sigma A covariance matrix for the random samples.
#' @param proportion The target proportion of positive outcomes (between 0 and 1).
#'
#' @return The computed offset that achieves the desired proportion.
#'
#' @importFrom mvnfast rmvn
#'
#' @examples
#' M <- 1000
#' gamma <- c(0.5, -0.2, 0.1)
#' Sigma <- diag(3)
#' proportion <- 0.7
#' offset <- offset.compute(M, gamma, Sigma, proportion)
#'
#' @export
offset.compute <- function(M, gamma, Sigma, proportion) {
  # Number of coefficients
  p <- length(gamma)
  
  # Computing proportion for given offset
  x.data <- rmvn(M, mu = rep(0, p), sigma = Sigma)
  coef.fit <- x.data %*% gamma
  
  proportion.difference <- function(offset, coef.fit, proportion) {
    prob.test <- sigmoid(coef.fit + offset)
    computed.proportion <- mean(round(prob.test, 0))
    return(abs(computed.proportion - proportion))
  }
  
  # Offset computation
  optimal.offset <- stats::optimize(
    f = proportion.difference, interval = c(-20, 20),
    coef.fit = coef.fit, proportion = proportion
  )$minimum
  return(optimal.offset)
}

#' Proximal Operator Function
#'
#' Computes the proximal operator for L1-regularization (soft thresholding).
#'
#' @param x A numeric value or vector of values.
#' @param lambda A regularization parameter. Must be a positive value.
#'
#' @return The result of applying the proximal operator to the input value(s).
#'
#' @examples
#' proximal_operator(c(1, -2, 3), 1) # returns c(0, -1, 2)
#' 
#' @export
proximal_operator <- function(x, lambda) {
  sign(x) * pmax(0, abs(x) - lambda)
}

#' Check if Two Values are Numerically Equal
#'
#' Checks if two values (or vectors) are equal within numerical precision.
#'
#' @param x A numeric value or vector of values.
#' @param y A numeric value or vector of values.
#' @param tolerance A numeric value indicating the allowable tolerance for equality. Default is the square root of machine precision.
#'
#' @return TRUE if the values are within the specified tolerance, otherwise FALSE.
#'
#' @examples
#' same(0.1 + 0.2, 0.3) # returns TRUE
#' same(c(0.1 + 0.2, 0.3), c(0.3, 0.3)) # returns TRUE
#'
#' @export
same <- function(x, y, tolerance = .Machine$double.eps^0.5) {
  abs(x - y) < tolerance
}
