#' Sigmoid function
#' @keywords internal
sigmoid <- function(t) {
  return(exp(t) / (1 + exp(t)))
}

#' Compute Offset for Desired Proportion
#'
#' @importFrom mvnfast rmvn
#' @keywords internal
offset.compute <- function(M, gamma, Sigma, proportion) {
  p <- length(gamma)
  x.data <- rmvn(M, mu = rep(0, p), sigma = Sigma)
  coef.fit <- x.data %*% gamma

  proportion.difference <- function(offset, coef.fit, proportion) {
    prob.test <- sigmoid(coef.fit + offset)
    return(abs(mean(round(prob.test, 0)) - proportion))
  }

  stats::optimize(
    f = proportion.difference, interval = c(-20, 20),
    coef.fit = coef.fit, proportion = proportion
  )$minimum
}

#' Soft Thresholding Operator
#' @keywords internal
soft_thres <- function(x, lambda) {
  sign(x) * pmax(0, abs(x) - lambda)
}

#' Check numerical equality
#' @keywords internal
same <- function(x, y, tolerance = .Machine$double.eps^0.5) {
  abs(x - y) < tolerance
}

#' LogSumExp (Stable)
#' @keywords internal
logsumexp <- function(x) {
  if (!is.numeric(x)) stop("x should be numeric")
  x_max <- max(x)
  x_max + log(sum(exp(x - x_max)))
}

#' Log1p (Stable)
#' @keywords internal
Log1p <- function(x) {
  if (!is.numeric(x)) stop("x should be numeric")
  y <- 1 + x
  z <- y - 1
  idx <- (z == 0)
  out <- rep(0, length(x))
  out[idx] <- x[idx]
  out[!idx] <- x[!idx] * log(y[!idx]) / z[!idx]
  return(out)
}

#' Interpolate Value based on Key mapping
#' @keywords internal
map_with_interpolation <- function(value) {
  mapping <- c(`5` = 0.26, `50` = 0.15, `150` = 0.07, `500` = .05)
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

  lower <- max(which(keys < value))
  upper <- min(which(keys > value))

  values[lower] + (values[upper] - values[lower]) * (value - keys[lower]) / (keys[upper] - keys[lower])
}

#' Find Intercept for Sigmoid Probability Target
#' @keywords internal
find_int_sigmoid <- function(lp, target) {
  f <- function(b) mean(sigmoid(b + lp)) - target
  tryCatch(uniroot(f, c(-100, 100))$root, error = function(e) 0)
}

#' Compute Intercept for General Link
#' @keywords internal
compute_intercept_sim <- function(linear_pred_main, target_prob, link_fun, linear_pred_other = NULL) {
  f <- function(b) {
    probs <- link_fun(b, linear_pred_main, linear_pred_other)
    mean(probs) - target_prob
  }
  tryCatch(uniroot(f, c(-20, 20))$root, error = function(e) -2.3)
}
