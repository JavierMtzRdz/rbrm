#' @description Sigmoid function (expit function)
#' @export
sigmoid <- function(t) {
    return(exp(t) / (1 + exp(t)))
}
#' @description function to return the offset of the sigmoid function
#' @importFrom mvnfast rmvn
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
    optimal.offset <- optimize(
        f = proportion.difference, interval = c(-20, 20),
        coef.fit = coef.fit, proportion = proportion
    )$minimum
    return(optimal.offset)
}


proximal_operator <- function(x, lambda) {
    sign(x) * pmax(0, abs(x) - lambda)
}

## Function for checking if two things are equal within numerical precision
same <- function(x, y, tolerance = .Machine$double.eps^0.5) {
    abs(x - y) < tolerance
}
