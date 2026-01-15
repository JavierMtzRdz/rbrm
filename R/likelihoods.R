#' Negative Log-Likelihood
#'
#' Computes the negative log-likelihood for the RBRM model given the coefficients and data.
#'
#' @param alpha Vector of alpha coefficients (Relative Risk model).
#' @param beta Vector of beta coefficients (Baseline Risk model).
#' @param va Matrix of predictors for alpha.
#' @param vb Matrix of predictors for beta.
#' @param x Treatment vector.
#' @param y Outcome vector.
#' @param prob_fun Function to calculate probabilities (default: `getProbRR.org`).
#' @param weights Optional weights (defaults to 1).
#' @param clipping Epsilon for probability clipping.
#'
#' @return The negative log-likelihood (scalar).
#' @export
nllh <- function(alpha, beta, va, vb, x, y,
                 prob_fun = getProbRR.org,
                 weights = rep(1, length(x)),
                 clipping = 1e-10) {
    n <- length(y)
    pa <- length(alpha) # Use length of coeff vector
    pb <- length(beta)
    # Safe matrix multiplication: result is 0 vector if no columns/coefficients
    if (pa > 0 && ncol(va) == pa) logrr <- va %*% alpha else logrr <- matrix(0, nrow = n, ncol = 1)
    if (pb > 0 && ncol(vb) == pb) logop <- vb %*% beta else logop <- matrix(0, nrow = n, ncol = 1)

    if ("clipping" %in% names(formals(prob_fun))) {
        ps <- prob_fun(as.vector(logrr), as.vector(logop), clipping = clipping)
    } else {
        ps <- prob_fun(as.vector(logrr), as.vector(logop))
    }
    p0 <- ps$p0
    p1 <- ps$p1

    idx0 <- which(x == 0)
    idx1 <- which(x == 1)
    nll <- 0
    # Calculate likelihood safely, avoiding issues if idx0 or idx1 are empty
    if (length(idx0) > 0) {
        nll <- nll - sum(y[idx0] * log(p0[idx0]) + (1 - y[idx0]) * log(1 - p0[idx0]))
    }
    if (length(idx1) > 0) {
        nll <- nll - sum(y[idx1] * log(p1[idx1]) + (1 - y[idx1]) * log(1 - p1[idx1]))
    }

    # size adjustment
    nll <- nll / n

    return(nll)
}


#' Penalized Negative Log-Likelihood
#'
#' Computes the penalized negative log-likelihood for alpha and beta, adding L1 penalties.
#'
#' @param alpha Vector of alpha coefficients.
#' @param beta Vector of beta coefficients.
#' @param va Covariate matrix for alpha.
#' @param vb Covariate matrix for beta.
#' @param x Treatment vector.
#' @param y Outcome vector.
#' @param lambda Regularization parameter for alpha.
#' @param lambda_beta Regularization parameter for beta (defaults to lambda * lambda_b_prop).
#' @param lambda_b_prop Ratio of beta penalty to alpha penalty.
#' @param intercept Logical. If TRUE, the first coefficient is unpenalized.
#' @param prob_fun Probability function.
#' @param nllh_fun Underlying NLL function (defaults to `nllh`).
#' @param clipping Clipping epsilon.
#'
#' @return The penalized negative log-likelihood (scalar).
#' @export
penalized_nllh <- function(alpha, beta, va, vb, x, y,
                           lambda,
                           lambda_beta = NULL,
                           lambda_b_prop = 1,
                           intercept = F,
                           prob_fun = getProbRR.org,
                           nllh_fun = nllh,
                           clipping = 1e-10) {
    if (is.null(lambda_beta)) lambda_beta <- lambda * lambda_b_prop

    unpenalized.nllh <- nllh_fun(alpha, beta, va, vb, x, y,
        prob_fun = prob_fun,
        clipping = clipping
    )

    # Calculate L1 penalty
    # Calculate L1 penalty
    val_alpha <- if (intercept) alpha[-1] else alpha
    val_beta <- if (intercept) beta[-1] else beta
    l1_norm_alpha <- sum(abs(val_alpha))
    l1_norm_beta <- sum(abs(val_beta))

    penalty <- lambda * l1_norm_alpha + lambda_beta * l1_norm_beta

    return(unpenalized.nllh + penalty)
}
