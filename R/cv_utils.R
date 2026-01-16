#' Standardize Design Matrices for RBRM
#'
#' @param va Matrix for alpha.
#' @param vb Matrix for beta.
#' @return List with standardized matrices and scaler info.
#' @export
standardize_data <- function(va, vb = NULL) {
    p_a <- ncol(va)
    scaler_a <- list(center = rep(0, p_a), scale = rep(1, p_a))

    # Standardize va
    if (p_a > 0) {
        # Check for intercept
        is_intercept_a <- apply(va, 2, function(x) var(x) == 0)

        # Scale non-intercept columns
        if (!all(is_intercept_a)) {
            va_sub <- va[, !is_intercept_a, drop = FALSE]
            center_a <- colMeans(va_sub, na.rm = TRUE)
            n <- nrow(va_sub)
            scale_a <- apply(va_sub, 2, sd, na.rm = TRUE) * sqrt((n - 1) / n)

            # Handle zero variance
            zero_var <- scale_a < .Machine$double.eps
            scale_a[zero_var] <- 1

            va_scaled_sub <- scale(va_sub, center = center_a, scale = scale_a)
            va[, !is_intercept_a] <- va_scaled_sub

            # Store info (mapping back to original indices)
            scaler_a$center[!is_intercept_a] <- center_a
            scaler_a$scale[!is_intercept_a] <- scale_a
        }
    }

    # Standardize vb
    scaler_b <- NULL
    if (!is.null(vb)) {
        p_b <- ncol(vb)
        scaler_b <- list(center = rep(0, p_b), scale = rep(1, p_b))

        is_intercept_b <- apply(vb, 2, function(x) var(x) == 0)

        if (!all(is_intercept_b)) {
            vb_sub <- vb[, !is_intercept_b, drop = FALSE]
            center_b <- colMeans(vb_sub, na.rm = TRUE)
            n <- nrow(vb_sub)
            scale_b <- apply(vb_sub, 2, sd, na.rm = TRUE) * sqrt((n - 1) / n)

            zero_var_b <- scale_b < .Machine$double.eps
            scale_b[zero_var_b] <- 1

            vb_scaled_sub <- scale(vb_sub, center = center_b, scale = scale_b)
            vb[, !is_intercept_b] <- vb_scaled_sub

            scaler_b$center[!is_intercept_b] <- center_b
            scaler_b$scale[!is_intercept_b] <- scale_b
        }
    }

    list(va = va, vb = vb, scaler_a = scaler_a, scaler_b = scaler_b)
}

#' Unstandardize RBRM Coefficients
#'
#' @export
unstandardize_coeffs <- function(alpha, beta, scaler_a, scaler_b) {
    # Unstandardize alpha
    if (!is.null(alpha) && !is.null(scaler_a)) {
        alpha_orig <- alpha / scaler_a$scale
        adj <- sum(alpha * scaler_a$center / scaler_a$scale)
        alpha_orig[1] <- alpha_orig[1] - adj
        alpha <- alpha_orig
    }

    # Unstandardize beta
    if (!is.null(beta) && !is.null(scaler_b)) {
        beta_orig <- beta / scaler_b$scale
        adj_b <- sum(beta * scaler_b$center / scaler_b$scale)
        beta_orig[1] <- beta_orig[1] - adj_b
        beta <- beta_orig
    }

    list(alpha = alpha, beta = beta)
}


#' Find Lambda Max for RBRM
#'
#' @export
find_lambda_max <- function(va, vb, x, y, alpha_start = NULL, beta_start = NULL, prob_fun = getProbRR.org, intercept = TRUE) {
    n <- length(y)

    if (is.null(alpha_start)) alpha_start <- rep(0, ncol(va))
    if (is.null(beta_start)) beta_start <- rep(0, ncol(vb))

    grads <- grad_nll(alpha_start, beta_start, y, x, va, vb, prob_fun, opt = "both")

    g_alpha <- abs(grads$grad_alpha)
    g_beta <- abs(grads$grad_beta)

    if (intercept) {
        # Remove gradient for intercept
        g_alpha <- g_alpha[-1]
        g_beta <- g_beta[-1]
    }

    max_grad <- max(c(g_alpha, g_beta), na.rm = TRUE)

    return(max_grad)
}

#' Create Lambda Grid
#'
#' @export
create_lambda_grid <- function(lambda_max, nlambda = 100, lambda.min.ratio = 1e-4) {
    if (is.null(lambda_max) || lambda_max == 0) lambda_max <- 1.0 # Fallback

    lambdas <- exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio), length.out = nlambda))
    return(lambdas)
}

#' Calculate RBRM Negative Log-Likelihood
#'
#' @export
calc_nll <- function(va, vb, x, y, alpha, beta, prob_fun = getProbRR.org) {
    theta <- as.vector(va %*% alpha)
    phi <- as.vector(vb %*% beta)

    ps <- prob_fun(theta, phi)
    p0 <- ps$p0
    p1 <- ps$p1

    # Clip for stability
    ep <- 1e-10
    p0 <- pmax(p0, ep)
    p1 <- pmax(p1, ep)
    p1 <- pmin(p1, 1 - ep)

    probs <- ifelse(x == 1, p1, p0)

    ll <- y * log(probs) + (1 - y) * log(1 - probs)

    return(-mean(ll, na.rm = TRUE))
}
