#' L-BFGS-B Optimization Wrapper for RBRM
#'
#' Drops parameters into optim usually for MLE (lambda=0).
#'
#' @param alpha_start Initial alpha coefficients.
#' @param beta_start Initial beta coefficients.
#' @param step_size_alpha Step size for alpha (ignored by L-BFGS).
#' @param step_size_beta Step size for beta (ignored by L-BFGS).
#' @param lambda Lasso penalty for alpha.
#' @param intercept Whether model includes intercept.
#' @param max_step Maximum iterations.
#' @param va Design matrix for alpha.
#' @param vb Design matrix for beta.
#' @param x Treatment vector.
#' @param y Outcome vector.
#' @param prob_fun Probability function (default getProbRR.org).
#' @param lambda_beta Penalty for beta (default NULL, uses lambda).
#' @param ... Additional arguments (e.g., clipping, save_history).
#' @export
optim_lbfgs <- function(alpha_start, beta_start,
                        step_size_alpha, step_size_beta, # Ignored
                        lambda,
                        intercept,
                        max_step,
                        va, vb, x, y,
                        prob_fun = getProbRR.org,
                        lambda_beta = NULL,
                        ...) {
    if (is.null(lambda_beta)) lambda_beta <- lambda

    if (lambda != 0 || lambda_beta != 0) {
        warning("optim_lbfgs is primarily for MLE (lambda=0). L1 penalty is non-smooth and L-BFGS-B may fail to select variables correctly.")
    }

    pa <- length(alpha_start)
    pb <- length(beta_start)
    start_par <- c(alpha_start, beta_start)

    # History storage
    save_history <- FALSE
    # Check if save_history is in ... (fit.rbrm passes save_opt as save_history if argument matches?)
    # fit.rbrm passes `save_history = save_opt` explicitly in `opt_args`.
    # So we need to add `save_history` to arguments of optim_lbfgs or capture it from ...
    args <- list(...)
    if ("save_history" %in% names(args)) save_history <- args$save_history

    history_alphas <- list()
    history_betas <- list()
    history_nll <- list()

    # Objective function
    fn <- function(par) {
        a <- par[1:pa]
        b <- par[(pa + 1):(pa + pb)]

        val <- penalized_nllh(a, b, va, vb, x, y, lambda, lambda_beta = lambda_beta, intercept = intercept, prob_fun = prob_fun)

        if (save_history) {
            # Capture history
            # Note: optim calls fn multiple times. We just append.
            history_alphas[[length(history_alphas) + 1]] <<- a
            history_betas[[length(history_betas) + 1]] <<- b
            history_nll[[length(history_nll) + 1]] <<- val
        }

        if (!is.finite(val)) {
            return(1e10)
        } # Soft barrier
        return(val)
    }

    # Gradient function
    gr <- function(par) {
        a <- par[1:pa]
        b <- par[(pa + 1):(pa + pb)]

        # Check if clipping is passed in ..., otherwise default
        args <- list(...)
        clipping <- if ("clipping" %in% names(args)) args$clipping else 1e-10

        g <- grad_nll(a, b, y, x, va, vb, prob_fun, opt = "both", method = "analytical", clipping = clipping)

        # Add L1 gradient (subgradient) if lambda > 0?
        # optim doesn't handle non-smooth. We'll just provide smooth gradient + "pseudo-gradient" of L1?
        # No, for L-BFGS on L1, we really shouldn't be here.
        # But we calculate gradient of the PENALIZED NLL.
        # subgrad |x| = sign(x).

        ga <- g$grad_alpha
        gb <- g$grad_beta

        if (lambda > 0 || lambda_beta > 0) {
            # Add L1 subgradient approximation
            # Exclude intercept from penalty if needed
            pen_idx_a <- if (intercept) 2:pa else 1:pa
            pen_idx_b <- if (intercept) 2:pb else 1:pb

            ga[pen_idx_a] <- ga[pen_idx_a] + lambda * sign(a[pen_idx_a])
            gb[pen_idx_b] <- gb[pen_idx_b] + lambda_beta * sign(b[pen_idx_b])
        }

        c(ga, gb)
    }

    # Run optim
    # Control factr: default is 1e7. For higher precision, use 1e1?
    # maxit from max_step
    res <- stats::optim(start_par, fn, gr, method = "L-BFGS-B", control = list(maxit = max_step, factr = 1e4))

    par_final <- res$par
    alpha_final <- par_final[1:pa]
    beta_final <- par_final[(pa + 1):(pa + pb)]

    alphas_ret <- if (save_history && length(history_alphas) > 0) do.call(rbind, history_alphas) else matrix(alpha_final, nrow = 1)
    betas_ret <- if (save_history && length(history_betas) > 0) do.call(rbind, history_betas) else matrix(beta_final, nrow = 1)
    nllh_ret <- if (save_history && length(history_nll) > 0) unlist(history_nll) else NULL

    return(list(
        alpha = alpha_final,
        beta = beta_final,
        convergence = (res$convergence == 0),
        step = res$counts[1], # function evaluations
        final_nll = res$value,
        # Compatibility returns
        alphas = alphas_ret,
        betas = betas_ret,
        grad_alphas = NULL, # Not saved for L-BFGS
        grad_betas = NULL, # Not saved for L-BFGS
        nllh_results = nllh_ret
    ))
}
