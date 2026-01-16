#' Refit Unpenalized Model on Active Set
#'
#' Helper to refit model without penalization on selected variables.
#'
#' @keywords internal
refit_unpenalized <- function(va, vb, x, y, idx_a, idx_b, ...) {
    # If sets are empty, return zeros (or intercepts if they were handled? rbrm structure assumes columns)
    p <- ncol(va)
    q <- ncol(vb)

    # Check if totally empty
    if (length(idx_a) == 0 && length(idx_b) == 0) {
        return(list(alpha = rep(0, p), beta = rep(0, q)))
    }

    # Subset
    va_sub <- va[, idx_a, drop = FALSE]
    vb_sub <- vb[, idx_b, drop = FALSE]

    # Fit with lambda = 0
    # Note: We assume standardize=FALSE because the passed va/vb are already standardized if needed
    fit <- fit.rbrm(va_sub, vb_sub, x, y, lambda = 0, standardize = FALSE, ...)

    # Expand
    alpha_full <- rep(0, p)
    beta_full <- rep(0, q)

    if (length(idx_a) > 0) alpha_full[idx_a] <- fit$alpha
    if (length(idx_b) > 0) beta_full[idx_b] <- fit$beta

    return(list(alpha = alpha_full, beta = beta_full))
}

#' Fit RBRM Regularization Path
#'
#' Fits the RBRM model for a sequence of lambda values using warm starts.
#' Can optionally refit unpenalized models on the active set from the penalized path.
#'
#' @param va Matrix of predictors for alpha.
#' @param vb Matrix of predictors for beta.
#' @param x Treatment vector.
#' @param y Outcome vector.
#' @param lambda_seq Vector of lambda values (decreasing). If NULL, automatically generated.
#' @param nlambda Number of lambda values if lambda_seq is NULL.
#' @param lambda_min_ratio Ratio of smallest to largest lambda if lambda_seq is NULL.
#' @param alpha_start Initial alpha.
#' @param beta_start Initial beta.
#' @param standardize Logical. Whether data should be standardized (if not already).
#'  If TRUE, internal standardization is applied and returned coefficients are on original scale.
#'  If FALSE, assumes data is prepared (intercepts added, scaling done).
#' @param adjusted Logical. If TRUE, re-estimates coefficients for active variables without penalization (relaxed fit).
#' @param ... Additional arguments to fit.rbrm.
#' @return An object of class `rbrm_path` (keeping class name for compatibility).
#' @export
rbrm <- function(va, vb, x, y, lambda_seq = NULL,
                 nlambda = 100, lambda_min_ratio = 1e-4,
                 alpha_start = NULL, beta_start = NULL,
                 standardize = TRUE,
                 adjusted = FALSE,
                 verbose = FALSE, ...) {
    if (is.null(vb)) vb <- va
    n <- length(y)

    if (nrow(va) != n) cli::cli_abort("{.arg va} rows ({nrow(va)}) must match length of {.arg y} ({n}).")
    if (nrow(vb) != n) cli::cli_abort("{.arg vb} rows ({nrow(vb)}) must match length of {.arg y} ({n}).")
    if (length(x) != n) cli::cli_abort("{.arg x} length ({length(x)}) must match length of {.arg y} ({n}).")
    if (!all(y %in% c(0, 1))) cli::cli_warn("{.arg y} should ideally be 0/1.")

    scaler_a <- NULL
    scaler_b <- NULL

    if (standardize) {
        std <- standardize_data(va, vb)
        va <- std$va
        vb <- std$vb
        scaler_a <- std$scaler_a
        scaler_b <- std$scaler_b
    }

    # Auto-generate lambda sequence if not provided
    if (is.null(lambda_seq)) {
        if (verbose) cli::cli_alert_info("Generating lambda sequence...")
        l_max <- find_lambda_max(va, vb, x, y, prob_fun = getProbRR.org)
        lambda_seq <- create_lambda_grid(l_max, nlambda, lambda_min_ratio)
        if (verbose) cli::cli_alert_success("Generated {length(lambda_seq)} lambdas (Max: {round(l_max, 4)})")
    }

    n_lambda <- length(lambda_seq)

    path_fits <- vector("list", n_lambda)
    path_alphas <- vector("list", n_lambda)
    path_betas <- vector("list", n_lambda)

    curr_alpha <- alpha_start
    curr_beta <- beta_start

    if (verbose) cli::cli_progress_bar("Fitting Path", total = n_lambda)

    for (i in seq_along(lambda_seq)) {
        lam <- lambda_seq[i]

        fit <- fit.rbrm(va, vb, x, y,
            alpha_start = curr_alpha,
            beta_start = curr_beta,
            lambda = lam,
            standardize = FALSE,
            ...
        )

        path_fits[[i]] <- fit

        # Determine what to store (adjusted or penalized)
        if (adjusted) {
            # Find active sets (using small threshold for robustness against numerical noise)
            idx_a <- which(abs(fit$alpha) > 1e-12)
            idx_b <- which(abs(fit$beta) > 1e-12)

            # Refit unpenalized (using current standardized data)
            refit <- refit_unpenalized(va, vb, x, y, idx_a, idx_b, ...)

            store_alpha <- refit$alpha
            store_beta <- refit$beta
        } else {
            store_alpha <- fit$alpha
            store_beta <- fit$beta
        }

        # Store unstandardized version for output
        if (standardize && !is.null(scaler_a)) {
            res_orig <- unstandardize_coeffs(store_alpha, store_beta, scaler_a, scaler_b)
            path_alphas[[i]] <- res_orig$alpha
            path_betas[[i]] <- res_orig$beta
        } else {
            path_alphas[[i]] <- store_alpha
            path_betas[[i]] <- store_beta
        }

        # Warm start ALWAYS updates from the PENALIZED solution
        # This preserves the regularization path characteristics
        curr_alpha <- fit$alpha
        curr_beta <- fit$beta

        if (verbose) cli::cli_progress_update()
    }

    res <- list(
        lambdas = lambda_seq,
        alphas = do.call(cbind, path_alphas), # p x n_lambda
        betas = do.call(cbind, path_betas),
        models = path_fits,
        scaler_a = scaler_a,
        scaler_b = scaler_b
    )

    class(res) <- "rbrm_path" # Keep attribute as is
    return(res)
}

#' Print RBRM Path
#' @export
print.rbrm_path <- function(x, ...) {
    cli::cat_rule(cli::style_bold("RBRM Regularization Path"), col = "#277DA1")
    cat("\n")

    if (!is.list(x) || is.null(x$lambdas)) {
        cli::cli_alert_danger("Invalid 'rbrm_path' object.", col = "#f94144")
        return(invisible(x))
    }

    n_lam <- length(x$lambdas)

    cli::cat_bullet("Lambdas: ", cli::col_cyan(n_lam), bullet = "info", bullet_col = "#F9C74F")
    cli::cat_bullet("Range: ", cli::col_cyan(sprintf("%.4f - %.4f", min(x$lambdas), max(x$lambdas))), bullet = "info", bullet_col = "#F9C74F")

    cat("\n")
    invisible(x)
}
