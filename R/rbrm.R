#' Refit Unpenalized Model on Active Set
#'
#' Helper to refit model without penalization on selected variables.
#'
#' @keywords internal
refit_unpenalized <- function(va, vb, x, y, idx_a, idx_b, intercept = FALSE, ...) {
    p <- ncol(va)
    q <- ncol(vb)

    # 1. Prepare Data Subsets logic
    if (intercept) {
        # Indices > 1 refer to columns of va
        vars_a <- idx_a[idx_a > 1] - 1
        vars_b <- idx_b[idx_b > 1] - 1
    } else {
        vars_a <- idx_a
        vars_b <- idx_b
    }

    # Subset
    va_sub <- va[, vars_a, drop = FALSE]
    vb_sub <- vb[, vars_b, drop = FALSE]

    # Fit with lambda = 0
    # Note: We assume standardize=FALSE because the passed va/vb are already standardized if needed
    fit <- fit.rbrm(va_sub, vb_sub, x, y, lambda = 0, standardize = FALSE, intercept = intercept, ...)

    # Expand Results
    if (intercept) {
        alpha_full <- rep(0, p + 1)
        beta_full <- rep(0, q + 1)

        # Assign Intercept
        alpha_full[1] <- fit$alpha[1]
        beta_full[1] <- fit$beta[1]

        # Assign Coefficients (map back to p+1 space)
        # fit$alpha has Intercept at 1, coeffs at 2..k+1
        # vars_a contains column indices (1..p)
        # alpha_full indices are vars_a + 1
        if (length(vars_a) > 0) alpha_full[vars_a + 1] <- fit$alpha[-1]
        if (length(vars_b) > 0) beta_full[vars_b + 1] <- fit$beta[-1]

        # Note: If idx_a contained 1 (Intercept), it was handled by alpha_full[1] assignment.
        # If idx_a did NOT contain 1, we still assigned the new Intercept (refitted).
        # This is generally desired behavior (refit intercept regardless of selection status, usually).
    } else {
        alpha_full <- rep(0, p)
        beta_full <- rep(0, q)

        if (length(vars_a) > 0) alpha_full[vars_a] <- fit$alpha
        if (length(vars_b) > 0) beta_full[vars_b] <- fit$beta
    }

    return(list(alpha = alpha_full, beta = beta_full))
}

#' Fit RBRM Regularization Path
#'
#' Fits the RBRM model for a sequence of lambda values using warm starts.
#' Can optionally refit unpenalized models on the active set from the penalized path.
#' If a single lambda is provided, returns a single 'rbrm' fit object.
#'
#' @param va Matrix of predictors for alpha.
#' @param vb Matrix of predictors for beta.
#' @param x Treatment vector.
#' @param y Outcome vector.
#' @param lambda Vector of lambda values (decreasing). If NULL, automatically generated.
#' @param nlambda Number of lambda values if lambda is NULL.
#' @param lambda_min_ratio Ratio of smallest to largest lambda if lambda is NULL.
#' @param alpha_start Initial alpha.
#' @param beta_start Initial beta.
#' @param standardize Logical. Whether data should be standardized (if not already).
#'  If TRUE, internal standardization is applied and returned coefficients are on original scale.
#'  If FALSE, assumes data is prepared (intercepts added, scaling done).
#' @param adjusted Logical. If TRUE, re-estimates coefficients for active variables without penalization (relaxed fit).
#' @param optimizer Optimization method: "fista" (default, C++), "lbfgs" (C++), "newton" (C++), "newton_active" (C++), or their R versions ("fista_R", "lbfgs_R", "newton_R", "newton_active_R").
#' @param ... Additional arguments to fit.rbrm.
#' @return An object of class `rbrm_path` (if multiple lambdas) or `rbrm` (if single lambda).
#' @export
rbrm <- function(va, vb, x, y, lambda = NULL,
                 nlambda = 50, lambda_min_ratio = 1e-3,
                 alpha_start = NULL, beta_start = NULL,
                 standardize = TRUE,
                 adjusted = FALSE,
                 intercept = TRUE,
                 optimizer = "fista",
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
    lambda_seq <- lambda
    if (is.null(lambda_seq)) {
        if (nlambda == 1) {
            # Special case: Auto-generate 1 lambda? Usually implies finding l_max?
            # Standard behavior might be to just generate the seq and take first?
            # Or treat as "path of length 1".
            # Let's generate as usual.
        }

        if (verbose) cli::cli_alert_info("Generating lambda sequence...")
        # Find lambda_max using standardized data (so effectively NO INTERCEPT in that matrix)
        l_max <- find_lambda_max(va, vb, x, y, prob_fun = getProbRR.org, intercept = FALSE)
        lambda_seq <- create_lambda_grid(l_max, nlambda, lambda_min_ratio)
        if (verbose) cli::cli_alert_success("Generated {length(lambda_seq)} lambdas (Max: {round(l_max, 4)})")
    }

    n_lambda <- length(lambda_seq)

    # --- SINGLE FIT MODE ---
    # Trigger if user provided single lambda OR asked for nlambda=1 (and got 1)
    if (n_lambda == 1) {
        lam <- lambda_seq[1]

        fit <- fit.rbrm(va, vb, x, y,
            alpha_start = alpha_start,
            beta_start = beta_start,
            lambda = lam,
            standardize = FALSE,
            intercept = intercept,
            optimizer = optimizer,
            ...
        )

        # Apply adjustment if requested
        if (adjusted) {
            idx_a <- which(abs(fit$alpha) > 1e-12)
            idx_b <- which(abs(fit$beta) > 1e-12)
            refit <- refit_unpenalized(va, vb, x, y, idx_a, idx_b, optimizer = optimizer, intercept = intercept, ...)

            fit$alpha <- refit$alpha
            fit$beta <- refit$beta
            fit$adjusted <- TRUE
        }

        # Unstandardize output
        if (standardize && !is.null(scaler_a)) {
            res_orig <- unstandardize_coeffs(fit$alpha, fit$beta, scaler_a, scaler_b)
            fit$alpha <- res_orig$alpha
            fit$beta <- res_orig$beta

            # Also update point.est for consistency
            fit$point.est <- c(fit$alpha, fit$beta)
            fit$va_scale_info <- scaler_a
            fit$vb_scale_info <- scaler_b
        }

        return(fit)
    }

    # --- PATH MODE ---
    path_fits <- vector("list", n_lambda)
    path_alphas <- vector("list", n_lambda)
    path_betas <- vector("list", n_lambda)

    curr_alpha <- alpha_start
    curr_beta <- beta_start

    if (verbose) params_pb <- cli::cli_progress_bar("Fitting Path", total = n_lambda)

    for (i in seq_along(lambda_seq)) {
        lam <- lambda_seq[i]

        fit <- fit.rbrm(va, vb, x, y,
            alpha_start = curr_alpha,
            beta_start = curr_beta,
            lambda = lam,
            standardize = FALSE,
            intercept = intercept,
            optimizer = optimizer,
            ...
        )

        path_fits[[i]] <- fit

        # Determine what to store (adjusted or penalized)
        if (adjusted) {
            # Find active sets
            idx_a <- which(abs(fit$alpha) > 1e-12)
            idx_b <- which(abs(fit$beta) > 1e-12)

            # Refit unpenalized
            refit <- refit_unpenalized(va, vb, x, y, idx_a, idx_b, optimizer = optimizer, intercept = intercept, ...)

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
        curr_alpha <- fit$alpha
        curr_beta <- fit$beta

        if (verbose) cli::cli_progress_update(id = params_pb)
    }

    res <- list(
        lambdas = lambda_seq,
        alphas = do.call(cbind, path_alphas), # p x n_lambda
        betas = do.call(cbind, path_betas),
        models = path_fits,
        scaler_a = scaler_a,
        scaler_b = scaler_b,
        adjusted = adjusted,
        intercept = intercept,
        dimensions = list(n = length(y), p_a = ncol(va), p_b = ncol(vb)),
        optimizer = optimizer
    )

    # Restore rownames
    if (!is.null(path_fits[[1]]$alpha)) {
        rownames(res$alphas) <- names(path_fits[[1]]$alpha)
    }
    if (!is.null(path_fits[[1]]$beta)) {
        rownames(res$betas) <- names(path_fits[[1]]$beta)
    }

    class(res) <- "rbrm_path"
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
    if (!is.null(x$adjusted) && x$adjusted) {
        cli::cat_bullet("Type: ", cli::col_cyan("Relaxed (Unpenalized Refit)"), bullet = "info", bullet_col = "#F9C74F")
    }
    if (!is.null(x$optimizer)) {
        cli::cat_bullet("Optimizer: ", cli::col_cyan(x$optimizer), bullet = "info", bullet_col = "#F9C74F")
    }

    if (!is.null(x$intercept)) {
        state <- if (x$intercept) "Included" else "Excluded"
        cli::cat_bullet("Intercept: ", cli::col_cyan(state), bullet = "info", bullet_col = "#F9C74F")
    }

    if (!is.null(x$dimensions)) {
        cli::cat_bullet("Features (Alpha): ", cli::col_cyan(x$dimensions$p_a), bullet = "info", bullet_col = "#F9C74F")
        cli::cat_bullet("Features (Beta): ", cli::col_cyan(x$dimensions$p_b), bullet = "info", bullet_col = "#F9C74F")
    }

    cat("\n")
    invisible(x)
}
