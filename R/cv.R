#' Cross-Validation for RBRM
#'
#' Performs k-fold cross-validation for the RBRM model.
#'
#' @param object Formula or matrix.
#' @param ... Additional arguments.
#' @export
cv_rbrm <- function(object, ...) {
    UseMethod("cv_rbrm")
}

#' @describeIn cv_rbrm Default method for matrices
#' @param measure Performance measure: "deviance" (default), "brier", "misclass", or "auc"
#' @param adjusted Logical. If TRUE (default), the final model is refitted without penalization on the active set.
#' @param optimizer Optimization method: "fista" (default, C++), "lbfgs" (C++), "newton" (C++), "newton_active" (C++), or their R versions ("fista_R", "lbfgs_R", "newton_R", "newton_active_R").
#' @export
#' @param vb Matrix for beta (default NULL, uses object).
#' @param x Treatment vector.
#' @param y Outcome vector.
#' @param nfold Number of folds (default 5).
#' @param nlambda Number of lambda values (default 50).
#' @param lambda_min_ratio Minimum lambda ratio (default based on n/p).
#' @param lambda_seq Custom lambda sequence (default NULL, auto-generated).
#' @param folds Custom fold assignments (default NULL, auto-generated).
#' @param alpha_start Initial alpha coefficients (default NULL).
#' @param beta_start Initial beta coefficients (default NULL).
#' @param intercept Whether to include intercept (default TRUE).
#' @param seed Random seed for fold generation (default NULL).
#' @param verbose Print progress messages (default TRUE).
#' @param data Data frame for formula interface.

#' @param vb Matrix for beta (default NULL, uses object).
#' @param x Treatment vector.
#' @param y Outcome vector.
#' @param nfold Number of folds (default 5).
#' @param nlambda Number of lambda values (default 50).
#' @param lambda_min_ratio Minimum lambda ratio (default based on n/p).
#' @param lambda_seq Custom lambda sequence (default NULL, auto-generated).
#' @param folds Custom fold assignments (default NULL, auto-generated).
#' @param alpha_start Initial alpha coefficients (default NULL).
#' @param beta_start Initial beta coefficients (default NULL).
#' @param intercept Whether to include intercept (default TRUE).
#' @param seed Random seed for fold generation (default NULL).
#' @param verbose Print progress messages (default TRUE).
#' @param data Data frame for formula interface.

cv_rbrm.default <- function(object, vb = NULL, x, y,
                            nfold = 5,
                            nlambda = 50,
                            lambda_min_ratio = ifelse(nrow(object) < ncol(object), 0.01, 0.001),
                            lambda_seq = NULL,
                            folds = NULL,
                            measure = "deviance",
                            adjusted = TRUE,
                            optimizer = "fista",
                            alpha_start = NULL, beta_start = NULL,
                            intercept = TRUE,
                            seed = NULL,
                            verbose = TRUE,
                            ...) {
    va <- object
    if (is.null(vb)) vb <- va

    if (!is.numeric(nfold) || length(nfold) != 1 || nfold < 2) {
        cli::cli_abort("{.arg nfold} must be an integer >= 2.")
    }
    if (!is.numeric(nlambda) || length(nlambda) != 1 || nlambda < 2) {
        cli::cli_abort("{.arg nlambda} must be an integer >= 2.")
    }


    n <- length(y)
    if (nrow(va) != n) cli::cli_abort("{.arg va} must have the same number of rows as length of {.arg y}.")
    if (nrow(vb) != n) cli::cli_abort("{.arg vb} must have the same number of rows as length of {.arg y}.")
    if (length(x) != n) cli::cli_abort("{.arg x} must have the same length as {.arg y}.")

    if (!all(y %in% c(0, 1))) cli::cli_abort("{.arg y} must contain only 0 and 1.")
    if (!all(x %in% c(0, 1))) cli::cli_warn("{.arg x} ideally contains only 0 and 1 (treatment indicator).")

    if (!is.null(lambda_seq)) {
        if (!is.numeric(lambda_seq) || any(lambda_seq < 0)) cli::cli_abort("{.arg lambda_seq} must be numeric and non-negative.")
    }

    if (!is.null(seed)) set.seed(seed)

    # Lambda Sequence Generation
    if (is.null(lambda_seq)) {
        if (verbose) cli::cli_alert_info("Generating lambda sequence...")
        std_tmp <- standardize_data(va, vb)

        # Find lambda_max using standardized data (so effectively NO INTERCEPT in that matrix)
        l_max <- find_lambda_max(std_tmp$va, std_tmp$vb, x, y, prob_fun = getProbRR.org, intercept = FALSE)

        lambda_seq <- create_lambda_grid(l_max, nlambda, lambda_min_ratio)
        if (verbose) cli::cli_alert_success("Generated {length(lambda_seq)} lambdas (Max: {round(l_max, 4)})")
    }

    if (is.null(folds)) {
        folds <- split(sample(seq(n)), rep(1:nfold, length = n))
    } else {
        nfold <- length(folds)
    }

    # Get performance measure name
    measure_name <- get_measure_name(measure)

    if (verbose) cli::cli_alert_info("Using {measure_name} as primary selection metric")

    n_lambda <- length(lambda_seq)

    # Storage for all metrics
    metrics <- c("deviance", "brier", "misclass", "auc", "f1")
    perf_arrays <- list()
    for (m in metrics) {
        perf_arrays[[m]] <- matrix(NA, nrow = nfold, ncol = n_lambda)
    }

    if (verbose) cli::cli_progress_bar("Running Cross-Validation", total = nfold)

    for (k in 1:nfold) {
        if (verbose) cli::cli_progress_update()

        idx_test <- folds[[k]]
        idx_train <- setdiff(1:n, idx_test)

        va_train <- va[idx_train, , drop = FALSE]
        vb_train <- vb[idx_train, , drop = FALSE]
        x_train <- x[idx_train]
        y_train <- y[idx_train]

        va_test <- va[idx_test, , drop = FALSE]
        vb_test <- vb[idx_test, , drop = FALSE]
        x_test <- x[idx_test]
        y_test <- y[idx_test]

        # Use rbrm (renamed from rbrm_path)
        # CV always runs on UNADJUSTED (penalized) models for selection speed/consistency
        path_fit <- rbrm(va_train, vb_train, x_train, y_train,
            lambda = lambda_seq,
            standardize = TRUE,
            adjusted = FALSE,
            intercept = intercept,
            optimizer = optimizer,
            verbose = FALSE, ...
        )

        # Evaluate on test data
        for (i in seq_along(lambda_seq)) {
            # Coefficients are p x 1 vectors
            a_est <- path_fit$alphas[, i]
            b_est <- path_fit$betas[, i]

            # Match Dimensions for Test Set (Add Intercept if needed)
            va_test_calc <- va_test
            vb_test_calc <- vb_test

            # If intercept=TRUE, refitted model has Intercept column.
            # calc_all_measures expects matrix matching coefficients.
            if (intercept) {
                # Add intercept column manually to test set
                if (!any(va_test_calc[, 1] == 1)) va_test_calc <- cbind(1, va_test_calc)
                if (!any(vb_test_calc[, 1] == 1)) vb_test_calc <- cbind(1, vb_test_calc)
            }

            # Calculate ALL metrics efficiently
            all_vals <- calc_all_measures(va_test_calc, vb_test_calc, x_test, y_test, a_est, b_est)

            for (m in metrics) {
                perf_arrays[[m]][k, i] <- all_vals[[m]]
            }
        }
    }

    result <- list()
    result$lambdas <- lambda_seq
    result$measure <- measure
    result$measure_name <- measure_name
    result$adjusted <- adjusted
    result$optimizer <- optimizer

    # Store full results for all metrics
    result$cv_results <- list()
    for (m in metrics) {
        mat <- perf_arrays[[m]]
        res <- list(
            mean = colMeans(mat, na.rm = TRUE),
            se = apply(mat, 2, sd, na.rm = TRUE) / sqrt(nfold)
        )
        result$cv_results[[m]] <- res
    }

    # Primary metric results (for compatibility and default print/plot)
    selected_res <- result$cv_results[[measure]]
    result$performance_mean <- selected_res$mean
    result$performance_se <- selected_res$se
    result$performance <- perf_arrays[[measure]] # Keep matrix for the primary one

    # Calculate Min and 1SE based on PRIMARY metric
    idx_min <- which.min(result$performance_mean)
    result$lambda_min <- lambda_seq[idx_min]

    min_perf <- result$performance_mean[idx_min]
    se_min <- result$performance_se[idx_min]

    idx_1se <- which(result$performance_mean <= min_perf + se_min)
    best_idx_1se <- min(idx_1se)
    result$lambda_1se <- lambda_seq[best_idx_1se]

    if (verbose) {
        cli::cli_alert_success("CV Complete. Min Lambda: {format(result$lambda_min, digits=4)}")
        status_msg <- if (adjusted) "Fitting adjusted final model (unpenalized refit)..." else "Fitting final model on full dataset..."
        cli::cli_alert_info("{status_msg}")
    }

    # Fit Final Model on Full Data
    # Apply adjustment if requested
    final_path <- rbrm(va, vb, x, y,
        lambda = lambda_seq,
        standardize = TRUE,
        adjusted = adjusted,
        intercept = intercept,
        optimizer = optimizer,
        verbose = FALSE, ...
    )

    # Construct "rbrm" object for optimal fit to mimic fit.rbrm output
    final_fit_obj <- list(
        call = match.call(),
        point.est = c(final_path$alphas[, idx_min], final_path$betas[, idx_min]),
        alpha = final_path$alphas[, idx_min],
        beta = final_path$betas[, idx_min],
        lambda = result$lambda_min,
        intercept = intercept,
        va_scale_info = final_path$scaler_a,
        vb_scale_info = final_path$scaler_b,
        dimensions = list(n = nrow(va), p_a = ncol(va), p_b = ncol(vb)),
        optimizer_details = NULL,
        time = NULL,
        convergence = TRUE,
        step = NA
    )
    class(final_fit_obj) <- "rbrm"

    result$path <- final_path
    result$final_fit <- final_fit_obj

    # Store 1se fit as simpler list (or full object if desired, sticking to list for now)
    result$fit_1se <- list(
        alpha = final_path$alphas[, best_idx_1se],
        beta = final_path$betas[, best_idx_1se],
        lambda = result$lambda_1se
    )

    class(result) <- "cv_rbrm"
    return(invisible(result))
}
#' @describeIn cv_rbrm Formula interface
#' @export
cv_rbrm.formula <- function(object, data, ...) {
    # TODO: Implement formula parsing to va, vb, x, y
    stop("Formula interface for cv_rbrm not fully implemented yet. Please use matrix interface.")
}

#' Print CV RBRM Object
#' @export
#' @param x A cv_rbrm object.
#' @param ... Additional arguments (unused).

print.cv_rbrm <- function(x, ...) {
    cli::cat_rule(cli::style_bold("RBRM Cross-Validation"), col = "#277DA1")
    cat("\n")

    n_folds <- nrow(x$performance)
    n_lam <- length(x$lambdas)
    measure_name <- if (!is.null(x$measure_name)) x$measure_name else "Performance"

    cli::cat_bullet("Folds: ", cli::col_cyan(n_folds), bullet = "info", bullet_col = "#F9C74F")
    cli::cat_bullet("Lambda Path Length: ", cli::col_cyan(n_lam), bullet = "info", bullet_col = "#F9C74F")
    cli::cat_bullet("Measure: ", cli::col_cyan(measure_name), bullet = "info", bullet_col = "#F9C74F")
    cli::cat_bullet("Refit Unpenalized: ", cli::col_cyan(if (x$adjusted) "Yes" else "No"), bullet = "info", bullet_col = "#F9C74F")
    if (!is.null(x$optimizer)) {
        cli::cat_bullet("Optimizer: ", cli::col_cyan(x$optimizer), bullet = "info", bullet_col = "#F9C74F")
    }

    cat("\n")
    cli::cat_rule("Optimal Lambdas", col = "#43AA8B")
    cli::cat_bullet("Min Lambda: ", cli::col_cyan(sprintf("%.4f", x$lambda_min)),
        " (", measure_name, ": ", sprintf("%.2f", min(x$performance_mean)), ")",
        bullet = "star", bullet_col = "#F9C74F"
    )
    cli::cat_bullet("1-SE Lambda: ", cli::col_cyan(sprintf("%.4f", x$lambda_1se)), bullet = "star", bullet_col = "#F9C74F")

    if (!is.null(x$final_fit)) {
        cat("\n")
        cli::cat_rule("Final Model (at Lambda Min)", col = "#43AA8B")

        # Get coeffs (handle rbrm object or legacy list)
        alpha <- if (!is.null(x$final_fit$alpha)) x$final_fit$alpha else x$final_fit$alpha_min
        beta <- if (!is.null(x$final_fit$beta)) x$final_fit$beta else x$final_fit$beta_min

        # Check intercept
        intercept <- if (!is.null(x$final_fit$intercept)) x$final_fit$intercept else x$intercept

        # Helper to count
        count_nz <- function(c, int) {
            # If intercept is TRUE, assume first element is intercept
            if (isTRUE(int) && length(c) > 0) c <- c[-1]
            sum(abs(c) > 1e-10)
        }

        nz_a <- count_nz(alpha, intercept)
        nz_b <- count_nz(beta, intercept)

        cli::cat_bullet("Active Va Coeffs (Penalized): ", cli::col_cyan(nz_a), bullet = "arrow_right", bullet_col = "#F9C74F")
        cli::cat_bullet("Active Vb Coeffs (Penalized): ", cli::col_cyan(nz_b), bullet = "arrow_right", bullet_col = "#F9C74F")
    }

    cat("\n")
    invisible(x)
}

