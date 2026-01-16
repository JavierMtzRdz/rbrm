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
#' @param optimizer Optimization method: "fista" (default), "lbfgs", "newton", or "newton_active".
#' @export
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

    # Validate optimizer
    optimizer <- rlang::arg_match(optimizer, c(
        "fista", "lbfgs", "newton", "newton_active",
        "fista_R", "lbfgs_R", "newton_R", "newton_active_R"
    ))

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

        l_max <- find_lambda_max(std_tmp$va, std_tmp$vb, x, y, prob_fun = getProbRR.org)

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
            optimizer = optimizer,
            verbose = FALSE, ...
        )

        # Evaluate on test data
        for (i in seq_along(lambda_seq)) {
            # Coefficients are p x 1 vectors
            a_est <- path_fit$alphas[, i]
            b_est <- path_fit$betas[, i]

            # Calculate ALL metrics efficiently
            all_vals <- calc_all_measures(va_test, vb_test, x_test, y_test, a_est, b_est)

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
        optimizer = optimizer,
        verbose = FALSE, ...
    )

    # Extract coefficients for lambda.min
    result$final_fit <- list(
        path = final_path,
        alpha_min = final_path$alphas[, idx_min],
        beta_min = final_path$betas[, idx_min],
        alpha_1se = final_path$alphas[, best_idx_1se],
        beta_1se = final_path$betas[, best_idx_1se]
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

        # Count non-zeros
        nz_a <- sum(abs(x$final_fit$alpha_min) > 1e-10)
        nz_b <- sum(abs(x$final_fit$beta_min) > 1e-10)

        cli::cat_bullet("Active Va Coeffs: ", cli::col_cyan(nz_a), bullet = "arrow_right", bullet_col = "#F9C74F")
        cli::cat_bullet("Active Vb Coeffs: ", cli::col_cyan(nz_b), bullet = "arrow_right", bullet_col = "#F9C74F")
    }

    cat("\n")
    invisible(x)
}
