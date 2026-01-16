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
#' @export
cv_rbrm.default <- function(object, vb = NULL, x, y,
                            nfold = 5,
                            nlambda = 100,
                            lambda_min_ratio = ifelse(nrow(object) < ncol(object), 0.01, 0.0001),
                            lambda_seq = NULL,
                            folds = NULL,
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

    nll_mat <- matrix(NA, nrow = nfold, ncol = length(lambda_seq))

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

        path_fit <- rbrm_path(va_train, vb_train, x_train, y_train,
            lambda_seq = lambda_seq,
            standardize = TRUE,
            verbose = FALSE, ...
        )

        # Evaluate on test data
        for (i in seq_along(lambda_seq)) {
            # Coefficients are p x 1 vectors
            a_est <- path_fit$alphas[, i]
            b_est <- path_fit$betas[, i]

            nll_val <- calc_nll(va_test, vb_test, x_test, y_test, a_est, b_est)
            nll_mat[k, i] <- nll_val
        }
    }

    result <- list()
    result$lambdas <- lambda_seq
    result$nll_fold <- nll_mat
    result$nll_mean <- colMeans(nll_mat, na.rm = TRUE)
    result$nll_se <- apply(nll_mat, 2, sd, na.rm = TRUE) / sqrt(nfold)

    idx_min <- which.min(result$nll_mean)
    result$lambda_min <- lambda_seq[idx_min]

    min_nll <- result$nll_mean[idx_min]
    se_min <- result$nll_se[idx_min]

    idx_1se <- which(result$nll_mean <= min_nll + se_min)
    best_idx_1se <- min(idx_1se)
    result$lambda_1se <- lambda_seq[best_idx_1se]

    if (verbose) {
        cli::cli_alert_success("CV Complete. Min Lambda: {format(result$lambda_min, digits=4)}")
    }

    if (verbose) cli::cli_alert_info("Fitting final path on full data...")
    result$fit <- rbrm_path(va, vb, x, y, lambda_seq = lambda_seq, standardize = TRUE, verbose = FALSE, ...)

    class(result) <- "cv_rbrm"
    return(result)
}

#' @describeIn cv_rbrm Formula interface
#' @export
cv_rbrm.formula <- function(object, data, ...) {
    # TODO: Implement formula parsing to va, vb, x, y
    # For now, placeholder or basic model.matrix
    # rbrm usually requires specific va, vb setup.
    # If formula is just y ~ x1 + ... it maps to va=vb=X.

    # Extract y
    # Extract X
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    x_mat <- model.matrix(mt, mf, contrasts)
    # Remove intercept if present? rbrm usually adds it or handles it.
    # fit.rbrm logic: check for intercept col.

    # Extract treatment 'x': Wait, rbrm needs x (treatment) separate from va/vb (covariates)?
    # fit.rbrm(va, vb, x, y).
    # Standard formula `y ~ x + covs` bundles treatment into X.
    # We need to know which var is `x` (treatment).
    # The formula interface needs to specify treatment variable. Or assumes first var?
    # Let's skip simplified formula interface for now unless simple standard case.

    stop("Formula interface for cv_rbrm not fully implemented yet. Please use matrix interface.")
}

#' Print CV RBRM Object
#' @export
print.cv_rbrm <- function(x, ...) {
    cli::cat_rule(cli::style_bold("RBRM Cross-Validation"), col = "#277DA1")
    cat("\n")

    n_folds <- nrow(x$nll_fold)
    n_lam <- length(x$lambdas)

    cli::cat_bullet("Folds: ", cli::col_cyan(n_folds), bullet = "info", bullet_col = "#F9C74F")
    cli::cat_bullet("Lambda Path Length: ", cli::col_cyan(n_lam), bullet = "info", bullet_col = "#F9C74F")

    cat("\n")
    cli::cat_rule("Optimal Lambdas", col = "#43AA8B")
    cli::cat_bullet("Min Lambda: ", cli::col_cyan(sprintf("%.4f", x$lambda_min)), " (NLL: ", sprintf("%.4f", min(x$nll_mean)), ")", bullet = "star", bullet_col = "#F9C74F")
    cli::cat_bullet("1-SE Lambda: ", cli::col_cyan(sprintf("%.4f", x$lambda_1se)), bullet = "star", bullet_col = "#F9C74F")

    cat("\n")
    invisible(x)
}
