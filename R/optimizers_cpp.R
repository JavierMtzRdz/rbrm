#' Newton-Coordinate Descent Optimizer (C++ Implementation)
#'
#' Standard Newton-CD using C++ backend.
#'
#' @inheritParams optim_newton_cd
#' @export
optim_newton_cd_cpp <- function(alpha_start, beta_start,
                                step_size_alpha, step_size_beta, # Ignored
                                lambda,
                                intercept,
                                max_step,
                                va, vb, x, y,
                                prob_fun = getProbRR.org,
                                lambda_beta = NULL,
                                tol = 1e-5,
                                save_history = FALSE,
                                ...) {
    # Determine prob_fun_selector
    prob_fun_name <- deparse(substitute(prob_fun))
    if (is.function(prob_fun)) {
        if (identical(prob_fun, getProbRR.alt)) {
            selector <- 1
        } else {
            selector <- 0
        }
    } else {
        selector <- 0
    }

    if (is.null(lambda_beta)) lambda_beta <- -1.0

    args <- list(...)
    clipping <- if ("clipping" %in% names(args)) args$clipping else 1e-10

    # Call C++
    # newton_cd_cpp is exported from RcppExports
    res <- newton_cd_cpp(
        alpha_start, beta_start, lambda, intercept, max_step,
        as.matrix(va), as.matrix(vb), x, y,
        selector, lambda_beta, tol, clipping
    )

    return(res)
}

#' Active Set Newton-Coordinate Descent Optimizer (C++ Implementation)
#'
#' Improves performance for sparse problems by updating only active set and checking KKT conditions periodically.
#'
#' @inheritParams optim_newton_cd
#' @inheritParams optim_newton_cd_active
#' @export
optim_newton_cd_active_cpp <- function(alpha_start, beta_start,
                                       step_size_alpha, step_size_beta, # Ignored
                                       lambda,
                                       intercept,
                                       max_step,
                                       va, vb, x, y,
                                       prob_fun = getProbRR.org,
                                       lambda_beta = NULL,
                                       tol = 1e-5,
                                       save_history = FALSE, # Not implemented in C++ yet (ignored)
                                       kkt_check_freq = 10,
                                       active_tol = 1e-6,
                                       ...) {
    # Determine prob_fun_selector
    prob_fun_name <- deparse(substitute(prob_fun))
    # Handle if prob_fun is passed as value
    if (is.function(prob_fun)) {
        if (identical(prob_fun, getProbRR.alt)) {
            selector <- 1
        } else {
            selector <- 0 # Default to Richardson
        }
    } else {
        selector <- 0
    }

    if (is.null(lambda_beta)) lambda_beta <- -1.0 # Signal to use lambda

    args <- list(...)
    clipping <- if ("clipping" %in% names(args)) args$clipping else 1e-10

    # Call C++
    res <- active_set_newton_cd_cpp(
        alpha_start, beta_start, lambda, intercept, max_step,
        as.matrix(va), as.matrix(vb), x, y,
        selector, lambda_beta, tol, clipping,
        kkt_check_freq, active_tol
    )

    return(res)
}

#' FISTA Optimizer (C++ Implementation)
#'
#' Fast Iterative Shrinkage-Thresholding Algorithm (FISTA) using C++ backend.
#'
#' @inheritParams optim_newton_cd
#' @param step_size_init Initial step size for backtracking line search.
#' @param armijo_c Armijo condition parameter (default 1e-4). Not used in current C++ implementation (quadratic approx used instead).
#' @param shrink_factor Factor to reduce step size during backtracking (default 0.5).
#' @export
optim_fista_cpp <- function(alpha_start, beta_start,
                            step_size_alpha = 0.5, step_size_beta = NULL,
                            lambda,
                            intercept,
                            max_step,
                            va, vb, x, y,
                            prob_fun = getProbRR.org,
                            lambda_beta = NULL,
                            tol = 1e-5,
                            save_history = FALSE,
                            step_size_init = 0.5,
                            armijo_c = 1e-4,
                            shrink_factor = 0.5,
                            ...) {
    # Determine prob_fun_selector
    prob_fun_name <- deparse(substitute(prob_fun))
    if (is.function(prob_fun)) {
        if (identical(prob_fun, getProbRR.alt)) {
            selector <- 1
        } else {
            selector <- 0
        }
    } else {
        selector <- 0
    }

    if (is.null(lambda_beta)) lambda_beta <- -1.0 # Signal to use lambda

    args <- list(...)
    clipping <- if ("clipping" %in% names(args)) args$clipping else 1e-10

    res <- fista_cpp(
        alpha_start, beta_start, lambda, intercept, max_step,
        as.matrix(va), as.matrix(vb), x, y,
        selector, lambda_beta, tol,
        step_size_init, armijo_c, shrink_factor, clipping
    )

    return(res)
}

#' L-BFGS Optimizer (C++ Implementation)
#'
#' L-BFGS optimization using R's internal C API (lbfgsb).
#'
#' @inheritParams optim_newton_cd
#' @export
optim_lbfgs_cpp <- function(alpha_start, beta_start,
                            step_size_alpha = NULL, step_size_beta = NULL,
                            lambda,
                            intercept,
                            max_step,
                            va, vb, x, y,
                            prob_fun = getProbRR.org,
                            lambda_beta = NULL,
                            tol = 1e-5,
                            save_history = FALSE,
                            ...) {
    if (is.null(lambda_beta)) lambda_beta <- lambda

    # Determine prob_fun_selector
    prob_fun_name <- deparse(substitute(prob_fun))
    if (is.function(prob_fun)) {
        if (identical(prob_fun, getProbRR.alt)) {
            selector <- 1
        } else {
            # Default to 0 for org or others
            selector <- 0
        }
    } else {
        selector <- 0
    }

    args <- list(...)
    clipping <- if ("clipping" %in% names(args)) args$clipping else 1e-10

    if (lambda != 0 || lambda_beta != 0) {
        warning("L-BFGS is primarily for smooth functions (MLE). Using it with L1 penalty may yield suboptimal results.")
    }

    res <- lbfgs_cpp(
        alpha_start, beta_start, lambda, intercept, max_step,
        as.matrix(va), as.matrix(vb), x, y,
        selector, lambda_beta, tol, clipping
    )

    return(res)
}
