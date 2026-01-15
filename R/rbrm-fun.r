#' Experimental RBRM Model Fitting Function
#'
#' Wraps an optimization function to fit the RBRM model with Lasso penalty.
#'
#' @param va Matrix of predictors for alpha.
#' @param vb Matrix of predictors for beta (defaults to `va` if `NULL`).
#' @param x Binary treatment vector (0/1).
#' @param y Binary outcome vector (0/1).
#' @param alpha_start Optional starting vector for alpha coefficients.
#' @param beta_start Optional starting vector for beta coefficients.
#' @param max_step Maximum iterations for the optimizer.
#' @param lambda Lasso penalty strength (non-negative scalar).
#' @param lambda_beta Optional separate penalty for beta.
#' @param lambda_b_prop Proportion of lambda to use for beta if lambda_beta is NULL.
#' @param intercept Logical. Does the model include an intercept term?
#' @param prob_fun Function to calculate probabilities (e.g., `getProbRR.org`).
#' @param opt_fun The optimization function to use (defaults to `fista_opt`).
#' @param save_opt Logical. If `TRUE`, include the full parameter history in results.
#' @param standardize Logical. scaling.
#' @param ... Additional arguments passed to the optimization function (e.g., `step_size_alpha`, `armijo_c`).
#'
#' @seealso \code{\link{fista_opt}}, \code{\link{optim_lbfgs}}, \code{\link{optim_newton_cd}}
#' @return An object of class "rbrm".
#' @export
fit.rbrm <- function(va, vb = NULL, x, y,
                     alpha_start = NULL, beta_start = NULL,
                     max_step = 1000, lambda = 0,
                     lambda_beta = NULL,
                     lambda_b_prop = 1,
                     intercept = FALSE,
                     prob_fun = getProbRR.org,
                     opt_fun = fista_opt,
                     save_opt = FALSE,
                     standardize = TRUE,
                     clipping = 1e-10,
                     ...) {
  tictoc::tic("rbrm_fit time")

  # Argument Setup
  if (is.null(lambda_beta)) lambda_beta <- lambda * lambda_b_prop
  if (lambda < 0) {
    cli::cli_warn("lambda is negative ({lambda}), using 0 instead.")
    lambda <- 0
  }
  if (lambda_beta < 0) {
    cli::cli_warn("lambda beta is negative ({lambda.beta}), using 0 instead.")
    lambda_beta <- 0
  }
  if (is.null(vb)) vb <- va

  va <- tryCatch(as.matrix(va), error = function(e) cli::cli_abort("Failed to coerce 'va' to matrix: {e$message}"))
  vb <- tryCatch(as.matrix(vb), error = function(e) cli::cli_abort("Failed to coerce 'vb' to matrix: {e$message}"))

  n <- length(y)

  # Add Intercept Automatically
  if (intercept) {
    # Check if intercept already exists (first column is all ones)
    has_intercept_va <- isTRUE(all(va[, 1] == 1))
    has_intercept_vb <- isTRUE(all(vb[, 1] == 1))

    if (!has_intercept_va) {
      va <- cbind(Intercept = 1, va)
    }
    if (!has_intercept_vb) {
      vb <- cbind(Intercept = 1, vb)
    }
  }

  pa <- ncol(va)
  pb <- ncol(vb)

  # Standardization
  va_scaled <- va
  vb_scaled <- vb
  va_scal_info <- NULL
  vb_scal_info <- NULL

  if (standardize) {
    intercept_col_va <- NULL
    predictors_va <- va
    if (intercept) {
      intercept_col_va <- va[, 1, drop = FALSE]
      predictors_va <- va[, -1, drop = FALSE]
    }

    intercept_col_vb <- NULL
    predictors_vb <- vb
    if (intercept) {
      intercept_col_vb <- vb[, 1, drop = FALSE]
      predictors_vb <- vb[, -1, drop = FALSE]
    }

    # Scale
    scaled_preds_va <- scale(predictors_va)
    scaled_preds_vb <- scale(predictors_vb)

    va_scal_info <- list(
      center = attr(scaled_preds_va, "scaled:center"),
      scale = attr(scaled_preds_va, "scaled:scale")
    )
    vb_scal_info <- list(
      center = attr(scaled_preds_vb, "scaled:center"),
      scale = attr(scaled_preds_vb, "scaled:scale")
    )

    # Recombine intercept
    va_scaled <- if (intercept) cbind(intercept_col_va, scaled_preds_va) else scaled_preds_va
    vb_scaled <- if (intercept) cbind(intercept_col_vb, scaled_preds_vb) else scaled_preds_vb
  }

  # Initialize Starting Values
  if (is.null(alpha_start)) alpha_start <- rep(0, pa)
  if (is.null(beta_start)) beta_start <- rep(0, pb)

  # Optimizer
  # Pass explict args + ... args
  opt_args <- list(
    alpha_start = alpha_start,
    beta_start = beta_start,
    max_step = max_step,
    lambda = lambda,
    intercept = intercept,
    va = va_scaled,
    vb = vb_scaled,
    x = x,
    y = y,
    prob_fun = prob_fun,
    lambda_beta = lambda_beta,
    save_history = save_opt,
    clipping = clipping
  )

  final_args <- c(opt_args, list(...))

  opt_result <- do.call(opt_fun, final_args)

  # Extract and Back-Transform Results
  alpha_std <- opt_result$alpha
  beta_std <- opt_result$beta

  # Only access history if it exists
  alphas_std <- if (save_opt) opt_result$alphas else NULL
  betas_std <- if (save_opt) opt_result$betas else NULL

  if (standardize) {
    alpha <- alpha_std
    beta <- beta_std

    alphas <- alphas_std
    betas <- betas_std

    # Back-transform slope coefficients
    slope_indices_a <- ifelse(intercept, list(2:pa), list(1:pa))[[1]]
    slope_indices_b <- ifelse(intercept, list(2:pb), list(1:pb))[[1]]

    alpha[slope_indices_a] <- as.vector(alpha_std[slope_indices_a, drop = FALSE] / va_scal_info$scale)
    beta[slope_indices_b] <- as.vector(beta_std[slope_indices_b, drop = FALSE] / vb_scal_info$scale)

    if (save_opt) {
      alphas[, slope_indices_a] <- alphas_std[, slope_indices_a] / va_scal_info$scale
      betas[, slope_indices_b] <- betas_std[, slope_indices_b] / vb_scal_info$scale
      opt_result$alphas <- alphas
      opt_result$betas <- betas
    }

    # Adjust intercept if it exists
    if (intercept) {
      intercept_adjustment_a <- sum((alpha_std[slope_indices_a, drop = FALSE] * va_scal_info$center) / va_scal_info$scale)
      intercept_adjustment_b <- sum((beta_std[slope_indices_b, drop = FALSE] * vb_scal_info$center) / vb_scal_info$scale)

      alpha[1] <- as.numeric(alpha_std[1] - intercept_adjustment_a)
      beta[1] <- as.numeric(beta_std[1] - intercept_adjustment_b)

      if (save_opt) {
        intercept_adjustment_as <- rowSums((alphas_std[, slope_indices_a] * va_scal_info$center) / va_scal_info$scale)
        intercept_adjustment_bs <- rowSums((betas_std[, slope_indices_b] * vb_scal_info$center) / vb_scal_info$scale)
        alphas[, 1] <- alphas_std[, 1] - intercept_adjustment_as
        betas[, 1] <- betas_std[, 1] - intercept_adjustment_bs
      }
    }

    # Ensure final results are vectors, not matrices
    alpha <- as.vector(alpha)
    beta <- as.vector(beta)
  } else {
    alpha <- alpha_std
    beta <- beta_std
  }

  # Structure Output
  step <- opt_result$step
  convergence <- opt_result$convergence

  time_info <- tictoc::toc(quiet = TRUE)
  run_time <- round(time_info$toc - time_info$tic, 4)

  if (!save_opt) opt_result <- NULL

  result <- list(
    call = match.call(), point.est = c(alpha, beta), alpha = alpha,
    beta = beta, convergence = convergence, step = step,
    optimizer_details = opt_result, lambda = lambda, intercept = intercept,
    va_scale_info = va_scal_info, vb_scale_info = vb_scal_info,
    dimensions = list(n = n, p_a = pa, p_b = pb), time = run_time
  )

  return(structure(result, class = c("rbrm")))
}
