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
#' @param optimizer Optimization method: "fista" (default, C++), "lbfgs" (C++), "newton" (C++), "newton_active" (C++), or their R versions ("fista_R", "lbfgs_R", "newton_R", "newton_active_R").
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
                     intercept = TRUE,
                     prob_fun = getProbRR.org,
                     optimizer = "fista",
                     save_opt = FALSE,
                     standardize = TRUE,
                     clipping = 1e-10,
                     ...) {
  tictoc::tic("rbrm_fit time")

  # Optimizer Selection
  optimizer <- rlang::arg_match(optimizer, c(
    "fista", "lbfgs", "newton", "newton_active",
    "fista_R", "lbfgs_R", "newton_R", "newton_active_R"
  ))

  opt_fun <- switch(optimizer,
    "fista" = optim_fista_cpp,
    "lbfgs" = optim_lbfgs_cpp,
    "newton" = optim_newton_cd_cpp,
    "newton_active" = optim_newton_cd_active_cpp,
    "fista_R" = optim_fista,
    "lbfgs_R" = optim_lbfgs,
    "newton_R" = optim_newton_cd,
    "newton_active_R" = optim_newton_cd_active
  )

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
    has_intercept_va <- if (ncol(va) > 0) isTRUE(all(va[, 1] == 1)) else FALSE
    has_intercept_vb <- if (ncol(vb) > 0) isTRUE(all(vb[, 1] == 1)) else FALSE

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
  grad_alphas_std <- if (save_opt) opt_result$grad_alphas else NULL
  grad_betas_std <- if (save_opt) opt_result$grad_betas else NULL

  if (standardize) {
    alpha <- alpha_std
    beta <- beta_std

    alphas <- alphas_std
    if (!is.null(alphas) && !is.matrix(alphas)) alphas <- matrix(alphas, nrow = 1)
    betas <- betas_std
    if (!is.null(betas) && !is.matrix(betas)) betas <- matrix(betas, nrow = 1)

    # Update std references to ensure they are matrices too for later use
    alphas_std <- alphas
    betas_std <- betas
    grad_alphas <- grad_alphas_std
    grad_betas <- grad_betas_std

    # Back-transform slope coefficients
    slope_indices_a <- ifelse(intercept, list(2:pa), list(1:pa))[[1]]
    slope_indices_b <- ifelse(intercept, list(2:pb), list(1:pb))[[1]]

    alpha[slope_indices_a] <- as.vector(alpha_std[slope_indices_a, drop = FALSE] / va_scal_info$scale)
    beta[slope_indices_b] <- as.vector(beta_std[slope_indices_b, drop = FALSE] / vb_scal_info$scale)

    if (save_opt && !is.null(alphas_std)) {
      alphas[, slope_indices_a] <- alphas_std[, slope_indices_a] / va_scal_info$scale
      betas[, slope_indices_b] <- betas_std[, slope_indices_b] / vb_scal_info$scale
      opt_result$alphas <- alphas
      opt_result$betas <- betas

      # Back-transform gradients if they exist
      if (!is.null(grad_alphas_std)) {
        if (!is.matrix(grad_alphas_std)) grad_alphas_std <- matrix(grad_alphas_std, nrow = 1)
        if (!is.matrix(grad_alphas)) grad_alphas <- matrix(grad_alphas, nrow = 1)
        grad_alphas[, slope_indices_a] <- grad_alphas_std[, slope_indices_a, drop = FALSE] * va_scal_info$scale
        opt_result$grad_alphas <- grad_alphas
      }
      if (!is.null(grad_betas_std)) {
        if (!is.matrix(grad_betas_std)) grad_betas_std <- matrix(grad_betas_std, nrow = 1)
        if (!is.matrix(grad_betas)) grad_betas <- matrix(grad_betas, nrow = 1)
        grad_betas[, slope_indices_b] <- grad_betas_std[, slope_indices_b, drop = FALSE] * vb_scal_info$scale
        opt_result$grad_betas <- grad_betas
      }
    }

    # Adjust intercept if it exists
    if (intercept) {
      intercept_adjustment_a <- sum((alpha_std[slope_indices_a, drop = FALSE] * va_scal_info$center) / va_scal_info$scale)
      intercept_adjustment_b <- sum((beta_std[slope_indices_b, drop = FALSE] * vb_scal_info$center) / vb_scal_info$scale)

      alpha[1] <- as.numeric(alpha_std[1] - intercept_adjustment_a)
      beta[1] <- as.numeric(beta_std[1] - intercept_adjustment_b)

      if (save_opt && !is.null(alphas_std)) {
        intercept_adjustment_as <- rowSums(alphas_std[, slope_indices_a, drop = FALSE] * (rep(1, nrow(alphas_std)) %*% t(va_scal_info$center / va_scal_info$scale)))

        adj_vec_a <- va_scal_info$center / va_scal_info$scale
        intercept_adjustment_as <- as.vector(alphas_std[, slope_indices_a, drop = FALSE] %*% adj_vec_a)

        adj_vec_b <- vb_scal_info$center / vb_scal_info$scale
        intercept_adjustment_bs <- as.vector(betas_std[, slope_indices_b, drop = FALSE] %*% adj_vec_b)

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

  if (!save_opt) {
    opt_result$alphas <- NULL
    opt_result$betas <- NULL
    opt_result$grad_alphas <- NULL
    opt_result$grad_betas <- NULL
  }

  # Assign Names (Robust)
  nms_a <- colnames(va)
  if (is.null(nms_a)) {
    seq_a <- seq_len(length(alpha))
    nms_a <- paste0("V", seq_a)
    if (intercept) nms_a[1] <- "Intercept"
  } else {
    missing_idx <- which(nms_a == "")
    if (length(missing_idx) > 0) nms_a[missing_idx] <- paste0("V", missing_idx)
    # Ensure intercept name is standard if it is first
    if (intercept && nms_a[1] != "Intercept" && nms_a[1] != "(Intercept)") nms_a[1] <- "Intercept"
  }
  names(alpha) <- nms_a

  nms_b <- colnames(vb)
  if (is.null(nms_b)) {
    seq_b <- seq_len(length(beta))
    nms_b <- paste0("V", seq_b)
    if (intercept) nms_b[1] <- "Intercept"
  } else {
    missing_idx <- which(nms_b == "")
    if (length(missing_idx) > 0) nms_b[missing_idx] <- paste0("V", missing_idx)
    if (intercept && nms_b[1] != "Intercept" && nms_b[1] != "(Intercept)") nms_b[1] <- "Intercept"
  }
  names(beta) <- nms_b

  result <- list(
    call = match.call(), point.est = c(alpha, beta), alpha = alpha,
    beta = beta, convergence = convergence, step = step,
    optimizer_details = opt_result, lambda = lambda, intercept = intercept,
    va_scale_info = va_scal_info, vb_scale_info = vb_scal_info,
    dimensions = list(n = n, p_a = pa, p_b = pb), time = run_time
  )

  return(structure(result, class = c("rbrm")))
}

#' Print RBRM Object
#'
#' @export
print.rbrm <- function(x, ...) {
  cli::cat_rule(cli::style_bold("RBRM Model Fit"), col = "#277DA1")
  cat("\n")

  if (!is.null(x$lambda)) {
    cli::cat_bullet("Lambda: ", cli::col_cyan(sprintf("%.4f", x$lambda)), bullet = "info", bullet_col = "#F9C74F")
  }


  if (!is.null(x$dimensions)) {
    cli::cat_bullet("Features (Alpha): ", cli::col_cyan(x$dimensions$p_a), bullet = "info", bullet_col = "#F9C74F")
    cli::cat_bullet("Features (Beta): ", cli::col_cyan(x$dimensions$p_b), bullet = "info", bullet_col = "#F9C74F")
  }

  # Helper to print coeffs
  print_coefs <- function(coefs, label, intercept_used, limit = 8) {
    idx <- which(abs(coefs) > 1e-10)

    # Handle Intercept
    has_intercept <- FALSE
    match_int <- grep("^(Intercept|\\(Intercept\\))$", names(coefs))

    # If standard 'intercept' flag is true, and we find it at index 1 usually
    if (intercept_used && length(match_int) > 0) {
      # Check if intercept is non-zero
      int_idx <- match_int[1] # Take first match
      if (int_idx %in% idx) {
        has_intercept <- TRUE
        val <- coefs[int_idx]
        cli::cat_bullet(paste0(label, " Intercept: "), cli::col_cyan(sprintf("%.3f", val)), bullet = "arrow_right", bullet_col = "#F94144")
        # Rmove from index list for counting/listing
        idx <- setdiff(idx, int_idx)
      }
    }

    n_nz <- length(idx)
    cli::cat_bullet(paste0(label, " Non-zero (Penalized): "), cli::col_cyan(n_nz), bullet = "arrow_right", bullet_col = "#43AA8B")

    if (n_nz > 0) {
      show_n <- min(n_nz, limit)
      show_idx <- idx[1:show_n]
      vals <- coefs[show_idx]

      nms <- names(coefs)[show_idx]
      if (is.null(nms)) nms <- paste0("[", show_idx, "]")

      items <- paste0(nms, "=", sprintf("%.3f", vals))
      out <- paste(items, collapse = ", ")
      if (n_nz > limit) out <- paste0(out, ", ...")

      cli::cat_line(paste0("   ", cli::col_grey(out)))
    }
  }

  cat("\n")
  cli::cat_line("Coefficients:")
  print_coefs(x$alpha, "Alpha", x$intercept)
  print_coefs(x$beta, "Beta", x$intercept)

  cat("\n")
  invisible(x)
}
