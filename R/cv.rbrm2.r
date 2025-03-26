#' Calculate Performance Metrics
#'
#' Calculates deviance, MAE, and MSE for given predictions and true outcomes.
#'
#' @param alpha Estimated alpha coefficients. Can be NULL if none selected.
#' @param beta Estimated beta coefficients. Can be NULL if none selected.
#' @param va_test Validation matrix for alpha.
#' @param vb_test Validation matrix for beta.
#' @param x_test Treatment assignment vector for the test set.
#' @param y_test Outcome vector for the test set.
#' @param prob_fun Function to calculate probabilities (e.g., brm::getProbRR). Must handle matrix inputs for logrr/logop.
#' @return A named vector containing deviance, mae, and mse. Returns Inf if inputs are inconsistent.
#' @keywords internal
#' @export
calculate_metrics <- function(alpha, beta, va_test, vb_test, x_test, y_test, prob_fun) {
  # Input validation
  if (is.null(prob_fun) || !is.function(prob_fun)) {
    cli::cli_abort("prob_fun must be a valid function for calculate_metrics.")
  }
  n_test <- length(y_test)
  if (nrow(va_test) != n_test || nrow(vb_test) != n_test || length(x_test) != n_test) {
    warning("Inconsistent input dimensions for calculate_metrics.")
    return(c(deviance = Inf, mae = Inf, mse = Inf))
  }
  
  p_a_test <- ncol(va_test)
  p_b_test <- ncol(vb_test)
  
  # Check coefficient compatibility
  len_alpha <- length(alpha) # Handles NULL (length 0)
  len_beta <- length(beta)   # Handles NULL (length 0)
  
  if (len_alpha != p_a_test || len_beta != p_b_test) {
    warning(sprintf("Coefficient length mismatch in calculate_metrics. Alpha: %d vs %d cols. Beta: %d vs %d cols.",
                    len_alpha, p_a_test, len_beta, p_b_test))
    return(c(deviance = Inf, mae = Inf, mse = Inf))
  }
  
  # Calculate logrr and logop safely
  logrr <- if (p_a_test > 0) va_test %*% alpha else matrix(0, nrow = n_test, ncol = 1)
  logop <- if (p_b_test > 0) vb_test %*% beta else matrix(0, nrow = n_test, ncol = 1)
  
  # Get probabilities
  ps <- tryCatch({
    prob_fun(logrr, logop)
  }, error = function(e) {
    warning(paste("prob_fun failed during metric calculation:", e$message))
    NULL
  })
  
  if (is.null(ps) || !is.list(ps) || is.null(ps$p0) || is.null(ps$p1) ||
      length(ps$p0) != n_test || length(ps$p1) != n_test) {
    warning("prob_fun did not return expected structure or dimensions.")
    return(c(deviance = Inf, mae = Inf, mse = Inf))
  }
  
  p0 <- ps$p0
  p1 <- ps$p1
  
  # Separate observations based on x_test
  fitted.prob <- numeric(n_test)
  idx0 <- which(x_test == 0)
  idx1 <- which(x_test == 1)
  if(length(idx0) > 0) fitted.prob[idx0] <- p0[idx0]
  if(length(idx1) > 0) fitted.prob[idx1] <- p1[idx1]
  
  # Avoid log(0) issues and ensure finite probabilities
  epsilon <- 1e-15
  fitted.prob <- pmax(epsilon, pmin(1 - epsilon, fitted.prob))
  fitted.prob[!is.finite(fitted.prob)] <- 0.5 # Fallback for NaNs
  
  true.y <- y_test
  
  # Deviance calculation
  dev <- tryCatch({
    (-2 / n_test) * (sum(log(fitted.prob[true.y == 1])) +
                       sum(log1p(-fitted.prob[true.y == 0])))
  }, warning = function(w) Inf, error = function(e) Inf) # Catch log errors
  
  # MAE calculation
  mae <- mean(abs(fitted.prob - true.y))
  # MSE calculation
  mse <- mean((fitted.prob - true.y)^2)
  
  # Ensure results are finite
  dev <- ifelse(is.finite(dev), dev, Inf)
  mae <- ifelse(is.finite(mae), mae, Inf)
  mse <- ifelse(is.finite(mse), mse, Inf)
  
  return(c(deviance = dev, mae = mae, mse = mse))
}

#' Fit Relaxed Lasso Model
#'
#' Fits the model using only the active set of variables identified by
#' an initial Lasso fit, applying a relaxed penalty.
#'
#' @param va Full validation matrix for alpha.
#' @param vb Full validation matrix for beta.
#' @param x Treatment assignment vector.
#' @param y Outcome vector.
#' @param initial_fit The result object from the initial Lasso fit (used to get active set).
#' @param lambda The lambda value used for the initial fit (and scaled for the relaxed fit).
#' @param relax_factor The relaxation factor (gamma).
#' @param implt The underlying model fitting function (e.g., rbrm).
#' @param prob_fun Function to calculate probabilities.
#' @param ... Additional arguments passed to \code{implt}.
#' @return A list containing the relaxed fit object (`fit`), the indices of selected variables
#'         relative to original `va`/`vb` (`selected_vars_indices`), and the active indices
#'         for alpha and beta separately (`active_alpha_indices`, `active_beta_indices`).
#'         Returns NULL for `fit` if fitting fails or dimensions mismatch.
#' @keywords internal
fit_relaxed_lasso <- function(va, vb, x, y, initial_fit, lambda, relax_factor, implt, prob_fun, ...) {
  
  p_a <- ncol(va)
  p_b <- ncol(vb)
  
  # Check if initial fit exists and has coefficients
  if (is.null(initial_fit) || is.null(initial_fit$point.est) || length(initial_fit$point.est) != (p_a + p_b)) {
    warning("Invalid initial_fit provided to fit_relaxed_lasso.")
    # Return structure indicating failure
    return(list(fit = NULL,
                selected_vars_indices = integer(0),
                active_alpha_indices = integer(0),
                active_beta_indices = integer(0)))
  }
  initial_coeffs <- initial_fit$point.est
  
  # Identify active variables from the initial fit
  active_alpha_indices <- which(abs(initial_coeffs[1:p_a]) > 1e-10)
  active_beta_indices <- active_alpha_indices +p_a
  # active_beta_indices <- which(abs(initial_coeffs[(p_a + 1):(p_a + p_b)]) > 1e-10)
  selected_vars_indices <- c(active_alpha_indices, p_a + active_beta_indices)
  
  
  # Create subsetted matrices
  va_relaxed <- va[, active_alpha_indices, drop = FALSE]
  # vb_relaxed <- vb[, active_beta_indices, drop = FALSE]
  vb_relaxed <- va_relaxed
  n_active_alpha <- ncol(va_relaxed)
  n_active_beta <- ncol(vb_relaxed)
  expected_relaxed_coeffs_len <- n_active_alpha + n_active_beta
  
  
  # Handle cases with no selected variables
  if (expected_relaxed_coeffs_len == 0) {
    warning("Relaxed Lasso: No variables selected in the initial fit. Returning trivial fit.")
    # Create a dummy fit structure indicating no variables
    dummy_fit <- list(
      point.est = numeric(0), # No coefficients
      step = 0,
      convergence = TRUE,
      time = 0
    )
    return(list(fit = dummy_fit,
                selected_vars_indices = integer(0),
                active_alpha_indices = integer(0),
                active_beta_indices = integer(0)))
  }
  
  
  # Calculate the lambda for the relaxed fit
  relaxed_lambda <- lambda * relax_factor
  
  # Prepare arguments for the implementation function
  fit_args <- list(
    va = va_relaxed,
    vb = vb_relaxed,
    x = x,
    y = y,
    lambda = relaxed_lambda,
    ...
  )
  if (!is.null(prob_fun)) {
    fit_args$prob_fun <- prob_fun
  }
  
  # Fit the model on the active set with the relaxed penalty
  relaxed_fit <- tryCatch({
    do.call(implt, fit_args)
  }, error = function(e) {
    warning(sprintf("Relaxed Lasso fit failed (lambda_rel=%.4f): %s", relaxed_lambda, e$message))
    NULL # Return NULL on error
  })
  
  # --- Check returned coefficient length ---
  if (is.null(relaxed_fit) || is.null(relaxed_fit$point.est) || length(relaxed_fit$point.est) != expected_relaxed_coeffs_len) {
    warning(sprintf("Relaxed Lasso: `implt` returned %d coefficients, but expected %d based on active set (%d alpha, %d beta). Relaxed fit considered invalid.",
                    length(relaxed_fit$point.est), expected_relaxed_coeffs_len, n_active_alpha, n_active_beta))
    relaxed_fit <- NULL # Invalidate the fit result
  }
  # --- End Check ---
  
  return(list(fit = relaxed_fit, # This will be NULL if check failed or fit errored
              selected_vars_indices = selected_vars_indices,
              active_alpha_indices = active_alpha_indices,
              active_beta_indices = active_beta_indices))
}

#' Cross-Validate Relax Factor for Relaxed Lasso
#'
#' Performs k-fold cross-validation to select the optimal relax_factor (gamma)
#' for a relaxed Lasso model, given a fixed lambda.
#'
#' @param va Full validation matrix for alpha.
#' @param vb Full validation matrix for beta.
#' @param x Treatment assignment vector.
#' @param y Outcome vector.
#' @param fold_ids A vector indicating fold membership for each observation.
#' @param selected_lambda The single lambda value chosen from the initial CV.
#' @param implt The underlying model fitting function (e.g., rbrm).
#' @param prob_fun Function to calculate probabilities.
#' @param relax_factors_grid A numeric vector of relax factors to try (e.g., c(0, 0.25, 0.5, 0.75, 1)).
#' @param type.measure The metric to optimize ("deviance" or "mae").
#' @param nfolds Number of folds.
#' @param ... Additional arguments passed to \code{implt}.
#' @return A list containing the best relax_factor and the CV results matrix.
#' @keywords internal
cv_relax_factor <- function(va, vb, x, y, fold_ids, selected_lambda, implt, prob_fun,
                            relax_factors_grid = c(0, 0.25, 0.5, 0.75, 1),
                            type.measure = "deviance", nfolds, ...) {
  
  # browser()
  cli::cli_alert_info("Starting cross-validation for relax_factor (gamma)...")
  n_relax_factors <- length(relax_factors_grid)
  relax_cv_metrics <- array(NA, dim = c(nfolds, n_relax_factors),
                            dimnames = list(Fold = 1:nfolds, RelaxFactor = relax_factors_grid))
  
  # --- Determine Probability Function (Consistent with cv_rbrm) ---
  if (is.null(prob_fun)) {
    implt_name <- deparse(substitute(implt))
    if (implt_name == "rbrm.experimental") {
      prob_fun_to_use <- tryCatch(get("getProbRR.org"), error = function(e) NULL)
      if(is.null(prob_fun_to_use)) cli::cli_warn("Cannot find getProbRR.org for rbrm.experimental in relax factor CV")
    } else if (implt_name == "rbrm" || grepl("brm::rbrm", implt_name)) {
      prob_fun_to_use <- tryCatch(brm::getProbRR, error = function(e) NULL)
      if(is.null(prob_fun_to_use)) cli::cli_warn("Cannot find brm::getProbRR in relax factor CV.")
    } else {
      prob_fun_to_use <- NULL
      # Warning only if needed later
    }
    if(!is.null(prob_fun_to_use)) {
      prob_fun <- prob_fun_to_use
    }
  }
  # Ensure prob_fun is valid if metrics calculation definitely needs it
  if (is.null(prob_fun) && type.measure %in% c("deviance", "mae", "mse")) {
    cli::cli_abort("prob_fun is required for calculating '%s' during relax factor CV, but it could not be determined or provided.", type.measure)
  }
  
  
  # --- Initialize Progress Bar ---
  cli::cli_progress_bar(name = "Relax Factor CV", total = nfolds * n_relax_factors)
  # ---
  
  for (fold in 1:nfolds) {
    
    test_idx <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    
    va_train <- va[train_idx, , drop = FALSE]; p_a_train <- ncol(va_train)
    vb_train <- vb[train_idx, , drop = FALSE]; p_b_train <- ncol(vb_train)
    x_train <- x[train_idx]; y_train <- y[train_idx]
    
    va_test <- va[test_idx, , drop = FALSE]
    vb_test <- vb[test_idx, , drop = FALSE]
    x_test <- x[test_idx]; y_test <- y[test_idx]
    
    # 1. Fit initial Lasso model for this fold using selected_lambda
    fit_args_initial <- list(va = va_train, vb = vb_train, x = x_train, y = y_train, lambda = selected_lambda, prob_fun = prob_fun, ...)
    # if (!is.null(prob_fun)) fit_args_initial$prob_fun <- prob_fun
    # initial_fold_fit <- tryCatch({ do.call(implt, fit_args_initial) }, error = function(e) NULL)
    initial_fold_fit <- do.call(implt, fit_args_initial) 
    
    
    # Get active set
    active_alpha_indices_fold <- integer(0)
    active_beta_indices_fold <- integer(0)
    if (!is.null(initial_fold_fit) && !is.null(initial_fold_fit$point.est) && length(initial_fold_fit$point.est) == (p_a_train + p_b_train)) {
      active_alpha_indices_fold <- which(abs(initial_fold_fit$point.est[1:p_a_train]) > 1e-10)
      # active_beta_indices_fold <- which(abs(initial_fold_fit$point.est[(p_a_train + 1):(p_a_train + p_b_train)]) > 1e-10)
      active_beta_indices_fold <- active_alpha_indices_fold +p_a_train
    } else {
      warning(paste("Initial fit failed or returned unexpected structure for fold", fold, "during relax factor CV. Skipping relax factors for this fold."))
      relax_cv_metrics[fold, ] <- Inf
      cli::cli_progress_update(inc = n_relax_factors, force = TRUE)
      next
    }
    
    n_active_alpha_fold <- length(active_alpha_indices_fold)
    n_active_beta_fold <- length(active_beta_indices_fold)
    expected_relaxed_coeffs_len_fold <- n_active_alpha_fold + n_active_beta_fold
    
    if (expected_relaxed_coeffs_len_fold == 0) {
      warning(paste("No variables selected in initial fit for fold", fold, "at lambda =", selected_lambda, ". Assigning Inf metric for all relax factors."))
      relax_cv_metrics[fold, ] <- Inf
      cli::cli_progress_update(inc = n_relax_factors, force = TRUE)
      next
    }
    
    # Subset data ONCE per fold
    va_train_relaxed <- va_train[, active_alpha_indices_fold, drop = FALSE]
    # vb_train_relaxed <- vb_train[, active_beta_indices_fold, drop = FALSE]
    vb_train_relaxed <- va_train_relaxed
    va_test_relaxed <- va_test[, active_alpha_indices_fold, drop = FALSE]
    vb_test_relaxed <- va_test_relaxed
    # vb_test_relaxed <- vb_test[, active_beta_indices_fold, drop = FALSE]
    
    # 2. Iterate through relax_factors
    for (j in 1:n_relax_factors) {
      current_relax_factor <- relax_factors_grid[j]
      relaxed_lambda_fold <- selected_lambda * current_relax_factor
      
      fit_args_relaxed <- list(
        va = va_train_relaxed, vb = vb_train_relaxed,
        x = x_train, y = y_train,
        lambda = relaxed_lambda_fold,
        prob_fun = prob_fun,...
      )
      # if (!is.null(prob_fun)) fit_args_relaxed$prob_fun <- prob_fun
      
      # relaxed_fold_fit_obj <- tryCatch({
      #   do.call(implt, fit_args_relaxed)
      # }, error = function(e) {
      #   warning(sprintf("Relaxed fit failed for fold %d, factor %.2f: %s", fold, current_relax_factor, e$message))
      #   NULL
      # })
      relaxed_fold_fit_obj <- do.call(implt, fit_args_relaxed)
      
      # --- Check coefficient length from relaxed fit ---
      valid_relaxed_fit <- FALSE
      if (!is.null(relaxed_fold_fit_obj) && !is.null(relaxed_fold_fit_obj$point.est)) {
        if (length(relaxed_fold_fit_obj$point.est) == expected_relaxed_coeffs_len_fold) {
          valid_relaxed_fit <- TRUE
        } else {
          warning(sprintf("Coeff length mismatch fold %d, factor %.2f. Expected %d, Got %d.",
                          fold, current_relax_factor, expected_relaxed_coeffs_len_fold, length(relaxed_fold_fit_obj$point.est)))
        }
      }
      # --- End Check ---
      
      if (!valid_relaxed_fit) {
        relax_cv_metrics[fold, j] <- Inf
        cli::cli_progress_update(inc = 1)
        next
      }
      
      # Extract coefficients (now guaranteed to be correct length relative to relaxed data)
      relaxed_coeffs <- relaxed_fold_fit_obj$point.est
      alpha_relaxed <- if(n_active_alpha_fold > 0) relaxed_coeffs[1:n_active_alpha_fold] else NULL
      beta_relaxed <- if(n_active_beta_fold > 0) relaxed_coeffs[(n_active_alpha_fold + 1):expected_relaxed_coeffs_len_fold] else NULL
      
      # Calculate metrics on TEST data using coefficients from relaxed fit
      fold_metrics <- calculate_metrics(
        alpha = alpha_relaxed,
        beta = beta_relaxed,
        va_test = va_test_relaxed, # Use subsetted test data
        vb_test = vb_test_relaxed, # Use subsetted test data
        x_test = x_test,
        y_test = y_test,
        prob_fun = prob_fun # Already checked if NULL isn't allowed
      )
      relax_cv_metrics[fold, j] <- fold_metrics[type.measure]
      
      cli::cli_progress_update(inc = 1)
      
    } # End loop over relax factors
  } # End loop over folds
  
  cli::cli_progress_done()
  
  # Aggregate results
  cv_mean_relax <- colMeans(relax_cv_metrics, na.rm = TRUE)
  cv_mean_relax[is.nan(cv_mean_relax) | !is.finite(cv_mean_relax)] <- Inf
  
  best_relax_factor_idx <- which.min(cv_mean_relax)
  if (length(best_relax_factor_idx) == 0 || !is.finite(cv_mean_relax[best_relax_factor_idx])) {
    warning("Could not determine best relax_factor from CV. Defaulting to 1.")
    best_relax_factor <- 1.0
  } else {
    best_relax_factor <- relax_factors_grid[best_relax_factor_idx]
  }
  
  cli::cli_alert_success(paste("Selected best relax_factor:", best_relax_factor))
  
  return(list(best_relax_factor = best_relax_factor,
              cv_results = relax_cv_metrics,
              cv_mean = cv_mean_relax))
}


# Ensure necessary libraries are loaded
# library(cli)
# library(tictoc)
# library(stats)
# Assuming 'rbrm' and probability functions are available
#' @export
cv_rbrm2 <- function(va, vb, x, y, lambda = NULL,
                    n_lambdas = 20, nfolds = 3,
                    implt = rbrm, # Make sure 'rbrm' is defined or loaded
                    prob_fun = NULL,
                    relax_lsso = FALSE,
                    relax_factor = NULL, # Default NULL triggers CV
                    relax_factors_grid = c(0, 0.25, 0.5, 0.75, 1), # Grid for CV
                    index = "min",       # "min" or "1se"
                    type.measure = "deviance", # "deviance", "mae", or "mse"
                    ...) {
  
  tictoc::tic("Total time")
  
  # --- Input Checks ---
  if (nfolds < 2) cli::cli_abort("nfolds must be at least 2")
  if (!type.measure %in% c("deviance", "mae", "mse")) {
    cli::cli_abort("type.measure must be one of 'deviance', 'mae', or 'mse'")
  }
  # Basic dimension checks
  n <- length(y)
  if(nrow(va) != n || nrow(vb) != n || length(x) != n) cli::cli_abort("Input dimensions mismatch (va, vb, x, y rows/lengths).")
  p_a <- ncol(va)
  p_b <- ncol(vb)
  if(is.null(colnames(va)) && p_a > 0) colnames(va) <- paste0("va_", 1:p_a)
  if(is.null(colnames(vb)) && p_b > 0) colnames(vb) <- paste0("vb_", 1:p_b)
  
  fold_ids <- sample(rep_len(1:nfolds, n))
  
  # --- Lambda Grid Setup ---
  # (Using simplified lambda max estimation - adjust if needed)
  if (is.null(lambda)) {
    lambda_grid <- tryCatch({
      X_combined <- cbind(va, vb)
      variances <- apply(X_combined, 2, var, na.rm = TRUE)
      X_combined_valid <- X_combined[, variances > 1e-8, drop = FALSE]
      if (ncol(X_combined_valid) == 0) stop("No variance in predictors.")
      # Ensure y has variation if it's binary
      if(length(unique(y)) < 2) stop("Outcome y has no variation.")
      
      # Use glmnet's approach for lambda max estimation (more robust for logistic-like)
      # Requires glmnet package: if (!requireNamespace("glmnet", quietly = TRUE)) stop("glmnet package needed for lambda grid estimation.")
      # fit_pilot <- glmnet::glmnet(X_combined_valid, y, family = "binomial", alpha = 1, nlambda = 1)
      # max_lambda <- fit_pilot$lambda[1]
      # Simple correlation as fallback (less ideal for logistic)
      cor_vals <- abs(stats::cor(X_combined_valid, y, use = "pairwise.complete.obs"))
      max_lambda_est <- max(cor_vals, na.rm = TRUE)
      if(!is.finite(max_lambda_est) || max_lambda_est == 0) max_lambda_est = 1
      
      max_lambda <- max_lambda_est # Use directly or add margin? Adjust heuristic.
      epsilon <- 0.001
      rev(exp(seq(log(epsilon * max_lambda), log(max_lambda), length.out = n_lambdas)))
    }, error = function(e) {
      warning("Automatic lambda grid estimation failed: ", e$message, ". Using default grid.")
      exp(seq(log(0.001), log(1), length.out = n_lambdas)) # Default fallback
    })
    n_lambdas <- length(lambda_grid)
  } else {
    lambda_grid <- sort(lambda, decreasing = TRUE)
    n_lambdas <- length(lambda_grid)
  }
  cli::cli_alert_info(paste("Using", n_lambdas, "lambdas. Range:", signif(min(lambda_grid), 3), "to", signif(max(lambda_grid), 3)))
  
  
  # --- Determine Probability Function ---
  prob_fun_determined <- NULL
  if (is.null(prob_fun)) {
    implt_name <- deparse(substitute(implt))
    if (implt_name == "rbrm.experimental") {
      prob_fun_determined <- tryCatch(get("getProbRR.org"), error = function(e) NULL)
      if(is.null(prob_fun_determined)) cli::cli_warn("Cannot find getProbRR.org for rbrm.experimental")
    } else if (implt_name == "rbrm" || grepl("brm::rbrm", implt_name)) {
      prob_fun_determined <- tryCatch(brm::getProbRR, error = function(e) NULL)
      if(is.null(prob_fun_determined)) cli::cli_warn("Cannot find brm::getProbRR.")
    } else {
      cli::cli_warn("Cannot automatically determine prob_fun for the provided implt.")
    }
    if(!is.null(prob_fun_determined)) {
      prob_fun <- prob_fun_determined # Assign if found
      cli::cli_alert_info(paste("Using determined prob_fun:", deparse(substitute(prob_fun_determined))))
    }
  } else {
    cli::cli_alert_info("Using user-provided prob_fun.")
    if (!is.function(prob_fun)) cli::cli_abort("Provided prob_fun is not a function.")
  }
  # Final check if prob_fun needed but unavailable
  if (is.null(prob_fun) && type.measure %in% c("deviance", "mae", "mse")) {
    cli::cli_abort("prob_fun is required for type.measure '%s', but it could not be determined or provided.", type.measure)
  }
  
  
  # --- Primary Cross-Validation for Lambda ---
  model_results <- vector("list", nfolds)
  cv_metrics_list <- vector("list", nfolds) # Stores 3xN matrices
  
  cli::cli_progress_bar("CV for Lambda", total = nfolds * n_lambdas)
  
  for (fold in 1:nfolds) {
    test_idx <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    
    va_train <- va[train_idx, , drop = FALSE]; p_a_train <- ncol(va_train)
    vb_train <- vb[train_idx, , drop = FALSE]; p_b_train <- ncol(vb_train)
    x_train <- x[train_idx]; y_train <- y[train_idx]
    
    va_test <- va[test_idx, , drop = FALSE]
    vb_test <- vb[test_idx, , drop = FALSE]
    x_test <- x[test_idx]; y_test <- y[test_idx]
    
    fold_models <- vector("list", n_lambdas)
    fold_metrics_raw <- matrix(NA, nrow = 3, ncol = n_lambdas,
                               dimnames = list(c("deviance", "mae", "mse"), lambda_grid))
    
    for (i in 1:n_lambdas) {
      current_lambda <- lambda_grid[i]
      fit_args <- list(va = va_train, vb = vb_train, x = x_train, y = y_train, lambda = current_lambda, ...)
      if (!is.null(prob_fun)) fit_args$prob_fun <- prob_fun
      
      fit <- tryCatch({ do.call(implt, fit_args) }, error = function(e) {
        warning(sprintf("Fold %d, lambda %.5f fit error: %s", fold, current_lambda, e$message))
        NULL
      })
      
      fold_models[[i]] <- fit
      valid_fit <- !is.null(fit) && !is.null(fit$point.est) && length(fit$point.est) == (p_a_train + p_b_train)
      
      if (valid_fit) {
        alpha_fit <- fit$point.est[1:p_a_train]
        beta_fit <- fit$point.est[(p_a_train + 1):(p_a_train + p_b_train)]
        metrics <- calculate_metrics(alpha_fit, beta_fit, va_test, vb_test, x_test, y_test, prob_fun)
        fold_metrics_raw[, i] <- metrics
      } else {
        if(!is.null(fit)) warning(sprintf("Fold %d, lambda %.5f coeff length mismatch. Expected %d, Got %d.",
                                          fold, current_lambda, p_a_train + p_b_train, length(fit$point.est)))
        fold_metrics_raw[, i] <- c(deviance = Inf, mae = Inf, mse = Inf)
      }
      cli::cli_progress_update()
    } # End lambda loop
    model_results[[fold]] <- fold_models
    cv_metrics_list[[fold]] <- fold_metrics_raw
  } # End fold loop
  cli::cli_progress_done()
  
  # --- Aggregate Lambda CV Results ---
  cv_results_matrix <- do.call(rbind, lapply(cv_metrics_list, function(m) m[type.measure, ]))
  cv_mean <- colMeans(cv_results_matrix, na.rm = TRUE)
  cv_sd   <- apply(cv_results_matrix, 2, sd, na.rm = TRUE)
  cv_mean[is.nan(cv_mean) | is.na(cv_mean)] <- Inf
  cv_sd[is.nan(cv_sd) | is.na(cv_sd)] <- 0
  cv_se <- cv_sd / sqrt(nfolds)
  
  # --- Select Lambda ---
  best_lambda_idx <- which.min(cv_mean)
  if (length(best_lambda_idx) == 0 || !is.finite(cv_mean[best_lambda_idx])) {
    cli::cli_abort("CV failed: No finite optimal lambda found for type.measure='%s'. Check model fits and metrics.", type.measure)
  }
  lambda_min <- lambda_grid[best_lambda_idx]
  threshold <- cv_mean[best_lambda_idx] + cv_se[best_lambda_idx]
  valid_idx <- which(cv_mean <= threshold + 1e-8) # Add tolerance
  lambda_1se <- max(lambda_grid[valid_idx], na.rm = TRUE) # Max lambda within 1SE
  
  if (tolower(index) == "min") {
    lambda_selected <- lambda_min
    cli::cli_alert_info(paste("Selected lambda (min):", signif(lambda_selected, 4)))
  } else if (tolower(index) == "1se") {
    lambda_selected <- lambda_1se
    cli::cli_alert_info(paste("Selected lambda (1se):", signif(lambda_selected, 4)))
  } else {
    warning("Invalid index '", index, "', defaulting to 'min'.")
    lambda_selected <- lambda_min
    cli::cli_alert_info(paste("Selected lambda (default min):", signif(lambda_selected, 4)))
  }
  
  # --- Refit Model(s) on Full Data ---
  cli::cli_alert_info("Refitting model on full data with selected lambda...")
  fit_args_full <- list(va = va, vb = vb, x = x, y = y, lambda = lambda_selected, ...)
  if (!is.null(prob_fun)) fit_args_full$prob_fun <- prob_fun
  initial_full_fit <- tryCatch({ do.call(implt, fit_args_full) }, error = function(e) NULL)
  
  if (is.null(initial_full_fit) || is.null(initial_full_fit$point.est) || length(initial_full_fit$point.est) != (p_a + p_b)) {
    cli::cli_abort("Failed to fit model on full data with lambda=%.4f. Check implt.", lambda_selected)
  }
  
  # --- Relaxed Lasso Step (Optional) ---
  final_fit <- initial_full_fit # Default to initial fit
  relax_info <- list( # Store relaxation details
    performed = FALSE,
    factor_used = NULL,
    factor_source = "N/A", # "Provided" or "CV"
    factor_cv_results = NULL,
    reconstruction_status = "N/A" # "Success", "Failed (Length Mismatch)", "Failed (Fit Error)"
  )
  
  if (relax_lsso) {
    relax_info$performed <- TRUE
    cli::cli_alert_info("Performing relaxed Lasso step...")
    actual_relax_factor <- NULL
    relax_cv_output <- NULL
    
    # Determine the relax factor
    if (is.null(relax_factor)) {
      relax_info$factor_source <- "CV"
      relax_cv_output <- cv_relax_factor(
        va = va, vb = vb, x = x, y = y, fold_ids = fold_ids,
        selected_lambda = lambda_selected, implt = implt, prob_fun = prob_fun,
        relax_factors_grid = relax_factors_grid, type.measure = type.measure,
        nfolds = nfolds, ...
      )
      actual_relax_factor <- relax_cv_output$best_relax_factor
      relax_info$factor_cv_results <- relax_cv_output # Store CV details
      cli::cli_alert_info(paste("Relax factor chosen by CV:", actual_relax_factor))
    } else {
      relax_info$factor_source <- "Provided"
      if (!is.numeric(relax_factor) || relax_factor < 0 || relax_factor > 1) {
        warning("Provided relax_factor invalid. Using 1.0.")
        actual_relax_factor <- 1.0
      } else {
        actual_relax_factor <- relax_factor
      }
      cli::cli_alert_info(paste("Using provided relax factor:", actual_relax_factor))
    }
    relax_info$factor_used <- actual_relax_factor
    
    
    # Fit the relaxed Lasso model using the initial full fit results
    relaxed_result <- fit_relaxed_lasso(
      va = va, vb = vb, x = x, y = y,
      initial_fit = initial_full_fit,
      lambda = lambda_selected,
      relax_factor = actual_relax_factor,
      implt = implt, prob_fun = prob_fun, ...
    )
    
    # --- Robust Coefficient Reconstruction ---
    # Check if relaxed fit itself was valid (returned by fit_relaxed_lasso)
    if (!is.null(relaxed_result$fit)) {
      # relaxed_result$fit$point.est length check was already done inside fit_relaxed_lasso
      relaxed_coeffs_vec <- relaxed_result$fit$point.est
      n_active_alpha <- length(relaxed_result$active_alpha_indices)
      n_active_beta <- length(relaxed_result$active_beta_indices)
      
      # Prepare final coefficient vector
      final_coeffs_relaxed <- vector("numeric", p_a + p_b)
      names(final_coeffs_relaxed) <- c(colnames(va), colnames(vb))
      
      # Assign coefficients
      if (n_active_alpha > 0) {
        final_coeffs_relaxed[relaxed_result$active_alpha_indices] <- relaxed_coeffs_vec[1:n_active_alpha]
      }
      if (n_active_beta > 0) {
        final_coeffs_relaxed[relaxed_result$active_beta_indices] <- relaxed_coeffs_vec[(n_active_alpha + 1):length(relaxed_coeffs_vec)]
      }
      
      # Update final_fit structure with relaxed results
      final_fit$point.est <- final_coeffs_relaxed
      final_fit$convergence <- relaxed_result$fit$convergence
      final_fit$step <- relaxed_result$fit$step
      final_fit$time <- initial_full_fit$time + ifelse(!is.null(relaxed_result$fit$time), relaxed_result$fit$time, 0)
      relax_info$reconstruction_status <- "Success"
      cli::cli_alert_success("Relaxed Lasso fitting and reconstruction complete.")
      
    } else {
      # Relaxed fit failed (either error or length mismatch reported by fit_relaxed_lasso)
      warning("Final relaxed Lasso step failed or produced inconsistent results. Using non-relaxed coefficients.")
      relax_info$reconstruction_status <- "Failed (See Previous Warnings)"
      # final_fit remains the initial_full_fit in this case
    }
    # --- End Robust Reconstruction ---
  } # End if(relax_lsso)
  
  # --- Prepare Output ---
  time_elapsed <- tictoc::toc(quiet = TRUE)
  
  # Final safety check on coefficient vector length before creating object

  if (length(final_fit$point.est) != p_a + p_b) {
    cli::cli_abort("FATAL: Final coefficient vector length (%d) inconsistent with original dimensions (%d). Check implt function and reconstruction logic.",
                   length(final_fit$point.est), p_a + p_b)
  }
  
  obj <- list(
    call = match.call(),
    lambda_grid = lambda_grid,
    cv_mean = cv_mean,
    cv_se = cv_se,
    type.measure = type.measure,
    lambda.min = lambda_min,
    lambda.1se = lambda_1se,
    lambda.selected = lambda_selected,
    index = index,
    alpha = final_fit$point.est[1:p_a],
    beta = final_fit$point.est[(p_a + 1):(p_a + p_b)],
    convergence = final_fit$convergence, # Reflects final model used
    step = final_fit$step,             # Reflects final model used
    relax_lsso_info = relax_info, # Contains all relaxation details
    # Store full CV details if needed for plotting etc.
    cv_metrics_all_folds = cv_metrics_list, # list of matrices (3 x n_lambdas)
    cv_results_matrix = cv_results_matrix, # nfolds x n_lambdas matrix for selected type.measure
    fold_ids = fold_ids,
    # model_details_per_fold = model_results, # List (folds) of lists (lambdas) of model fits (can be very large!) - uncomment if needed
    final_fit_object = final_fit, # Contains point.est etc. of the FINAL model used
    time = round(time_elapsed$toc - time_elapsed$tic, 4)
  )
  
  class(obj) <- "cv_rbrm2"
  cli::cli_alert_success(paste("Total cv_rbrm execution time:", obj$time, "seconds"))
  return(obj)
}

