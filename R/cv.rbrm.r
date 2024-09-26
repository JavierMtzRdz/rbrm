
#' Cross-Validation for Regularized Binary Regression Model (RBRM)
#'
#' This function performs k-fold cross-validation to tune the regularization parameter (lambda) 
#' for the RBRM model. The optimal lambda is selected based on the minimum deviance.
#'
#' @param va A matrix of independent variables (without an intercept) for alpha.
#' @param vb A matrix of independent variables (without an intercept) for beta.
#' @param x A vector indicating the group assignment (1 or 0) for each observation.
#' @param y A vector of binary outcomes (0/1).
#' @param lambda Optional. A vector of lambda values for regularization. If NULL, a grid of values will be generated.
#' @param n_lambdas The number of lambda values to generate if lambda is NULL. Default is 50.
#' @param nfolds The number of folds for cross-validation. Default is 3. Must be at least 2.
#' @param ... Additional arguments passed to the `rbrm` function.
#'
#' @return A list containing:
#' \item{lambda}{The best lambda value selected based on cross-validation.}
#' \item{alpha}{The estimated alpha coefficients from the best fit model.}
#' \item{beta}{The estimated beta coefficients from the best fit model.}
#' \item{convergence}{Logical indicating whether convergence was achieved.}
#' \item{step}{The number of steps taken for convergence.}
#' \item{cv_mean_deviance}{The mean deviance across all folds for each lambda.}
#' \item{fold_deviances}{A list of deviances for each fold.}
#' \item{best_lambda_idx}{The index of the best lambda in the lambda grid.}
#' \item{model_history}{An array storing the coefficients for each lambda across folds.}
#' \item{lambda_grid}{The grid of lambda values used in cross-validation.}
#'
#' @examples
#' # Example usage:
#' #va <- matrix(c(1, 1, 1, 1, 0, 0, 0, 0), ncol = 2)
#' #vb <- matrix(c(1, 1, 0, 0, 1, 1, 0, 0), ncol = 2)
#' #y <- c(1, 0, 1, 0)
#' #x <- c(1, 1, 0, 0)
#' #cv_result <- cv_rbrm(va, vb, x, y, nfolds = 3)
#'
#' @export
cv_rbrm <- function(va, vb, x, y, lambda = NULL, n_lambdas = 50, nfolds = 3, ...) {
  # nfolds need to be at least 2
  if (nfolds < 2) {
    # stop("nfolds must be at least 2")
    cli::cli_alert_danger("nfolds must be at least 2")
  }
  n <- length(y)
  fold_ids <- sample(rep_len(1:nfolds, n))
  
  if (is.null(lambda)) {
    epsilon <- 0.001
    X <- cbind(va)
    max_lambda <- max(abs(stats::cor(X, y)))
    lambda_grid <- rev(exp(seq(log(epsilon * max_lambda),
                               log(max_lambda),
                               length.out = n_lambdas
    ))) # descending (large to small)
  } else {
    lambda_grid <- lambda
  }
  
  model_history <- array(0, dim = c(ncol(va) * 2, n_lambdas, nfolds))
  
    cli::cli_progress_bar("Processing", total = nfolds*n_lambdas)
    
  cv_results <- lapply(1:nfolds, function(fold) {
    # print(paste("Processing fold", fold, "of", nfolds, sep = " "))
    cli::cli_alert_info(paste0("Processing fold ", fold, " of ", nfolds))
    test_idx <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    
    va_train <- va[train_idx, ]
    vb_train <- vb[train_idx, ]
    x_train <- x[train_idx]
    y_train <- y[train_idx]
    
    va_test <- va[test_idx, ]
    vb_test <- vb[test_idx, ]
    x_test <- x[test_idx]
    y_test <- y[test_idx]
    
    fold_models <- lapply(1:length(lambda_grid),
                          function(x) {
      current_lambda <- lambda_grid[x]
      # cat("Processing lambda:", current_lambda, "\n")
      cli::cli_alert(paste0("Processing lambda: ", round(current_lambda, 5)))
      fit <- rbrm(va_train, vb_train, x_train, y_train, lambda = current_lambda, ...)
      model_history[, x, fold] <- fit$point.est
      cli::cli_progress_update(.envir = parent.frame(4))
      return(fit)
    })
    
    test_deviance <- sapply(fold_models, function(fit) {
      alpha <- fit$point.est[1:(length(fit$point.est) / 2)]
      beta <- fit$point.est[((length(fit$point.est) / 2) + 1):length(fit$point.est)]
      
      logrr <- va_test %*% alpha
      logop <- vb_test %*% beta
      
      p0 <- brm::getProbRR(logrr, logop)[, 1]
      p1 <- brm::getProbRR(logrr, logop)[, 2]
      
      fitted.prob <- c(p0[x_test == 0], p1[x_test == 1])
      true.y <- c(y_test[x_test == 0], y_test[x_test == 1])
      
      return((-2 * sum(log(fitted.prob[true.y == 1])) - 2 * sum(log1p(-fitted.prob[true.y == 0]))) / length(true.y))
    })
    
    return(test_deviance)
  }) 
  
  cv_mean_deviance <- rowMeans(do.call(cbind, cv_results))
  
  best_lambda_idx <- which.min(cv_mean_deviance)
  best_fit <- rbrm(va, vb, x, y, lambda = lambda_grid[best_lambda_idx], ...)
  
  list(
    lambda = lambda_grid[best_lambda_idx],
    alpha = best_fit$point.est[1:ncol(va)],
    beta = best_fit$point.est[(ncol(va) + 1):(ncol(va) + ncol(vb))],
    convergence = best_fit$convergence,
    step = best_fit$step,
    cv_mean_deviance = cv_mean_deviance,
    fold_deviances = cv_results,
    best_lambda_idx = best_lambda_idx,
    model_history = model_history,
    lambda_grid = lambda_grid
  )
}

