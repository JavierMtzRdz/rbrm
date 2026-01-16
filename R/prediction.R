

#' @title Calculate Prediction Errors for Log Risk Ratios
#' @description Computes the root mean square error (RMSE) and median absolute error (MAE) between true and estimated log risk ratios based on input covariates and coefficients.
#' @param va A matrix or data frame of covariates with dimensions \eqn{n \times p}, where each row represents an observation and each column represents a covariate.
#' @param true_alpha A numeric vector of true coefficients with length equal to the number of columns in `va`.
#' @param est_alpha A numeric vector of estimated coefficients with length equal to the number of columns in `va`.
#' @details The function calculates two error metrics between the true and estimated log risk ratios (LRR):
#' \itemize{
#'   \item{RMSE: Root mean square error of the differences between true and estimated LRRs.}
#'   \item{MAE: Median absolute error of the differences between true and estimated LRRs.}
#' }
#' These metrics help assess the accuracy of estimated coefficients in predicting risk ratios.
#' @return A list containing:
#' \describe{
#'   \item{mae}{Median absolute error of the predicted LRR.}
#'   \item{rmse}{Root mean square error of the predicted LRR.}
#' }
#' 
#' @export
get_rr_predict_err <- function(va, true_alpha, est_alpha) {
  true_lrr <- va %*% true_alpha # dim = n x 1
  est_lrr <- va %*% est_alpha # dim = n x 1
  
  rmse <- sqrt(mean((true_lrr - est_lrr)^2)) # one value for every dataset/replicate
  # mape <- mean(abs((true_rr - est_rr) / true_rr)) # one value for every dataset/replicate
  # mae <- mean(abs(true_lrr - est_lrr)) # one value for every dataset/replicate
  mae <- stats::median(abs(true_lrr - est_lrr)) # one value for every dataset/replicate
  return(list(mae = mae, rmse = rmse))
}
