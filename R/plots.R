#' Plot Cross-Validation Results for Regularized Models
#'
#' Generates a plot to visualize cross-validation results for a regularized model. 
#' Displays the mean deviance with standard errors for a grid of regularization 
#' parameters (`lambda`) and marks the optimal `lambda` value with a dashed line.
#'
#' @param .model A list containing:
#'   - `fold_deviances`: A structure (e.g., list of numeric vectors) containing deviance values for each fold.
#'   - `lambda_grid`: A numeric vector specifying the grid of `lambda` values.
#'   - `lambda`: The optimal `lambda` value selected based on minimum deviance.
#' @return A `ggplot` object visualizing:
#'   - Mean deviance values across folds for each `lambda`.
#'   - Standard error bars for deviance values.
#'   - A vertical dashed line indicating the selected `lambda`.
#'
#' @examples
#' #library(rbrm)
#' ## Example model
#' #model <- list(
#' #  fold_deviances = list(
#' #    c(1.2, 1.35, 1.3, 1.3, 1.1),
#' #    c(1.4, 1.20, 1.2, 1.2, 1.3),
#' #    c(1.1, 1.05, 1.0, 1.05, 1.2)
#' #  ),
#' #  lambda_grid = c(0.01, 0.03, 0.1, 0.3, 1),
#' #  lambda = 0.1
#' #)
#' #
#' ## Plot cross-validation results
#' #plot.cv_rbrm(model)
#'
#' @importFrom magrittr %>%
#' @export
plot.cv_rbrm <- function(.model, type.measure = "mae") {
  result <- do.call(cbind, lapply(.model$cv_results,
                                   function(m) m[type.measure, ])) 
  
  lambda_grid <- .model$lambda_grid
  # Choose which measure to use
  
  nfolds <- length(.model$cv_results)
  cv_mean <- rowMeans(result)
  cv_sd   <- apply(result, 1, sd)
  cv_se <- cv_sd / sqrt(nfolds)
  
  # Lambda selection: "min" or "1se"
  best_lambda_idx <- which.min(cv_mean)
  threshold <- cv_mean[best_lambda_idx] + cv_se[best_lambda_idx]
  valid_idx <- which(cv_mean <= threshold)
  lambda_1se <- max(lambda_grid[valid_idx])
  
  tibble::as_tibble(result) %>%
    dplyr::mutate(lambda = .model$lambda_grid) %>%
    tidyr::pivot_longer(-lambda, names_to = "fold") %>%
    ggplot2::ggplot(ggplot2::aes(x = lambda, y = value, group = lambda)) +
    ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar",
                          width = 0.05, colour = "#277DA1") +
    ggplot2::stat_summary(fun = mean, geom = "point",
                          colour = "#277DA1", alpha = 1) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = lambda_grid[best_lambda_idx], ,
                                     linetype = "Lambda min"),
                        colour = "#f94144") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = lambda_1se,
                                     linetype = "Lambda 1SE"),
                        colour = "#f94144") +
    ggplot2::scale_x_log10(n.breaks = 6) +
    ggplot2::scale_linetype_manual(values = c("dotted", "dashed")) +
    ggplot2::labs(x = "Lambda", y = type.measure,
                  linetype = ggplot2::element_blank(),
                  colour = ggplot2::element_blank()) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")
}
