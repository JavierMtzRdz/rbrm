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
plot.cv_rbrm <- function(.model, type.measure = c("deviance", "auc", "accuracy",
                                                  "sensitivity", "specificity",
                                                  "precision",
                                                  "f1_score")) {
  type.measure <- rlang::arg_match(type.measure)
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
                          width = 0.01, colour = "#277DA1") +
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
#' @export
plot.cv_rbrm2 <- function(.model,
                          type.measure = NULL, opt = "lambda") {
  
  if(tolower(opt) == "lambda") {
    if(ifelse(is.null(type.measure), TRUE, tolower(type.measure) == .model$type.measure)) {
      data <- tibble::as_tibble(.model$cv_results_matrix)
      min <- .model$lambda.min
    } else {
      type.measure <- tolower(type.measure)
      data <- purrr::map_df(.model$cv_metrics_all_folds, ~.x[type.measure,])
      means <- rowMeans(data)
      min <- names(data)[which.min(means)]
    }}
  
  type.measure <- rlang::arg_match0(type.measure,
                                    c("deviance", "auc", "accuracy",
                                      "sensitivity", "specificity",
                                      "precision",
                                      "f1_score"))
  
  if(tolower(opt) %in% c("relax_factor", 
                         "gamma")) {
    data <- tibble::as_tibble(.model$relax_lsso_info$factor_cv_results$cv_results)
    min <- .model$relax_lsso_info$factor_cv_results$best_relax_factor
    
    if(!is.null(type.measure)) if (tolower(type.measure) != .model$type.measure) cli::cli_abort("{type.measure} is not available for gamma.")
    
  }
  opt_uppper <- paste(toupper(substring(opt, 1,1)), tolower(substring(opt, 2)),
                      sep="")
  
  data %>%
    tidyr::pivot_longer(everything(), names_to = "opt") %>%
    dplyr::mutate(opt = as.numeric(opt)) %>% 
    ggplot2::ggplot(ggplot2::aes(x = opt, y = value, group = opt)) +
    ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar",
                          width = 0.01, colour = "#277DA1") +
    ggplot2::stat_summary(fun = mean, geom = "point",
                          colour = "#277DA1", alpha = 1) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = min,
                                     linetype = paste(opt_uppper, "min")),
                        colour = "#f94144") +
    {if(ifelse(tolower(opt) == "lambda", TRUE, 
               tolower(type.measure) == .model$type.measure)) {
      ggplot2::geom_vline(
      ggplot2::aes(xintercept = .model$lambda.1se,
                   linetype = "Lambda 1SE"),
      colour = "#f94144")}} +
    {if(tolower(opt) == "lambda") {ggplot2::scale_x_log10(n.breaks = 6)}} +
    {if(tolower(opt) == "gamma") {ggplot2::scale_x_continuous(limits = c(-0.1, 1.1))}} +
    ggplot2::scale_linetype_manual(values = c("dotted", "dashed")) +
    ggplot2::labs(x = opt_uppper, y = .model$type.measure,
                  linetype = ggplot2::element_blank(),
                  colour = ggplot2::element_blank()) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")
}


#' Prepare Data for Plotting cv_rbrm2 Object
#'
#' Extracts and summarizes cross-validation results for plotting.
#'
#' @param x A fitted object of class \code{cv_rbrm2}.
#' @param opt Character string: "lambda" (default) or "gamma"/"relax_factor".
#' @param type.measure Optional character string: metric to plot if different
#'   from the one used in fitting (only applies if opt = "lambda").
#'
#' @return A list containing: `plot_data` (tibble with parameter, mean, se, etc.),
#'   `param_min` (value of lambda.min or best gamma), `param_1se` (value of
#'   lambda.1se or NA), `param_lab` (label for x-axis), `measure_lab` (label
#'   for y-axis), `opt_lower` (validated lower-case opt). Returns NULL elements
#'   if data is unavailable.
#' @keywords internal
prepare_plot_data_cv_rbrm2 <- function(x, opt = "lambda", type.measure = NULL) {
  # --- Validate Inputs ---
  if (!inherits(x, "cv_rbrm2")) {
    cli::cli_abort("Input object must be of class 'cv_rbrm2'.")
  }
  opt_lower <- tolower(opt)
  if (!(opt_lower %in% c("lambda", "gamma", "relax_factor"))) {
    cli::cli_abort("Argument 'opt' must be 'lambda', 'gamma', or 'relax_factor'.")
  }
  
  # --- Initialize Outputs ---
  plot_data <- NULL
  cv_data_matrix <- NULL
  param_grid <- NULL
  param_min <- NA_real_
  param_1se <- NA_real_
  param_lab <- ""
  measure_lab <- x$type.measure # Default to original measure
  
  # --- Extract Data Based on `opt` ---
  if (opt_lower == "lambda") {
    param_lab <- "Lambda"
    param_grid <- x$lambda_grid
    
    # Determine which metric to use
    plot_measure <- x$type.measure # Default
    if (!is.null(type.measure)) {
      user_measure <- tolower(type.measure)
      # Check if user-requested measure exists in the stored detailed metrics
      if (!user_measure %in% rownames(x$cv_metrics_all_folds[[1]])) {
        cli::cli_warn("Requested type.measure '{user_measure}' not found in cv_metrics_all_folds. Using original '{x$type.measure}'.")
      } else if (user_measure != x$type.measure) {
        cli::cli_alert_info("Plotting user-specified type.measure: '{user_measure}'.")
        plot_measure <- user_measure
        # Need to recalculate means/SEs for this measure
        cv_data_matrix <- tryCatch({
          do.call(rbind, lapply(x$cv_metrics_all_folds, function(m) m[plot_measure, ]))
        }, error = function(e){ NULL })
        if(is.null(cv_data_matrix)) cli::cli_warn("Could not extract data for type.measure '{plot_measure}'.")
        # Also recalculate min/1se based on this measure
        temp_lambda_sel <- select_lambda(cv_data_matrix, param_grid, plot_measure, nrow(cv_data_matrix), "min") # Use helper
        param_min <- temp_lambda_sel$lambda.min
        param_1se <- temp_lambda_sel$lambda.1se # select_lambda already calculates this
        
      } else {
        # User requested the same measure as original, use stored results
        cv_data_matrix <- x$cv_results_matrix
        param_min <- x$lambda.min
        param_1se <- x$lambda.1se
      }
    } else {
      # type.measure is NULL, use original stored results
      cv_data_matrix <- x$cv_results_matrix
      param_min <- x$lambda.min
      param_1se <- x$lambda.1se
    }
    measure_lab <- plot_measure # Set label to the measure actually being plotted
    
  } else { # opt is "gamma" or "relax_factor"
    param_lab <- "Relax Factor (Gamma)"
    
    # Check if relaxation was performed and results exist
    if (!x$relax_lsso_info$performed || is.null(x$relax_lsso_info$factor_cv_results) || is.null(x$relax_lsso_info$factor_cv_results$cv_results)) {
      cli::cli_alert_warning("Relaxation factor CV results not found in the object. Cannot plot '{opt}'.")
      return(list(plot_data=NULL)) # Return NULL list elements
    }
    
    # Check consistency of type.measure
    if (!is.null(type.measure) && tolower(type.measure) != x$type.measure) {
      cli::cli_abort("Cannot plot gamma/relax_factor for a different 'type.measure' ('{type.measure}') than used during fitting ('{x$type.measure}').")
    }
    measure_lab <- x$type.measure # Use original measure
    
    cv_data_matrix <- x$relax_lsso_info$factor_cv_results$cv_results
    # Ensure colnames exist and are numeric for param_grid
    if(is.null(colnames(cv_data_matrix))) {
      cli::cli_warn("Relax factor CV results matrix missing column names. Cannot plot.")
      return(list(plot_data=NULL))
    }
    param_grid <- as.numeric(colnames(cv_data_matrix))
    param_min <- x$relax_lsso_info$factor_cv_results$best_relax_factor
    param_1se <- NA_real_ # No 1SE rule for relax factor CV here
  }
  
  # --- Check if Data Matrix is Valid ---
  if (is.null(cv_data_matrix) || !is.matrix(cv_data_matrix) || nrow(cv_data_matrix) == 0 || ncol(cv_data_matrix) == 0) {
    cli::cli_alert_warning("No valid CV data matrix found for plotting.")
    return(list(plot_data=NULL))
  }
  if (length(param_grid) != ncol(cv_data_matrix)) {
    cli::cli_alert_warning("Parameter grid length does not match CV results matrix columns.")
    return(list(plot_data=NULL))
  }
  colnames(cv_data_matrix) <- as.character(param_grid) # Ensure names match grid
  
  # --- Calculate Summary Statistics ---
  nfolds_actual <- nrow(cv_data_matrix)
  plot_data <- cv_data_matrix %>%
    as.data.frame() %>%
    # Add fold ID if needed: tibble::rownames_to_column(var = "fold") %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(), # Pivot all columns (parameters)
      names_to = "parameter_str",
      values_to = "metric_value",
      values_drop_na = TRUE # Drop NA metric values from folds/lambdas that failed
    ) %>%
    dplyr::mutate(
      parameter = as.numeric(.data$parameter_str) # Convert param string back to numeric
    ) %>%
    dplyr::filter(is.finite(.data$metric_value)) %>% # Exclude Inf values before summarizing
    dplyr::group_by(.data$parameter) %>%
    dplyr::summarise(
      mean = mean(.data$metric_value), # Already removed NA/Inf
      sd = sd(.data$metric_value),
      n = dplyr::n(), # Number of non-NA/Inf folds for this parameter
      se = ifelse(.data$n > 1, .data$sd / sqrt(.data$n), 0), # Calculate SE, handle n=1
      upper = .data$mean + .data$se,
      lower = .data$mean - .data$se,
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$n > 0) # Keep only parameters with at least one valid fold result
  
  if(nrow(plot_data) == 0){
    cli::cli_alert_warning("No finite results available to plot after summary.")
    return(list(plot_data=NULL))
  }
  
  # Return prepared data and parameters
  list(
    plot_data = plot_data,
    param_min = param_min,
    param_1se = param_1se,
    param_lab = param_lab,
    measure_lab = measure_lab,
    opt_lower = opt_lower
  )
}

# --- Refactored plot.cv_rbrm2 ---

#' Plot Cross-Validation Curve from cv_rbrm2 Object
#'
#' Creates a plot of the cross-validation curve(s) showing the mean performance
#' metric versus the tuning parameter (lambda or relax factor).
#'
#' @param x A fitted object of class \code{cv_rbrm2}.
#' @param type.measure Optional character string: metric to plot if different
#'   from the one used in fitting (only applies if \code{opt = "lambda"}).
#' @param opt Character string: Parameter to plot against. Either "lambda" (default)
#'   or "gamma"/"relax_factor".
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{ggplot} object.
#'
#' @method plot cv_rbrm2
#' @import ggplot2
#' @importFrom dplyr %>% mutate group_by summarise filter n
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom rlang .data
#' @export
#' @examples
#' \dontrun{
#' # Assuming 'cv_fit' is a result from cv_rbrm2
#' plot(cv_fit)
#' plot(cv_fit, opt = "gamma") # If relaxation CV was performed
#' plot(cv_fit, type.measure = "mae") # Plot MAE instead of original metric (if lambda)
#' }
plot.cv_rbrm2 <- function(x, type.measure = NULL, opt = "lambda", ...) {
  
  # --- 1. Prepare Data ---
  prep_results <- prepare_plot_data_cv_rbrm2(x = x, opt = opt, type.measure = type.measure)
  
  plot_data <- prep_results$plot_data
  param_min <- prep_results$param_min
  param_1se <- prep_results$param_1se
  param_lab <- prep_results$param_lab
  measure_lab <- prep_results$measure_lab
  opt_lower <- prep_results$opt_lower
  
  # Check if data preparation was successful
  if (is.null(plot_data)) {
    cli::cli_alert_info("Cannot generate plot due to missing or invalid CV data for selected options.")
    return(invisible(NULL)) # Return NULL invisibly
  }
  # browser()
  # --- 2. Base Plot ---
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$parameter, y = .data$mean)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      width = 0.01, # Adjust width as needed, maybe relative?
      colour = "#277DA1"
    ) +
    ggplot2::geom_point(colour = "#277DA1")
  
  # --- 3. Add Vertical Lines for Optimal Values ---
  plot_captions <- character() # Store captions for lines
  
  if (!is.na(param_min)) {
    p <- p + ggplot2::geom_vline(
      aes(linetype = paste0(param_lab, ".min"),
          xintercept = param_min),
      # linetype = "dashed", # Dashed for minimum
      colour = "#F94144")
  }
  
  # Add 1se line only if plotting lambda and value is valid
  if (opt_lower == "lambda" && !is.na(param_1se)) {
    p <- p + ggplot2::geom_vline(
      aes(linetype = paste0(param_lab, ".1se"),
          xintercept = param_1se),
      # linetype = "dotted", # Dotted for 1se
      colour = "#F94144" # Different color
    )
  }
  
  
  # --- 4. Apply Scales (Log for Lambda) ---
  if (opt_lower == "lambda") {
    # Use pseudo-log if data might be near zero? No, lambda > 0.
    # Add more breaks, customize labels if needed
    p <- p + ggplot2::scale_x_log10()
  } else { # relax_factor / gamma
    # Ensure axis covers 0 to 1 nicely
    limits_x <- c(min(0, min(plot_data$parameter, na.rm=TRUE) - 0.05),
                  max(1, max(plot_data$parameter, na.rm=TRUE) + 0.05))
    p <- p + ggplot2::scale_x_continuous(limits = limits_x, breaks = seq(0, 1, by=0.25))
  }
  
  # --- 5. Labels and Theme ---
  ltp <- setNames(c("dashed", "dotted"),
                  c(paste0(param_lab, ".min"), paste0(param_lab, ".1se")))
  p <- p + ggplot2::labs(
    x = param_lab,
    y = tools::toTitleCase(measure_lab), # Capitalize metric name
    title = paste("Cross-Validation Results for", param_lab),
    linetype = element_blank()) +
    ggplot2::scale_linetype_manual(values = ltp)  +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "t0p", # No legend needed for fixed linetypes/colors
      plot.caption = ggplot2::element_text(hjust = 0, size = 8, colour="gray30") # Align caption left
    )
  
  return(p)
}
