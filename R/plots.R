#' Plot Cross-Validation Results
#'
#' @param x A 'cv_rbrm' object.
#' @param measure Optional character string to plot a specific metric (e.g., "auc", "brier").
#'                If NULL, plots the primary metric used in CV.
#' @param ... Additional arguments (unused).
#' @return A ggplot object.
#' @export
#' @param x A cv_rbrm object.
#' @param measure Performance measure to plot (default NULL).
#' @param ... Additional arguments (unused).

#' @param x A cv_rbrm object.
#' @param measure Performance measure to plot (default NULL).
#' @param ... Additional arguments (unused).

plot.cv_rbrm <- function(x, measure = NULL, ...) {
  if (!inherits(x, "cv_rbrm")) {
    cli::cli_abort("Object must be of class 'cv_rbrm'")
  }

  # Determine which metric to plot
  if (is.null(measure)) {
    meas_key <- x$measure
    meas_name <- if (!is.null(x$measure_name)) x$measure_name else "Performance"
    mean_val <- x$performance_mean
    se_val <- x$performance_se
    # Use the primary min/1se
    l_min <- x$lambda_min
    l_1se <- x$lambda_1se
  } else {
    measure <- tolower(measure)
    if (is.null(x$cv_results) || is.null(x$cv_results[[measure]])) {
      cli::cli_warn("Metric '{measure}' not found in CV results. Using default.")
      meas_key <- x$measure
      meas_name <- x$measure_name
      mean_val <- x$performance_mean
      se_val <- x$performance_se
      l_min <- x$lambda_min
      l_1se <- x$lambda_1se
    } else {
      res <- x$cv_results[[measure]]
      mean_val <- res$mean
      se_val <- res$se
      meas_name <- get_measure_name(measure)

      # We need to find the lambda min/1se for THIS metric to plot lines correctly
      idx_min <- which.min(mean_val)
      l_min <- x$lambdas[idx_min]

      min_perf <- mean_val[idx_min]
      se_min <- se_val[idx_min]

      idx_1se <- which(mean_val <= min_perf + se_min)
      best_idx_1se <- min(idx_1se)
      l_1se <- x$lambdas[best_idx_1se]
    }
  }

  # Extract variable counts
  fit_obj <- if (!is.null(x$path)) x$path else if (!is.null(x$final_fit) && !is.null(x$final_fit$path)) x$final_fit$path else x$fit

  if (!is.null(fit_obj)) {
    # Helper to count variables excluding intercept
    count_vars <- function(coefs, has_int) {
      if (has_int && nrow(coefs) > 0) {
        # Exclude first row (Intercept)
        colSums(abs(coefs[-1, , drop = FALSE]) > 1e-10)
      } else {
        colSums(abs(coefs) > 1e-10)
      }
    }

    has_int <- isTRUE(fit_obj$intercept)
    n_vars_a <- count_vars(fit_obj$alphas, has_int)
    n_vars_b <- count_vars(fit_obj$betas, has_int)
    n_vars <- n_vars_a + n_vars_b
  } else {
    n_vars <- rep(NA, length(x$lambdas))
  }

  plot_data <- data.frame(
    lambda = x$lambdas,
    mean_perf = mean_val,
    upper = mean_val + se_val,
    lower = mean_val - se_val,
    n_vars = round(n_vars)
  )

  # Prepare axis labels for n_vars
  n_total <- nrow(plot_data)
  n_labels <- min(n_total, 10)
  label_indices <- round(seq(1, n_total, length.out = n_labels))
  axis_labels_data <- plot_data[label_indices, ]

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = lambda, y = mean_perf)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
      width = 0.05,
      color = "grey80"
    ) +
    ggplot2::geom_point(
      color = "#f94144",
      alpha = 1
    ) +
    ggplot2::geom_vline(ggplot2::aes(
      xintercept = l_min,
      color = "Lambda.min",
      linetype = "Lambda.min"
    )) +
    ggplot2::geom_vline(ggplot2::aes(
      xintercept = l_1se,
      color = "Lambda.1se",
      linetype = "Lambda.1se"
    )) +
    ggplot2::scale_color_manual(
      name = NULL,
      values = c("Lambda.min" = "#277DA1", "Lambda.1se" = "#264653")
    ) +
    ggplot2::scale_linetype_manual(
      name = NULL,
      values = c("Lambda.min" = "dashed", "Lambda.1se" = "dotted")
    ) +
    ggplot2::labs(
      x = "Lambda",
      y = meas_name,
      title = paste("CV Performance:", meas_name),
      color = "",
      linetype = ""
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")

  if (!all(is.na(n_vars))) {
    p <- p + ggplot2::scale_x_log10(
      sec.axis = ggplot2::sec_axis(
        trans = ~.,
        name = "Number of Selected Variables",
        breaks = axis_labels_data$lambda,
        labels = axis_labels_data$n_vars
      )
    )
  } else {
    p <- p + ggplot2::scale_x_log10()
  }
#' @param x An rbrm object.
#' @param measure Whether to plot intercept (default FALSE).
#' @param ... Additional arguments (unused).


  return(p)
}

#' @param x An rbrm object.
#' @param measure Whether to plot intercept (default FALSE).
#' @param ... Additional arguments (unused).

#' Plot RBRM Coefficients (Single Fit)
#'
#' @param x An object of class `rbrm`.
#' @param measure Logical. Include intercept in plot?
#' @export
plot.rbrm <- function(x, intercept = FALSE, ...) {
  if (!inherits(x, "rbrm")) cli::cli_abort("Object must be of class 'rbrm'")

  # helper
  prep_df <- function(coefs, type) {
    df <- data.frame(
      variable = names(coefs),
      value = as.numeric(coefs),
      type = type,
      stringsAsFactors = FALSE
    )
    if (is.null(df$variable)) df$variable <- paste0(substr(type, 1, 1), 1:nrow(df))
    df
  }

  df_a <- prep_df(x$alpha, "Alpha")
  df_b <- prep_df(x$beta, "Beta")
  df <- rbind(df_a, df_b)

  if (!intercept) {
    df <- df[!grepl("^(Intercept|\\(Intercept\\))$", df$variable), ]
  }

  # Filter non-zero
  df <- df[abs(df$value) > 1e-10, , drop = FALSE]

  if (nrow(df) == 0) {
    cli::cli_warn("No non-zero coefficients to plot.")
    return(invisible(NULL))
  }

  ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(variable, value), y = value, fill = type)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::facet_wrap(~type, scales = "free", ncol = 1) +
    ggplot2::coord_flip() +
    ggplot2::labs(
#' @param x An rbrm_path object.
#' @param plot_intercept Whether to plot intercept (default FALSE).
#' @param ... Additional arguments (unused).

      title = "RBRM Coefficients",
      x = "Variable",
      y = "Value"
    ) +
    ggplot2::theme_minimal()
}


#' @param x An rbrm_path object.
#' @param plot_intercept Whether to plot intercept (default FALSE).
#' @param ... Additional arguments (unused).

#' Plot Coefficient Paths from rbrm_path Object
#'
#' @param x An object of class `rbrm_path`.
#' @param plot_intercept Logical. Include intercept?
#' @export
plot.rbrm_path <- function(x, plot_intercept = FALSE, ...) {
  if (!inherits(x, "rbrm_path")) {
    cli::cli_abort("Object must be of class 'rbrm_path'")
  }
  if (is.null(x$alphas) || is.null(x$betas) || is.null(x$lambdas)) {
    cli::cli_abort("Invalid 'rbrm_path' object: missing components.")
  }

  df_a <- as.data.frame(x$alphas)
  df_a$variable <- rownames(x$alphas)
  if (is.null(df_a$variable)) df_a$variable <- paste0("A", 1:nrow(df_a))
  df_a$type <- "Alpha"

  df_b <- as.data.frame(x$betas)
  df_b$variable <- rownames(x$betas)
  if (is.null(df_b$variable)) df_b$variable <- paste0("B", 1:nrow(df_b))
  df_b$type <- "Beta"

  colnames(df_a)[1:length(x$lambdas)] <- as.character(x$lambdas)
  colnames(df_b)[1:length(x$lambdas)] <- as.character(x$lambdas)

  df_all <- rbind(df_a, df_b)

  plot_data <- df_all %>%
    tidyr::pivot_longer(
      cols = -c(variable, type),
      names_to = "lambda",
      values_to = "coefficient"
    ) %>%
    dplyr::mutate(lambda = as.numeric(lambda))

  if (!plot_intercept) {
    plot_data <- dplyr::filter(plot_data, !grepl("^(Intercept|\\(Intercept\\))$", variable, ignore.case = TRUE))
  }

  ggplot2::ggplot(plot_data, ggplot2::aes(x = lambda, y = coefficient, group = interaction(variable, type), color = type)) +
    ggplot2::geom_line(alpha = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Lambda",
      y = "Coefficient Value",
      title = "Coefficient Regularization Paths",
      color = "Parameter Type"
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::facet_wrap(~type, scales = "free_y", ncol = 1)
}
