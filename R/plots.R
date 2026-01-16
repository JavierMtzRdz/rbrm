#' Plot Cross-Validation Results
#'
#' @return A ggplot object.
#' @export
plot.cv_rbrm <- function(x, ...) {
  if (!inherits(x, "cv_rbrm")) {
    cli::cli_abort("Object must be of class 'cv_rbrm'")
  }
  if (is.null(x$nll_mean) || is.null(x$lambdas)) {
    cli::cli_abort("Invalid 'cv_rbrm' object: missing components.")
  }

  n_vars_a <- colSums(abs(x$fit$alphas) > 1e-10)
  n_vars_b <- colSums(abs(x$fit$betas) > 1e-10)
  n_vars <- n_vars_a + n_vars_b

  plot_data <- data.frame(
    lambda = x$lambdas,
    mean_dev = x$nll_mean,
    upper = x$nll_mean + x$nll_se,
    lower = x$nll_mean - x$nll_se,
    n_vars = round(n_vars)
  )

  n_total <- nrow(plot_data)
  n_labels <- min(n_total, 10)
  label_indices <- round(seq(1, n_total, length.out = n_labels))
  axis_labels_data <- plot_data[label_indices, ]

  ggplot2::ggplot(plot_data, ggplot2::aes(x = lambda, y = mean_dev)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
      width = 0.05,
      color = "grey80"
    ) +
    ggplot2::geom_point(
      color = "#f94144",
      alpha = 1
    ) +
    ggplot2::geom_vline(ggplot2::aes(
      xintercept = x$lambda_min,
      color = "Lambda.min",
      linetype = "Lambda.min"
    )) +
    ggplot2::geom_vline(ggplot2::aes(
      xintercept = x$lambda_1se,
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
      y = "Negative Log-Likelihood",
      title = "Cross-Validation Performance",
      color = "",
      linetype = ""
    ) +
    ggplot2::scale_x_log10(
      sec.axis = ggplot2::sec_axis(
        trans = ~.,
        name = "Number of Selected Variables",
        breaks = axis_labels_data$lambda,
        labels = axis_labels_data$n_vars
      )
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")
}


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
    plot_data <- dplyr::filter(plot_data, !variable %in% c("(Intercept)", "Intercept"))
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
