#' Analyze Optimization Results
#'
#' Generates diagnostic plots for the optimization process
#' (trace of parameters, gradients, and objective).
#' Requires the model to be fitted with `save_history = TRUE`.
#'
#' @param model A fitted `rbrm` model object.
#' @param true_alpha Optional vector of true alpha coefficients for comparison.
#' @param true_beta Optional vector of true beta coefficients for comparison.
#' @return A list containing data frames and ggplot objects for analysis.
#' @import dplyr tidyr ggplot2 tibble
#' @importFrom stats sd uniroot var
#' @export
analyze_optimization <- function(model, true_alpha = NULL, true_beta = NULL) {
    # Handle structure where history is nested in optimizer_details
    opt_source <- model
    if (is.null(model$alphas) && !is.null(model$optimizer_details)) {
        opt_source <- model$optimizer_details
    }

    if (is.null(opt_source$alphas) || is.null(opt_source$betas)) {
        stop("Model does not contain optimization history. Refit with `save_opt = TRUE` (or `save_history=TRUE`).")
    }

    # Extract history
    # Function to tidy history matrices
    tidy_history <- function(mat, type, param_names) {
        # Ensure it's a matrix (not vector from drop=TRUE somewhere)
        if (!is.matrix(mat)) mat <- matrix(mat, nrow = 1)
        colnames(mat) <- param_names
        tibble::as_tibble(mat) %>%
            dplyr::mutate(step = dplyr::row_number()) %>%
            tidyr::pivot_longer(-step, names_to = "name", values_to = "value") %>%
            dplyr::mutate(type = type)
    }

    pa <- ncol(opt_source$alphas)
    pb <- ncol(opt_source$betas)
    alpha_names <- paste0("alpha_", 0:(pa - 1))
    beta_names <- paste0("beta_", 0:(pb - 1))

    # Helper to parse labels
    parse_labels <- function(strings) {
        # Convert "alpha_0" to "alpha[0]"
        new_strings <- gsub("_([0-9]+)", "[\\1]", strings)
        parse(text = new_strings)
    }

    # Parameter Trace
    history_alpha <- tidy_history(opt_source$alphas, "alpha", alpha_names)
    history_beta <- tidy_history(opt_source$betas, "beta", beta_names)
    history_params <- dplyr::bind_rows(history_alpha, history_beta)

    # Add ground truth if available
    if (!is.null(true_alpha) || !is.null(true_beta)) {
        truth_df <- dplyr::bind_rows(
            if (!is.null(true_alpha)) tibble::tibble(name = alpha_names, true_val = true_alpha, type = "alpha"),
            if (!is.null(true_beta)) tibble::tibble(name = beta_names, true_val = true_beta, type = "beta")
        )
        history_params <- dplyr::left_join(history_params, truth_df, by = c("name", "type"))
    }

    # Create plot_params with conditional color mapping
    if ("true_val" %in% names(history_params)) {
        # Color by true value when available
        plot_params <- ggplot(history_params, aes(step, value, color = as.factor(true_val), group = name)) +
            geom_path() +
            geom_hline(aes(yintercept = true_val, color = as.factor(true_val)), linetype = "dashed") +
            facet_wrap(~type, scales = "free_y", ncol = 1, labeller = as_labeller(c(alpha = "alpha", beta = "beta"))) +
            labs(
                title = expression(paste("Parameter Evolution ", theta)),
                y = "Value",
                x = "Iteration",
                color = expression(theta^"*")
            ) +
            theme_minimal()
    } else {
        # Default: color by parameter name
        plot_params <- ggplot(history_params, aes(step, value, color = name, group = name)) +
            geom_path() +
            facet_wrap(~type, scales = "free_y", ncol = 1, labeller = as_labeller(c(alpha = expression(alpha), beta = expression(beta)), default = label_parsed)) +
            labs(
                title = expression(paste("Parameter Evolution ", theta)),
                y = "Value",
                x = "Iteration",
                color = "Parameter"
            ) +
            theme_minimal() +
            scale_color_discrete(labels = function(x) parse_labels(x))
    }

    # Gradient Trace
    has_grads <- !is.null(opt_source$grad_alphas) && !is.null(opt_source$grad_betas)
    plot_grads <- NULL
    plot_grad_norm <- NULL

    if (has_grads) {
        grad_alpha <- tidy_history(opt_source$grad_alphas, "alpha", alpha_names)
        grad_beta <- tidy_history(opt_source$grad_betas, "beta", beta_names)
        history_grads <- dplyr::bind_rows(grad_alpha, grad_beta)

        plot_grads <- ggplot(history_grads, aes(step, value, color = name, group = name)) +
            geom_path(alpha = 0.5) +
            facet_wrap(~type, scales = "free_y", ncol = 1, labeller = as_labeller(c(alpha = "alpha", beta = "beta"))) +
            labs(title = expression("Gradient Evolution " ~ nabla ~ L), y = "Gradient", x = "Iteration") +
            theme_minimal() +
            theme(legend.position = "right") +
            scale_color_discrete(labels = function(x) parse_labels(x))

        # Gradient Norms
        grad_norms <- history_grads %>%
            dplyr::group_by(step, type) %>%
            dplyr::summarise(norm_l2 = sqrt(sum(value^2)), .groups = "drop")

        plot_grad_norm <- ggplot(grad_norms, aes(step, norm_l2, color = type)) +
            geom_line() +
            scale_y_log10() +
            labs(title = expression("Gradient Norm " ~ "||" ~ nabla ~ L ~ "||"[2]), y = expression(log[10] ~ "||" ~ nabla ~ L ~ "||"[2]), x = "Iteration") +
            theme_minimal() +
            scale_color_discrete(labels = c(alpha = expression(nabla[alpha]), beta = expression(nabla[beta])))
    }

    # Objective Trace & Optimality Gap
    has_nll <- !is.null(opt_source$nllh_results)
    plot_nll <- NULL
    plot_opt_gap <- NULL

    if (has_nll) {
        nll_vals <- opt_source$nllh_results
        if (all(nll_vals == 0)) warning("NLL History is all zeros. Check if optimization recorded NLL correctly.")

        nll_df <- tibble::tibble(step = 1:length(nll_vals), value = nll_vals)

        # Calculates "Optimality Gap" assuming min(history) is close to f*
        min_nll <- min(nll_vals[nll_vals != 0], na.rm = TRUE)
        # Compute gap and clip at small positive value to handle numerical precision
        # (sometimes NLL can be slightly below min due to rounding)
        nll_df <- nll_df %>% dplyr::mutate(gap = pmax(value - min_nll, 1e-16))

        plot_nll <- ggplot(nll_df, aes(step, value)) +
            geom_line() +
            labs(title = "Penalized Negative Log-Likelihood", y = expression(L(theta)), x = "Iteration") +
            theme_minimal()

        # Optimality Gap with Rates
        # Add reference lines: 1/k and 1/k^2
        # Match scale intercept at step 1 or 10
        ref_start <- nll_df$gap[1]
        if (is.na(ref_start) || ref_start <= 0) ref_start <- max(nll_df$gap)

        nll_df <- nll_df %>%
            dplyr::mutate(
                ref_1k = ref_start / step,
                ref_1k2 = ref_start / (step^2)
            )

        plot_opt_gap <- ggplot(nll_df, aes(step, gap)) +
            geom_line(aes(color = "Gap"), linewidth = 1) +
            geom_line(aes(y = ref_1k, color = "O(1/k)"), linetype = "dashed", alpha = 0.5) +
            geom_line(aes(y = ref_1k2, color = "O(1/k^2)"), linetype = "dotted", alpha = 0.5) +
            scale_y_log10() +
            # scale_x_log10() +
            labs(
                title = expression("Optimality Gap " ~ f(x[k]) - f^"*"),
                y = expression(log[10](Gap)),
                x = "Iteration",
                color = "Rate"
            ) +
            theme_minimal()
    }

    # Error from Truth (Norm)
    plot_error_norm <- NULL
    if (!is.null(true_alpha) && !is.null(true_beta)) {
        error_norms <- history_params %>%
            dplyr::filter(!is.na(true_val)) %>%
            dplyr::group_by(step, type) %>%
            dplyr::summarise(error_l2 = sqrt(sum((value - true_val)^2)), .groups = "drop")

        plot_error_norm <- ggplot(error_norms, aes(step, error_l2, color = type)) +
            geom_line() +
            scale_y_log10() +
            labs(
                title = expression("Parameter Error " ~ "||" ~ hat(theta) - theta^"*" ~ "||"[2]),
                y = expression(log[10] ~ Error),
                x = "Iteration"
            ) +
            theme_minimal() +
            scale_color_discrete(labels = c(alpha = expression(alpha), beta = expression(beta)))
    }

    list(
        data_params = history_params,
        plot_params = plot_params,
        plot_grads = plot_grads,
        plot_grad_norm = plot_grad_norm,
        plot_nll = plot_nll,
        plot_opt_gap = plot_opt_gap,
        plot_error_norm = plot_error_norm
    )
}

utils::globalVariables(c("step", "value", "name", "type", "true_val", "norm_l2", "gap", "ref_1k", "ref_1k2", "error_l2", "coefficient", "lambda", "lower", "mean_perf", "upper", "variable"))
