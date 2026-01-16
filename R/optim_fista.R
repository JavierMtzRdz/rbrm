#' Internal helper for FISTA update step (Line Search or BB)
#' @keywords internal
perform_fista_step <- function(param, y_param, grad, step_size, lambda, intercept,
                               cost_fun, current_cost,
                               use_line_search, armijo_c, line_search_shrink, line_search_max_iter,
                               last_grad = NULL, last_param = NULL, iter = 1) {
    p <- length(param)
    param_next <- param # default
    eta <- step_size

    if (use_line_search) {
        for (ls_iter in 1:line_search_max_iter) {
            z_trial <- y_param - eta * grad
            param_trial <- soft_thres(z_trial, lambda * eta)
            if (intercept && p > 0) param_trial[1] <- z_trial[1]

            f_trial <- cost_fun(param_trial)

            # Armijo condition roughly: f(new) <= f(current) - c * eta * ||grad||^2
            # We use f(current) as approximation for f(y) upper bound or similar
            armijo_rhs <- current_cost - armijo_c * eta * sum(grad^2)

            if (!is.na(f_trial) && (f_trial <= armijo_rhs || eta < 1e-12)) {
                param_next <- param_trial
                break
            }
            eta <- eta * line_search_shrink
        }
        # Return updated step size? Or keep fixed? function usually adapts locally.
        # Return eta as next step?
        return(list(param_next = param_next, step_size_next = eta))
    } else {
        # Barzilai-Borwein
        bb_step <- step_size
        if (!is.null(last_grad) && iter > 1) {
            s_param <- param - last_param
            y_grad <- grad - last_grad
            sy <- sum(s_param * y_grad)
            if (sy > 1e-8) {
                bb_step <- sum(s_param * s_param) / sy
                # Bound the step size
                bb_step <- max(min(bb_step, step_size * 2), step_size * 0.1)
            }
        }
        z_param <- y_param - bb_step * grad
        param_next <- soft_thres(z_param, lambda * bb_step)
        if (intercept && p > 0) param_next[1] <- z_param[1]

        return(list(param_next = param_next, step_size_next = bb_step))
    }
}

#' Optimized FISTA Algorithm
#'
#' Performs FISTA optimization for RBRM with inlined updates and efficiency improvements.
#'
#' @param alpha_start Initial alpha vector.
#' @param beta_start Initial beta vector.
#' @param step_size_alpha Fixed step size for alpha.
#' @param step_size_beta Fixed step size for beta.
#' @param lambda L1 penalty.
#' @param intercept Logical, whether to exclude first coef from penalty.
#' @param max_step Maximum iterations.
#' @param va Covariate matrix for alpha.
#' @param vb Covariate matrix for beta.
#' @param x Treatment vector.
#' @param y Outcome vector.
#' @param prob_fun Probability function.
#' @param lambda_beta Optional separate lambda for beta.
#' @param eval_grad Boolean, whether to evaluate gradients for stopping.
#' @param save_history Boolean, whether to store full parameter history (memory intensive).
#' @param thres Convergence threshold.
#' @param use_line_search Boolean, use backtracking line search?
#' @param armijo_c Armijo condition constant.
#' @param line_search_shrink Shrinkage factor for line search.
#' @param line_search_max_iter Maximum line search iterations.
#' @param clipping Probability clipping threshold.
#' @param ... Additional arguments.
#'
#' @export
optim_fista <- function(alpha_start, beta_start,
                        step_size_alpha = NULL, step_size_beta = NULL,
                        lambda,
                        intercept,
                        max_step,
                        va, vb, x, y,
                        prob_fun = getProbRR.org,
                        lambda_beta = NULL,
                        eval_grad = TRUE,
                        save_history = FALSE,
                        thres = 1e-8,
                        use_line_search = TRUE, # Enabled by default for robustness
                        armijo_c = 1e-4,
                        line_search_shrink = 0.5,
                        line_search_max_iter = 20,
                        clipping = 1e-10,
                        ...) {
    if (is.null(lambda_beta)) lambda_beta <- lambda

    p_a <- length(alpha_start)
    p_b <- length(beta_start)

    # Adaptive step size initialization for high dimensions
    if (is.null(step_size_alpha)) {
        # Conservative: 0.1 for high-dim, 0.5 for low-dim
        step_size_alpha <- ifelse(p_a > 20, 0.1, 0.5)
    }
    if (is.null(step_size_beta)) {
        step_size_beta <- ifelse(p_b > 20, 0.1, 0.5)
    }
    # p_a, p_b already defined above

    alpha <- alpha_start
    beta <- beta_start

    # Momentum variables
    last_alpha <- alpha
    last_beta <- beta

    t_cur <- 1

    alpha_best <- alpha
    beta_best <- beta
    f_best <- Inf

    # Adaptive restart tracking
    restart_count <- 0
    f_prev <- Inf

    # Barzilai-Borwein step size tracking
    last_g_alpha <- NULL
    last_g_beta <- NULL

    current_step_alpha <- step_size_alpha
    current_step_beta <- step_size_beta

    if (save_history) {
        alphas <- matrix(0, max_step, p_a)
        betas <- matrix(0, max_step, p_b)
        g_alphas <- matrix(0, max_step, p_a)
        g_betas <- matrix(0, max_step, p_b)
        nllh_results <- vector("double", max_step)
    } else {
        alphas <- NULL
        betas <- NULL
        g_alphas <- NULL
        g_betas <- NULL
        nllh_results <- NULL
    }

    # Constants for relative change check
    stabilizer <- sqrt(.Machine$double.eps)

    # Step Loop
    step <- 0
    convergence <- FALSE

    for (iter in 1:max_step) {
        step <- iter

        # Adaptive Restart & Monotonicity Check
        f_current <- penalized_nllh(alpha, beta, va, vb, x, y, lambda, lambda_beta = lambda_beta, lambda_b_prop = 1, intercept = intercept, prob_fun = prob_fun, clipping = clipping)

        # Monotonicity / Restart Check
        do_restart <- FALSE
        if (iter > 1 && f_current > f_prev) {
            do_restart <- TRUE
        }

        if (do_restart) {
            restart_count <- restart_count + 1
            # Revert to previous
            alpha <- last_alpha
            beta <- last_beta
            f_current <- f_prev

            # Reset momentum
            t_cur <- 1
        }

        # Update best solution tracking
        if (f_current < f_best) {
            alpha_best <- alpha
            beta_best <- beta
            f_best <- f_current
        }
        f_prev <- f_current

        # Momentum Step
        t_next <- (1 + sqrt(1 + 4 * t_cur^2)) / 2
        momentum <- (t_cur - 1) / t_next

        y_alpha <- alpha + momentum * (alpha - last_alpha)

        # Gradient calculation at y_alpha
        g_alpha <- grad_nll(y_alpha, beta, y, x, va, vb, prob_fun, opt = "alpha", method = "analytical", clipping = clipping)

        # Check Gradient Restart condition (y - x) * g > 0 (Angle between momentum and gradient step)
        if (sum((y_alpha - alpha) * g_alpha) > 0) {
            t_cur <- 1
            t_next <- (1 + sqrt(1 + 4 * t_cur^2)) / 2
            momentum <- 0
            y_alpha <- alpha # Reset y to x
            g_alpha <- grad_nll(y_alpha, beta, y, x, va, vb, prob_fun, opt = "alpha", method = "analytical", clipping = clipping)
        }

        # Alpha Update
        # Closure for cost function
        cost_alpha <- function(a) {
            penalized_nllh(a, beta, va, vb, x, y, lambda, lambda_beta = lambda_beta, lambda_b_prop = 1, intercept = intercept, prob_fun = prob_fun, clipping = clipping)
        }

        res_alpha <- perform_fista_step(
            param = alpha, y_param = y_alpha, grad = g_alpha,
            step_size = current_step_alpha, lambda = lambda, intercept = intercept,
            cost_fun = cost_alpha, current_cost = f_current,
            use_line_search = use_line_search, armijo_c = armijo_c,
            line_search_shrink = line_search_shrink, line_search_max_iter = line_search_max_iter,
            last_grad = last_g_alpha, last_param = last_alpha, iter = iter
        )
        alpha_next <- res_alpha$param_next
        if (use_line_search) {
            current_step_alpha <- res_alpha$step_size_next
        } else {
            current_step_alpha <- res_alpha$step_size_next # BB step
        }

        last_g_alpha <- g_alpha

        # Beta Update
        # Apply momentum
        y_beta <- beta + momentum * (beta - last_beta)

        # Calculate gradient at y_beta
        g_beta <- grad_nll(alpha_next, y_beta, y, x, va, vb, prob_fun, opt = "beta", method = "analytical", clipping = clipping)

        if (sum((y_beta - beta) * g_beta) > 0) {
            t_cur <- 1
            t_next <- (1 + sqrt(1 + 4 * t_cur^2)) / 2
            momentum <- 0
            y_beta <- beta
            g_beta <- grad_nll(alpha_next, y_beta, y, x, va, vb, prob_fun, opt = "beta", method = "analytical", clipping = clipping)
        }

        cost_beta <- function(b) {
            penalized_nllh(alpha_next, b, va, vb, x, y, lambda, lambda_beta = lambda_beta, lambda_b_prop = 1, intercept = intercept, prob_fun = prob_fun, clipping = clipping)
        }

        res_beta <- perform_fista_step(
            param = beta, y_param = y_beta, grad = g_beta,
            step_size = current_step_beta, lambda = lambda_beta, intercept = intercept,
            cost_fun = cost_beta, current_cost = f_current, # Using f_current as approx base.
            use_line_search = use_line_search, armijo_c = armijo_c,
            line_search_shrink = line_search_shrink, line_search_max_iter = line_search_max_iter,
            last_grad = last_g_beta, last_param = last_beta, iter = iter
        )
        beta_next <- res_beta$param_next
        if (use_line_search) {
            current_step_beta <- res_beta$step_size_next
        } else {
            current_step_beta <- res_beta$step_size_next
        }

        last_g_beta <- g_beta

        # Check Monotonicity
        f_next <- penalized_nllh(alpha_next, beta_next, va, vb, x, y, lambda, lambda_beta = lambda_beta, lambda_b_prop = 1, intercept = intercept, prob_fun = prob_fun, clipping = clipping)

        if (f_next > f_current) {
            # Restart momentum and reject step
            alpha_next <- alpha
            beta_next <- beta
            t_next <- 1 # Reset momentum for next iter
            f_next <- f_current
            restart_count <- restart_count + 1
        }

        # Check Convergence
        # Proximal Gradient Norm: ||(x - x_next) / step||

        pg_a <- sqrt(sum(((alpha - alpha_next) / current_step_alpha)^2))
        pg_b <- sqrt(sum(((beta - beta_next) / current_step_beta)^2))

        prox_conv <- (pg_a < thres && pg_b < thres)

        denom_a <- max(sum(alpha^2), sum(alpha_next^2), stabilizer)
        rel_a <- sum((alpha_next - alpha)^2) / denom_a

        denom_b <- max(sum(beta^2), sum(beta_next^2), stabilizer)
        rel_b <- sum((beta_next - beta)^2) / denom_b

        rel_conv <- (rel_a < thres && rel_b < thres)

        if (save_history) {
            alphas[iter, ] <- alpha_next
            betas[iter, ] <- beta_next
            g_alphas[iter, ] <- g_alpha
            g_betas[iter, ] <- g_beta
            nllh_results[iter] <- penalized_nllh(alpha_next, beta_next, va, vb, x, y, lambda, lambda_beta = lambda_beta, intercept = intercept, prob_fun = prob_fun, clipping = clipping)
        }

        last_alpha <- alpha
        last_beta <- beta
        alpha <- alpha_next
        beta <- beta_next
        t_cur <- t_next

        if (iter >= 2) {
            if (eval_grad) {
                if (prox_conv) {
                    convergence <- TRUE
                    break
                }
            } else {
                if (rel_conv) {
                    convergence <- TRUE
                    break
                }
            }
        }

        if (!all(is.finite(alpha)) || !all(is.finite(beta))) {
            warning("Non-finite parameters detected in FISTA. Stopping.")
            convergence <- FALSE # Failed
            break
        }
    }

    if (save_history && step < max_step) {
        alphas <- alphas[1:step, , drop = FALSE]
        betas <- betas[1:step, , drop = FALSE]
        g_alphas <- g_alphas[1:step, , drop = FALSE]
        g_betas <- g_betas[1:step, , drop = FALSE]
        nllh_results <- nllh_results[1:step]
    }

    return(list(
        alpha = alpha_best,
        beta = beta_best,
        step = step,
        convergence = convergence,
        alphas = alphas,
        betas = betas,
        grad_alphas = g_alphas,
        grad_betas = g_betas,
        nllh_results = nllh_results,
        restart_count = restart_count,
        final_nll = f_best
    ))
}
