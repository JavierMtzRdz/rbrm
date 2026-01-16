#' Active Set Newton-Coordinate Descent Optimizer for RBRM
#'
#' Implements a Proximal Newton algorithm with Active Set strategy and Cyclic Coordinate Descent.
#' More efficient than full Newton-CD for sparse problems (L1-penalized).
#'
#' @param kkt_check_freq How often to check KKT conditions for inactive variables (default: 10 iterations)
#' @param active_tol Tolerance for determining if a coefficient is active (default: 1e-6)
#' @inheritParams optim_newton_cd
#' @export
optim_newton_cd_active <- function(alpha_start, beta_start,
                                   step_size_alpha, step_size_beta,
                                   lambda,
                                   intercept,
                                   max_step,
                                   va, vb, x, y,
                                   prob_fun = getProbRR.org,
                                   lambda_beta = NULL,
                                   tol = 1e-5,
                                   save_history = FALSE,
                                   kkt_check_freq = 10,
                                   active_tol = 1e-6,
                                   ...) {
    # Setup
    if (is.null(lambda_beta)) lambda_beta <- lambda

    pa <- length(alpha_start)
    pb <- length(beta_start)

    alpha <- alpha_start
    beta <- beta_start

    # Extract clipping once
    args <- list(...)
    clip <- if ("clipping" %in% names(args)) args$clipping else 1e-10

    # Helper to compute penalized objective
    calc_obj <- function(a, b) {
        penalized_nllh(a, b, va, vb, x, y, lambda, lambda_beta = lambda_beta, intercept = intercept, prob_fun = prob_fun, clipping = clip)
    }

    # Initialize active sets
    # Start with all coefficients active
    active_alpha <- rep(TRUE, pa)
    active_beta <- rep(TRUE, pb)

    # Intercept is always active if present
    if (intercept) {
        active_alpha[1] <- TRUE
        active_beta[1] <- TRUE
    }

    # History storage
    if (save_history) {
        alphas <- matrix(0, max_step, pa)
        betas <- matrix(0, max_step, pb)
        grad_alphas <- matrix(0, max_step, pa)
        grad_betas <- matrix(0, max_step, pb)
        nllh_results <- vector("numeric", max_step)
    } else {
        alphas <- NULL
        betas <- NULL
        grad_alphas <- NULL
        grad_betas <- NULL
        nllh_results <- NULL
    }
    saved_steps <- 0

    # Newton Loop
    for (iter in 1:max_step) {
        obj_prev <- calc_obj(alpha, beta)

        # Compute Gradient and Diagonal Hessian
        grads <- grad_nll(alpha, beta, y, x, va, vb, prob_fun, opt = "both", method = "analytical", clipping = clip)
        g_alpha <- grads$grad_alpha
        g_beta <- grads$grad_beta

        # Compute Weights
        theta <- as.vector(va %*% alpha)
        phi <- as.vector(vb %*% beta)
        ps <- prob_fun(theta, phi)
        p0 <- ps$p0
        p1 <- ps$p1

        pA <- rep(0, length(y))
        pA[x == 0] <- p0[x == 0]
        pA[x == 1] <- p1[x == 1]

        weights <- pmax(pA * (1 - pA), 1e-4)

        # Diagonal Hessian
        h_alpha <- colSums(va^2 * weights) / length(y)
        h_beta <- colSums(vb^2 * weights) / length(y)

        h_alpha <- h_alpha + 1e-6
        h_beta <- h_beta + 1e-6

        # KKT Check for Inactive Variables
        if (iter %% kkt_check_freq == 0 || iter == 1) {
            if (lambda > 0 && sum(!active_alpha) > 0) {
                inactive_idx <- which(!active_alpha)
                if (intercept && 1 %in% inactive_idx) {
                    inactive_idx <- setdiff(inactive_idx, 1)
                }

                if (length(inactive_idx) > 0) {
                    # Compute gradient at zero for inactive variables
                    violated <- abs(z_alpha[inactive_idx]) > lambda + active_tol

                    if (any(violated)) {
                        active_alpha[inactive_idx[violated]] <- TRUE
                    }
                }
            }

            # Beta
            if (lambda_beta > 0 && sum(!active_beta) > 0) {
                inactive_idx <- which(!active_beta)
                if (intercept && 1 %in% inactive_idx) {
                    inactive_idx <- setdiff(inactive_idx, 1)
                }

                if (length(inactive_idx) > 0) {
                    z_beta <- beta * h_beta - g_beta
                    violated <- abs(z_beta[inactive_idx]) > lambda_beta + active_tol

                    if (any(violated)) {
                        active_beta[inactive_idx[violated]] <- TRUE
                    }
                }
            }
        }

        # Coordinate Descent on Active Set Only
        alpha_new <- alpha
        beta_new <- beta

        # Alpha Update
        active_idx_a <- which(active_alpha)
        if (length(active_idx_a) > 0) {
            z_alpha <- alpha * h_alpha - g_alpha
            lam_seq_alpha <- rep(lambda, pa)
            if (intercept) lam_seq_alpha[1] <- 0

            # Update only active coefficients
            alpha_new[active_idx_a] <- soft_thres(z_alpha[active_idx_a], lam_seq_alpha[active_idx_a]) / h_alpha[active_idx_a]
        }

        # Beta Update
        active_idx_b <- which(active_beta)
        if (length(active_idx_b) > 0) {
            z_beta <- beta * h_beta - g_beta
            lam_seq_beta <- rep(lambda_beta, pb)
            if (intercept) lam_seq_beta[1] <- 0

            beta_new[active_idx_b] <- soft_thres(z_beta[active_idx_b], lam_seq_beta[active_idx_b]) / h_beta[active_idx_b]
        }

        # Shrink inactive variables to exactly zero
        alpha_new[!active_alpha] <- 0
        beta_new[!active_beta] <- 0

        # Line Search (Backtracking)
        d_alpha <- alpha_new - alpha
        d_beta <- beta_new - beta

        # Norm of change
        if (sum(d_alpha^2) + sum(d_beta^2) < tol^2) {
            if (iter > 1) break
        }

        step_ls <- 1
        accepted <- FALSE
        for (ls in 1:30) {
            a_cand <- alpha + step_ls * d_alpha
            b_cand <- beta + step_ls * d_beta
            obj_cand <- calc_obj(a_cand, b_cand)

            if (obj_cand <= obj_prev) {
                alpha <- a_cand
                beta <- b_cand
                accepted <- TRUE
                break
            }
            step_ls <- step_ls * 0.5
        }

        if (!accepted) {
            break
        }

        # Update active sets based on current solution
        if (lambda > 0) {
            active_alpha <- abs(alpha) > active_tol
            if (intercept) active_alpha[1] <- TRUE
        }
        if (lambda_beta > 0) {
            active_beta <- abs(beta) > active_tol
            if (intercept) active_beta[1] <- TRUE
        }

        if (save_history) {
            alphas[iter, ] <- alpha
            betas[iter, ] <- beta
            grad_alphas[iter, ] <- g_alpha
            grad_betas[iter, ] <- g_beta
            nllh_results[iter] <- obj_cand
            saved_steps <- iter
        }
    }

    alphas_ret <- if (save_history && saved_steps > 0) alphas[1:saved_steps, , drop = FALSE] else matrix(alpha, nrow = 1)
    betas_ret <- if (save_history && saved_steps > 0) betas[1:saved_steps, , drop = FALSE] else matrix(beta, nrow = 1)
    grad_alphas_ret <- if (save_history && saved_steps > 0) grad_alphas[1:saved_steps, , drop = FALSE] else NULL
    grad_betas_ret <- if (save_history && saved_steps > 0) grad_betas[1:saved_steps, , drop = FALSE] else NULL
    nllh_ret <- if (save_history && saved_steps > 0) nllh_results[1:saved_steps] else NULL

    return(list(
        alpha = alpha,
        beta = beta,
        convergence = (iter < max_step),
        step = iter,
        final_nll = calc_obj(alpha, beta),
        alphas = alphas_ret,
        betas = betas_ret,
        grad_alphas = grad_alphas_ret,
        grad_betas = grad_betas_ret,
        nllh_results = nllh_ret,
        n_active_alpha = sum(abs(alpha) > active_tol),
        n_active_beta = sum(abs(beta) > active_tol)
    ))
}
