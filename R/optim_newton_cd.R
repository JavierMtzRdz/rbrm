#' Newton-Coordinate Descent Optimizer for RBRM
#'
#' Implements a Proximal Newton algorithm with Cyclic Coordinate Descent updates.
#' Suitable for high-dimensional, L1-penalized problems.
#'
#' @export
optim_newton_cd <- function(alpha_start, beta_start,
                            step_size_alpha, step_size_beta, # Ignored
                            lambda,
                            intercept,
                            max_step,
                            va, vb, x, y,
                            prob_fun = getProbRR.org,
                            lambda_beta = NULL,
                            tol = 1e-5,
                            save_history = FALSE,
                            ...) {
    # Setup
    if (is.null(lambda_beta)) lambda_beta <- lambda

    pa <- length(alpha_start)
    pb <- length(beta_start)

    alpha <- alpha_start
    beta <- beta_start

    # History storage
    if (save_history) {
        alphas <- matrix(0, max_step, pa)
        betas <- matrix(0, max_step, pb)
        # Gradient history?
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
    converged <- FALSE

    # Extract clipping once
    args <- list(...)
    clip <- if ("clipping" %in% names(args)) args$clipping else 1e-10

    # Helper to compute penalized objective
    calc_obj <- function(a, b) {
        penalized_nllh(a, b, va, vb, x, y, lambda, lambda_beta = lambda_beta, intercept = intercept, prob_fun = prob_fun, clipping = clip)
    }

    # Newton Loop
    for (iter in 1:max_step) {
        obj_prev <- calc_obj(alpha, beta)

        # Compute Gradient and Diagonal Hessian
        # Analytical Gradient
        grads <- grad_nll(alpha, beta, y, x, va, vb, prob_fun, opt = "both", method = "analytical", clipping = clip)
        g_alpha <- grads$grad_alpha
        g_beta <- grads$grad_beta

        # Diagonal Hessian Approximation
        # We can use a simple numerical approximation: (g(x+h) - g(x)) / h
        # Or strict diagonal of Hessian if available.
        # For robustness, let's use a "Trust Region" style positive diagonal approx
        # or just fixed step size if Hessian is ill-conditioned.
        # To be high-dim efficient, we need a vector.
        # Let's approximate diagonal H approx by finite difference of gradients on each dimension?
        # No, that's O(p) function evaluations (expensive for huge p).
        # Better: Use variance of X * weights?
        # W_ii ~ p(1-p).
        # h_j ~ sum_i X_ij^2 * W_ii.
        # Let's compute weights W_ii.

        # Compute Weights O(n)
        theta <- as.vector(va %*% alpha)
        phi <- as.vector(vb %*% beta)
        ps <- prob_fun(theta, phi)
        p0 <- ps$p0
        p1 <- ps$p1

        # Weights for NLL (approximate upper bound weights)
        # Binary NLL Hessian weights are roughly p(1-p).
        # For RBRM, it's more complex, but p(1-p) is a decent proxy for convexity.
        # Let's use pA (predicted probs)
        pA <- rep(0, length(y))
        pA[x == 0] <- p0[x == 0]
        pA[x == 1] <- p1[x == 1]

        weights <- pmax(pA * (1 - pA), 1e-4) # Safety floor

        # Diagonal Hessian O(np)
        # h_a_j = sum (va[,j]^2 * weights)
        # h_b_j = sum (vb[,j]^2 * weights)
        h_alpha <- colSums(va^2 * weights) / length(y)
        h_beta <- colSums(vb^2 * weights) / length(y)

        # Add Ridge diagonal for stability?
        h_alpha <- h_alpha + 1e-6
        h_beta <- h_beta + 1e-6

        # Coordinate Descent (Quadratic model)
        # Minimize Q(d) = g'd + 0.5 d'H d + lambda |x+d|
        # x_new = soft(x * H - g, lambda) / H

        alpha_new <- alpha
        beta_new <- beta

        # Update active set first?
        # For now, full update cycle (or active set if p is huge).

        # Alpha Update
        z_alpha <- alpha * h_alpha - g_alpha
        lam_seq_alpha <- rep(lambda, pa)
        if (intercept) lam_seq_alpha[1] <- 0
        alpha_new <- soft_thres(z_alpha, lam_seq_alpha) / h_alpha

        # Beta Update
        z_beta <- beta * h_beta - g_beta
        lam_seq_beta <- rep(lambda_beta, pb)
        if (intercept) lam_seq_beta[1] <- 0
        beta_new <- soft_thres(z_beta, lam_seq_beta) / h_beta

        # Line Search (Backtracking)
        # Direction d
        d_alpha <- alpha_new - alpha
        d_beta <- beta_new - beta

        # Norm of change
        if (max(abs(d_alpha), abs(d_beta)) < tol) {
            if (iter > 1) {
                converged <- TRUE
                break
            }
        }

        step_ls <- 1
        accepted <- FALSE
        for (ls in 1:30) {
            a_cand <- alpha + step_ls * d_alpha
            b_cand <- beta + step_ls * d_beta
            obj_cand <- calc_obj(a_cand, b_cand)

            # Simple decrease check (Armijo sufficient decrease usually better but this suffices for convex-ish)
            if (obj_cand <= obj_prev && is.finite(obj_cand)) {
                alpha <- a_cand
                beta <- b_cand
                accepted <- TRUE
                break
            }
            step_ls <- step_ls * 0.5
        }

        if (!accepted) {
            # If step 1e-3 fails, maybe converged or stuck.
            # For CD, we should perhaps take the step anyway if it's CD?
            # But this is Proximal Newton.
            # If line search fails, stop.
            break
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
        convergence = converged,
        step = iter,
        final_nll = calc_obj(alpha, beta),
        alphas = alphas_ret,
        betas = betas_ret,
        grad_alphas = grad_alphas_ret,
        grad_betas = grad_betas_ret,
        nllh_results = nllh_ret
    ))
}
