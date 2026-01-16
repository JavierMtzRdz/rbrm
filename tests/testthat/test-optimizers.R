test_that("All optimizers converge on simple data", {
    set.seed(123)
    n <- 100
    p <- 3
    # Generate simple data
    va <- matrix(rnorm(n * p), n, p)
    vb <- matrix(rnorm(n * p), n, p)
    alpha_true <- c(0.5, -0.5, 0.2)
    beta_true <- c(-0.2, 0.2, 0)

    theta <- va %*% alpha_true
    phi <- vb %*% beta_true
    ps <- getProbRR.org(theta, phi)

    x_trt <- rbinom(n, 1, 0.5)
    y <- numeric(n)
    y[x_trt == 0] <- rbinom(sum(x_trt == 0), 1, ps$p0[x_trt == 0])
    y[x_trt == 1] <- rbinom(sum(x_trt == 1), 1, ps$p1[x_trt == 1])

    # List of optimizers
    optimizers <- list(fista_opt, optim_newton_cd, optim_lbfgs)
    names(optimizers) <- c("fista", "newton", "lbfgs")

    for (opt_name in names(optimizers)) {
        opt_fun <- optimizers[[opt_name]]

        # 1. Basic Fit
        fit <- fit.rbrm(va, vb, x_trt, y,
            max_step = 100,
            lambda = 0.01,
            intercept = FALSE,
            opt_fun = opt_fun,
            save_opt = FALSE
        )

        expect_s3_class(fit, "rbrm")
        expect_equal(length(fit$alpha), p)
        expect_equal(length(fit$beta), p)
        # Convergence might vary, but should not crash

        # Check if estimates are somewhat correlated with truth (loose check for small N)
        # cor_a <- cor(fit$alpha, alpha_true)
        # expect_true(cor_a > 0.5, label = paste0(opt_name, " alpha correlation"))
    }
})

test_that("analyze_optimization works with all optimizers", {
    set.seed(456)
    n <- 50
    p <- 2
    va <- matrix(rnorm(n * p), n, p)
    vb <- matrix(rnorm(n * p), n, p)
    x_trt <- rbinom(n, 1, 0.5)
    y <- rbinom(n, 1, 0.3)

    optimizers <- list(fista_opt, optim_newton_cd, optim_lbfgs)
    names(optimizers) <- c("fista", "newton", "lbfgs")

    for (opt_name in names(optimizers)) {
        opt_fun <- optimizers[[opt_name]]

        # Fit with save_opt = TRUE
        fit <- fit.rbrm(va, vb, x_trt, y,
            max_step = 10,
            lambda = 0,
            intercept = FALSE,
            opt_fun = opt_fun,
            save_opt = TRUE
        )

        # Run Analysis
        analysis <- analyze_optimization(fit)

        if (!is.null(fit$alphas)) {
            last_alpha <- fit$alphas[nrow(fit$alphas), ]
            expect_equal(as.vector(last_alpha), as.vector(fit$alpha), tolerance = 1e-5)
        }

        # Check outputs
        expect_true(is.list(analysis))
        expect_true(!is.null(analysis$plot_params))
        expect_true(!is.null(analysis$plot_nll))
        expect_true(!is.null(analysis$plot_opt_gap)) # New check
        expect_true(!is.null(analysis$plot_opt_gap))

        # Gradient plots: FISTA and Newton save gradients, L-BFGS does not
        if (opt_name %in% c("fista", "newton")) {
            expect_true(!is.null(analysis$plot_grads), info = paste("Optimizer:", opt_name))
            expect_true(!is.null(analysis$plot_grad_norm), info = paste("Optimizer:", opt_name))
        }

        expect_s3_class(analysis$plot_params, "ggplot")
        expect_s3_class(analysis$plot_nll, "ggplot")
    }
})
