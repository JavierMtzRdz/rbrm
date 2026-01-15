test_that("fit.rbrm converges on simple separable data", {
    set.seed(123)
    n <- 50
    p <- 2
    va <- matrix(rnorm(n * p), n, p)
    vb <- matrix(rnorm(n * p), n, p)

    # Simple true model
    alpha_true <- c(1, -1)
    beta_true <- c(-0.5, 0.5)

    theta <- va %*% alpha_true
    phi <- vb %*% beta_true
    ps <- getProbRR.org(theta, phi)

    # Simulate data with no noise to ensure easy convergence
    # But binary data always has 'noise' due to sampling.
    # We just want reasonable separation.
    y <- rbinom(n, 1, ps$p0)
    x_trt <- rbinom(n, 1, 0.5)
    y <- ifelse(x_trt == 1, rbinom(n, 1, ps$p1), y)

    fit <- fit.rbrm(va, vb, x_trt, y, max_step = 200, lambda = 0, intercept = FALSE, thres = 1e-3)

    expect_true(fit$convergence)
    expect_equal(length(fit$alpha), p)
    expect_equal(length(fit$beta), p)
})

test_that("fit.rbrm recovers parameters in large sample", {
    skip_on_cran()
    set.seed(42)
    n <- 3000 # Large n for recovery
    pa <- 3

    # True parameters (including intercept)
    # Intercept, A1, A2, A3
    alpha_true <- c(0.2, 0.5, -0.5, 0)
    beta_true <- c(-0.5, -1, 0.5, 0)
    gamma_true <- c(0.5, -0.5, 0) # For Propensity Score (no intercept in input, calculated by func)

    dat <- generate_data(
        pa = pa, pb = pa, n = n, n_test = 100,
        alpha = alpha_true, beta = beta_true, gamma = gamma_true,
        treatment_prob = 0.5
    )

    va <- dat$v.train
    vb <- dat$v.train # generate_data uses same v for both parts currently
    x_trt <- dat$x.train
    y <- dat$y.train

    # Fit with intercept=TRUE (since generated data has intercepts)
    fit <- fit.rbrm(va, vb, x_trt, y, max_step = 3000, lambda = 0, intercept = TRUE, thres = 1e-4)

    if (!fit$convergence) warning("Fit did not converge within max_step")

    # Compare coefficients
    rmse_alpha <- sqrt(mean((fit$alpha - alpha_true)^2))
    rmse_beta <- sqrt(mean((fit$beta - beta_true)^2))

    expect_lt(rmse_alpha, 0.5)
    expect_lt(rmse_beta, 0.8)
})

test_that("fit.rbrm handles standardization and intercepts correctly", {
    set.seed(888)
    n <- 100
    p <- 2
    va <- matrix(rnorm(n * p, mean = 10, sd = 2), n, p) # Non-centered
    vb <- matrix(rnorm(n * p, mean = -5, sd = 2), n, p)

    x_trt <- rbinom(n, 1, 0.5)
    y <- rbinom(n, 1, 0.3) # Random y

    # Fit with intercept and standardization (default)
    fit <- fit.rbrm(va, vb, x_trt, y, max_step = 50, intercept = TRUE, standardize = TRUE)

    # Output alpha should have length p+1 (intercept)
    expect_equal(length(fit$alpha), p + 1)
    expect_equal(length(fit$beta), p + 1)

    # Fit without standardization
    fit_ns <- fit.rbrm(va, vb, x_trt, y, max_step = 50, intercept = TRUE, standardize = FALSE)
    expect_equal(length(fit_ns$alpha), p + 1)

    # Values won't be identical due to optimization path differences on scaled vs unscaled landscape,
    # but they should be in similar magnitude range if converged.
})

test_that("L1 regularization (Lambda) shrinks coefficients", {
    set.seed(99)
    n <- 50
    p <- 5
    va <- matrix(rnorm(n * p), n, p)
    vb <- matrix(rnorm(n * p), n, p)
    x_trt <- rbinom(n, 1, 0.5)
    y <- rbinom(n, 1, 0.5)

    # Fit with no lambda
    fit_0 <- fit.rbrm(va, vb, x_trt, y, lambda = 0, max_step = 100, intercept = FALSE)
    l1_0 <- sum(abs(fit_0$alpha)) + sum(abs(fit_0$beta))

    # Fit with large lambda
    fit_High <- fit.rbrm(va, vb, x_trt, y, lambda = 0.5, max_step = 100, intercept = FALSE)
    l1_High <- sum(abs(fit_High$alpha)) + sum(abs(fit_High$beta))

    expect_lt(l1_High, l1_0)

    # Sparse result check
    expect_true(sum(fit_High$alpha == 0) + sum(fit_High$beta == 0) >
        sum(fit_0$alpha == 0) + sum(fit_0$beta == 0))
})
