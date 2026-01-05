test_that("Analytical gradients (Richardson) match numerical gradients", {
    set.seed(123)
    n <- 50
    p <- 3

    # Generate synthetic data
    data <- rbrm::generate_data(pa = p, pb = p, n = n, n_test = 10, alpha = rep(0, p), beta = rep(0, p), gamma = rep(0, p))
    va <- data$v.train
    vb <- data$v.train
    x <- data$x.train
    y <- data$y.train

    # Random params
    alpha_eval <- rnorm(p)
    beta_eval <- rnorm(p)

    # 1. Test Alpha Gradient
    nll_alpha <- function(a) {
        rbrm::nllh(alpha = a, beta = beta_eval, va = va, vb = vb, x = x, y = y, prob_fun = rbrm::getProbRR.org)
    }

    num_grad_alpha <- numDeriv::grad(nll_alpha, alpha_eval)

    # Analytical Gradient (via grad_nll which uses dp0_theta/dp0_phi)
    ana_grad_full <- rbrm::grad_nll(alpha = alpha_eval, beta = beta_eval, va = va, vb = vb, x = x, y = y, prob_fun = rbrm::getProbRR.org)
    ana_grad_alpha <- as.vector(ana_grad_full$grad_alpha)

    expect_equal(ana_grad_alpha, num_grad_alpha, tolerance = 1e-4, label = "Alpha Gradient (Richardson)")

    # 2. Test Beta Gradient
    nll_beta <- function(b) {
        rbrm::nllh(alpha = alpha_eval, beta = b, va = va, vb = vb, x = x, y = y, prob_fun = rbrm::getProbRR.org)
    }

    num_grad_beta <- numDeriv::grad(nll_beta, beta_eval)
    ana_grad_beta <- as.vector(ana_grad_full$grad_beta)

    expect_equal(ana_grad_beta, num_grad_beta, tolerance = 1e-4, label = "Beta Gradient (Richardson)")
})

test_that("Analytical gradients (Pozza) match numerical gradients", {
    set.seed(456)
    n <- 50
    p <- 3

    data <- rbrm::generate_data(pa = p, pb = p, n = n, n_test = 10, alpha = rep(0, p), beta = rep(0, p), gamma = rep(0, p))
    va <- data$v.train
    vb <- data$v.train
    x <- data$x.train
    y <- data$y.train

    alpha_eval <- rnorm(p)
    beta_eval <- rnorm(p)

    # 1. Test Alpha Gradient
    nll_alpha <- function(a) {
        rbrm::nllh(alpha = a, beta = beta_eval, va = va, vb = vb, x = x, y = y, prob_fun = rbrm::getProbRR.alt)
    }

    num_grad_alpha <- numDeriv::grad(nll_alpha, alpha_eval)
    ana_grad_full <- rbrm::grad_nll(alpha = alpha_eval, beta = beta_eval, va = va, vb = vb, x = x, y = y, prob_fun = rbrm::getProbRR.alt)
    ana_grad_alpha <- as.vector(ana_grad_full$grad_alpha)

    expect_equal(ana_grad_alpha, num_grad_alpha, tolerance = 1e-4, label = "Alpha Gradient (Pozza)")

    # 2. Test Beta Gradient
    nll_beta <- function(b) {
        rbrm::nllh(alpha = alpha_eval, beta = b, va = va, vb = vb, x = x, y = y, prob_fun = rbrm::getProbRR.alt)
    }

    num_grad_beta <- numDeriv::grad(nll_beta, beta_eval)
    ana_grad_beta <- as.vector(ana_grad_full$grad_beta)

    expect_equal(ana_grad_beta, num_grad_beta, tolerance = 1e-4, label = "Beta Gradient (Pozza)")
})

test_that("Gradients handle near-zero phi correctly", {
    # This tests the stability logic for phi ~ 0
    n <- 10
    p <- 1
    va <- matrix(rnorm(n), ncol = 1)
    vb <- matrix(rnorm(n), ncol = 1)
    x <- rbinom(n, 1, 0.5)
    y <- rbinom(n, 1, 0.5)

    alpha <- 0.5
    # Choose beta such that vb %*% beta is close to 0 for some observations
    # Here we just force beta to be small
    beta <- 1e-8

    nll_beta <- function(b) {
        rbrm::nllh(alpha = alpha, beta = b, va = va, vb = vb, x = x, y = y, prob_fun = rbrm::getProbRR.org)
    }

    num_grad_beta <- as.vector(numDeriv::grad(nll_beta, beta))
    ana_grad_full <- rbrm::grad_nll(alpha = alpha, beta = beta, va = va, vb = vb, x = x, y = y, prob_fun = rbrm::getProbRR.org)
    ana_grad_beta <- as.vector(ana_grad_full$grad_beta)

    expect_equal(ana_grad_beta, num_grad_beta, tolerance = 1e-2, label = "Beta Gradient (Small Phi)")
})

test_that("Gradients handle large phi boundary correctly", {
    n <- 10
    p <- 1
    va <- matrix(rnorm(n), ncol = 1)
    vb <- matrix(rep(1, n), ncol = 1) # All 1s
    x <- rbinom(n, 1, 0.5)
    y <- rbinom(n, 1, 0.5)

    alpha <- 0.5
    # Large beta to trigger phi > 12
    beta <- 15

    nll_beta <- function(b) {
        rbrm::nllh(alpha = alpha, beta = b, va = va, vb = vb, x = x, y = y, prob_fun = rbrm::getProbRR.org)
    }

    num_grad_beta <- as.vector(numDeriv::grad(nll_beta, beta))
    ana_grad_full <- rbrm::grad_nll(alpha = alpha, beta = beta, va = va, vb = vb, x = x, y = y, prob_fun = rbrm::getProbRR.org)
    ana_grad_beta <- as.vector(ana_grad_full$grad_beta)

    expect_equal(ana_grad_beta, num_grad_beta, tolerance = 1e-3, label = "Beta Gradient (Large Phi)")
})
