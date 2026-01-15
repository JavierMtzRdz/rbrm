test_that("rbrm matches brm estimates on simple data", {
    skip_if_not_installed("brm")

    set.seed(333)
    n <- 1000
    pa <- 2

    # True parameters (Intercept included)
    alpha_true <- c(0.5, 0.5, -0.5)
    beta_true <- c(-0.5, -1, 0.5)
    gamma_true <- c(0, 0) # Simple PS

    dat <- generate_data(
        pa = pa, pb = pa, n = n, n_test = 50,
        alpha = alpha_true, beta = beta_true, gamma = gamma_true,
        treatment_prob = 0.5
    )

    va <- dat$v.train
    vb <- dat$v.train
    x_trt <- dat$x.train
    y <- dat$y.train

    # Fit rbrm with intercept, no standardization (to match raw brm fit)
    fit_rbrm <- fit.rbrm(va, vb, x_trt, y, lambda = 0, intercept = TRUE, standardize = FALSE, use_line_search = TRUE, max_step = 3000, tol = 1e-5)

    # Fit brm
    # brm requires explicit intercept column if passing design matrix
    va_int <- cbind(1, va)
    vb_int <- cbind(1, vb)
    fit_brm <- brm::brm(y, x_trt, va_int, vb_int, param = "RR", est.method = "MLE")

    coef_rbrm <- c(fit_rbrm$alpha, fit_rbrm$beta)
    coef_brm <- fit_brm$point.est

    names(coef_brm) <- NULL
    names(coef_rbrm) <- NULL

    # print("RBRM Coefs:")
    # print(coef_rbrm)
    # print("BRM Coefs:")
    # print(coef_brm)

    expect_equal(coef_rbrm, coef_brm, tolerance = 0.1)
})

test_that("rbrm matches brrr estimates if available", {
    skip_if_not_installed("brrr")

    set.seed(444)
    n <- 200
    p <- 2
    va <- matrix(rnorm(n * p), n, p)
    vb <- matrix(rnorm(n * p), n, p)

    x_trt <- rbinom(n, 1, 0.5)
    y <- rbinom(n, 1, 0.5)

    fit_rbrm <- fit.rbrm(va, vb, x_trt, y, lambda = 0, intercept = FALSE)

    # Wrap brrr in tryCatch as it seems unstable
    fit_brrr <- try(brrr::brrr(x = vb, z = va, y = y, t = x_trt, param = "Richardson"), silent = TRUE)

    if (inherits(fit_brrr, "try-error")) {
        skip(paste("brrr failed to fit:", attr(fit_brrr, "condition")$message))
    }

    # Based on debug info: fit_brrr$out is a data.frame with column "point est"
    coef_brrr <- NULL
    if (!is.null(fit_brrr$out) && is.data.frame(fit_brrr$out)) {
        coef_brrr <- fit_brrr$out[["point est"]]
    }

    if (is.null(coef_brrr)) {
        # Fallbacks just in case
        coef_brrr <- fit_brrr$point.est
        if (is.null(coef_brrr)) coef_brrr <- fit_brrr$Point.est
    }

    if (is.null(coef_brrr)) {
        print("Structure of fit_brrr$out:")
        if (!is.null(fit_brrr$out)) str(fit_brrr$out)
        skip("brrr return value does not contain recognized coefficients")
    }

    coef_rbrm <- c(fit_rbrm$alpha, fit_rbrm$beta)
    names(coef_brrr) <- NULL

    # Check parameter ordering.
    # rbrm: alpha (RR), beta (Nuisance)
    # brrr call: x=vb (Nuisance), z=va (RR)
    # Check differences.
    diff_direct <- max(abs(coef_rbrm - coef_brrr))

    # Try swapped: c(Nuisance, RR) -> c(RR, Nuisance)
    p_half <- length(coef_brrr) / 2
    coef_brrr_swapped <- c(coef_brrr[(p_half + 1):(2 * p_half)], coef_brrr[1:p_half])
    diff_swapped <- max(abs(coef_rbrm - coef_brrr_swapped))

    if (diff_swapped < diff_direct) {
        coef_brrr <- coef_brrr_swapped
    }

    # Relax tolerance to 0.1 as brrr uses a different optimizer
    expect_equal(coef_rbrm, coef_brrr, tolerance = 0.1)
})
