test_that("predict.rbrm handles intercepts correctly", {
    set.seed(123)
    n <- 50
    p <- 2
    va <- matrix(rnorm(n * p), n, p)
    vb <- matrix(rnorm(n * p), n, p)
    x_trt <- rbinom(n, 1, 0.5)
    y <- rbinom(n, 1, 0.5)

    # Fit WITH intercept
    fit_int <- fit.rbrm(va, vb, x_trt, y, intercept = TRUE, max_step = 10)

    # Predict with raw matrix (p columns) -> Should work (auto-pads intercept)
    preds_int <- predict(fit_int, new_va = va, new_vb = vb, type = "link")
    expect_equal(length(preds_int$logrr), n)
    expect_equal(length(preds_int$logop), n)

    # Predict with padded matrix (p+1 columns) -> Should work
    va_pad <- cbind(1, va)
    vb_pad <- cbind(1, vb)
    preds_pad <- predict(fit_int, new_va = va_pad, new_vb = vb_pad, type = "link")
    expect_equal(preds_int$logrr, preds_pad$logrr)

    # Fit WITHOUT intercept
    fit_noint <- fit.rbrm(va, vb, x_trt, y, intercept = FALSE, max_step = 10)

    # Predict with raw matrix -> Should work
    preds_noint <- predict(fit_noint, new_va = va, new_vb = vb, type = "link")
    expect_equal(length(preds_noint$logrr), n)

    # Predict with bad dimensions -> Should error
    expect_error(predict(fit_noint, new_va = va_pad, new_vb = vb_pad), "Dimension mismatch")
})

test_that("predict.rbrm returns probabilities", {
    set.seed(456)
    n <- 20
    p <- 1
    va <- matrix(rnorm(n), n, 1)
    vb <- matrix(rnorm(n), n, 1)
    x_trt <- rbinom(n, 1, 0.5)
    y <- rbinom(n, 1, 0.5)

    fit <- fit.rbrm(va, vb, x_trt, y, intercept = TRUE, max_step = 10)

    probs <- predict(fit, new_va = va, new_vb = vb, type = "response")

    expect_true("p0" %in% names(probs))
    expect_true("p1" %in% names(probs))
    expect_true(all(probs$p0 >= 0 & probs$p0 <= 1))
    expect_true(all(probs$p1 >= 0 & probs$p1 <= 1))
})
