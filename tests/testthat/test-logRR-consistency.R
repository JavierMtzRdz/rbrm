test_that("getProbRR.org maintains p1 = p0 * exp(logrr) relationship", {
    # Test across a grid of values
    logrr_vals <- seq(-5, 5, length.out = 20)
    logop_vals <- seq(-5, 5, length.out = 20)

    grid <- expand.grid(logrr = logrr_vals, logop = logop_vals)

    for (i in 1:nrow(grid)) {
        theta <- grid$logrr[i]
        phi <- grid$logop[i]

        probs <- rbrm::getProbRR.org(theta, phi)
        p0 <- probs$p0
        p1 <- probs$p1

        # Expected relationship: p1 = p0 * exp(theta)
        # Check within tolerance (allowing for potential floating point/clipping)
        expected_p1 <- p0 * exp(theta)

        expect_equal(p1, expected_p1,
            tolerance = 1e-5,
            label = paste0("Consistency at theta=", round(theta, 2), ", phi=", round(phi, 2))
        )
    }
})

test_that("getProbRR.org boundary conditions (South)", {
    # South Region: phi large negative (e.g., -15), logrr negative but > phi
    phi <- -15
    theta <- -5 # > phi, < 0

    probs <- rbrm::getProbRR.org(theta, phi)
    p0 <- probs$p0
    p1 <- probs$p1

    # Previous bug clamped p1 to 0 here. Now it should be consistent.
    expect_gt(p1, 0)
    expect_equal(p1, p0 * exp(theta), tolerance = 1e-5)
})


test_that("getProbRR.alt maintains p1 = p0 * exp(logrr) relationship", {
    logrr_vals <- seq(-5, 5, length.out = 10)
    logop_vals <- seq(-5, 5, length.out = 10)
    grid <- expand.grid(logrr = logrr_vals, logop = logop_vals)

    for (i in 1:nrow(grid)) {
        theta <- grid$logrr[i]
        phi <- grid$logop[i]

        probs <- rbrm::getProbRR.alt(theta, phi)
        p0 <- probs$p0
        p1 <- probs$p1

        expect_equal(p1, p0 * exp(theta), tolerance = 1e-4) # Slightly looser tolerance for alt due to different calculation path potentially
    }
})
