test_that("cv_rbrm runs and returns correct class/structure", {
    set.seed(123)
    n <- 50
    p <- 3
    data <- rbrm::generate_data(pa = p, pb = p, n = n, n_test = 10, alpha = rep(0, p), beta = rep(0, p), gamma = rep(0, p))

    # Run CV with a small grid
    lambda_grid <- c(1, 0.5, 0.1)

    cv_fit <- rbrm::cv_rbrm(
        va = data$v.train,
        vb = data$v.train,
        x = data$x.train,
        y = data$y.train,
        lambda = lambda_grid,
        nfolds = 3,
        prob_fun = rbrm::getProbRR.org,
        n_lambdas = 3
    )

    expect_s3_class(cv_fit, "cv_rbrm")
    expect_named(cv_fit, c(
        "call", "lambda_grid", "cv_mean", "cv_se", "type.measure",
        "lambda.min", "lambda.1se", "lambda.selected", "index",
        "alpha", "beta", "convergence", "step", "relax_lsso_info",
        "relaxation_performed", "relax_factor_effective",
        "cv_metrics_all_folds", "cv_results_matrix", "fold_ids",
        "final_fit_object", "time"
    ))

    expect_true(cv_fit$lambda.min %in% lambda_grid)
    expect_true(cv_fit$lambda.selected %in% lambda_grid)

    # Check coefficients length
    expect_length(cv_fit$alpha, p)
    expect_length(cv_fit$beta, p)
})

test_that("cv_rbrm with relaxed lasso works", {
    set.seed(123)
    n <- 50
    p <- 3
    data <- rbrm::generate_data(pa = p, pb = p, n = n, n_test = 10, alpha = rep(0, p), beta = rep(0, p), gamma = rep(0, p))

    cv_fit_relaxed <- rbrm::cv_rbrm(
        va = data$v.train,
        vb = data$v.train,
        x = data$x.train,
        y = data$y.train,
        lambda = c(0.5, 0.1),
        nfolds = 3,
        relax_lsso = TRUE,
        relax_factors_grid = c(0.5, 1), # Small grid for speed
        prob_fun = rbrm::getProbRR.org
    )

    expect_true(cv_fit_relaxed$relax_lsso_info$performed)
    expect_s3_class(cv_fit_relaxed, "cv_rbrm")
})

test_that("cv_rbrm selection criteria (min vs 1se)", {
    set.seed(123)
    n <- 50
    p <- 3
    data <- rbrm::generate_data(pa = p, pb = p, n = n, n_test = 10, alpha = rep(0, p), beta = rep(0, p), gamma = rep(0, p))

    # Ensure we have enough lambdas to potentially differ
    lambda_grid <- seq(1, 0.1, length.out = 5)

    cv_fit_1se <- rbrm::cv_rbrm(
        va = data$v.train,
        vb = data$v.train,
        x = data$x.train,
        y = data$y.train,
        lambda = lambda_grid,
        nfolds = 3,
        index = "1se",
        prob_fun = rbrm::getProbRR.org
    )

    expect_equal(cv_fit_1se$index, "1se")
    expect_equal(cv_fit_1se$lambda.selected, cv_fit_1se$lambda.1se)

    cv_fit_min <- rbrm::cv_rbrm(
        va = data$v.train,
        vb = data$v.train,
        x = data$x.train,
        y = data$y.train,
        lambda = lambda_grid,
        nfolds = 3,
        index = "min",
        prob_fun = rbrm::getProbRR.org
    )
    expect_equal(cv_fit_min$index, "min")
    expect_equal(cv_fit_min$lambda.selected, cv_fit_min$lambda.min)
})
