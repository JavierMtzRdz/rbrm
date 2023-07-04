file_path <- here::here("code", "scripts")
source(file.path(file_path, "rbrm.r"))
cv_rbrm <- function(va, vb, x, y, lambda = NULL, n_lambdas = 50, nfolds = 3, ...) {
    # nfolds need to be at least 2
    if (nfolds < 2) {
        stop("nfolds must be at least 2")
    }
    n <- length(y)
    fold_ids <- sample(rep_len(1:nfolds, n))

    if (is.null(lambda)) {
        epsilon <- 0.001
        X <- cbind(va) # should we include x?
        # y_bin <- factor(ifelse(y == 0, -1, 1))
        # prp <- prop.table(table(y_bin))
        # max_lambda <- (max(abs(colSums(X * ifelse(y == 0, prp[[1]], prp[[2]])))) / 100)
        max_lambda <- max(abs(cor(X, y)))
        lambda_grid <- rev(exp(seq(log(epsilon * max_lambda),
            log(max_lambda),
            length.out = n_lambdas
        ))) # descending (large to small)
    } else {
        lambda_grid <- lambda
    }
    # print(lambda_grid)

    # array to store the model coefficients for each lambda for each fold
    model_history <- array(0, dim = c(ncol(va) * 2, n_lambdas, nfolds))

    cv_results <- lapply(1:nfolds, function(fold) {
        print(paste("Processing fold", fold, "of", nfolds, sep = " "))
        test_idx <- which(fold_ids == fold)
        train_idx <- which(fold_ids != fold)

        va_train <- va[train_idx, ]
        vb_train <- vb[train_idx, ]
        x_train <- x[train_idx]
        y_train <- y[train_idx]

        va_test <- va[test_idx, ]
        vb_test <- vb[test_idx, ]
        x_test <- x[test_idx]
        y_test <- y[test_idx]

        fold_models <- lapply(c(1:length(lambda_grid)), function(x) {
            current_lambda <- lambda_grid[x]
            cat("Processing lambda:", current_lambda, "\n")
            fit <- rbrm(va_train, vb_train, x_train, y_train, lambda = current_lambda, ...)
            # print(fit$value)
            model_history[, x, fold] <- fit$point.est
            return(fit)
        })

        test_deviance <- sapply(fold_models, function(fit) {
            alpha <- fit$point.est[1:(length(fit$point.est) / 2)]
            beta <- fit$point.est[((length(fit$point.est) / 2) + 1):length(fit$point.est)]
            # print("alpha")
            # print(alpha)
            # print("beta")
            # print(beta)

            logrr <- va_test %*% alpha
            logop <- vb_test %*% beta

            p0 <- brm::getProbRR(logrr, logop)[, 1]
            p1 <- brm::getProbRR(logrr, logop)[, 2]

            # p0[p0 == 1] <- 1 - 1e-15
            # p0[p0 <= 0] <- 1e-15

            # p1[p1 == 1] <- 1 - 1e-15
            # p1[p1 <= 0] <- 1e-15

            fitted.prob <- c(p0[x_test == 0], p1[x_test == 1])
            # print(paste0("fitted.prob: ", fitted.prob))
            true.y <- c(y_test[x_test == 0], y_test[x_test == 1])

            return((-2 * sum(log(fitted.prob[true.y == 1])) - 2 * sum(log1p(-fitted.prob[true.y == 0]))) / length(true.y))
        })

        return(test_deviance)
    }) # Closing bracket for lapply

    cv_mean_deviance <- rowMeans(do.call(cbind, cv_results))
    # print(cv_mean_deviance)

    # also make sure its not choosing the null model
    best_lambda_idx <- which.min(cv_mean_deviance)
    # print(best_lambda_idx)
    # while (all(model_history[[best_lambda_idx]]$point.est[1:ncol(va)] == 0)) {
    #     cat("Lambda", lambda_grid[best_lambda_idx], "results in a null model. Choosing the next lambda.\n")
    #     best_lambda_idx <- best_lambda_idx + 1
    #     if (best_lambda_idx > length(lambda_grid)) {
    #         stop("All lambda values result in a null model.")
    #     }
    # }
    best_fit <- rbrm(va, vb, x, y, lambda = lambda_grid[best_lambda_idx], ...)


    list(
        lambda = lambda_grid[best_lambda_idx],
        alpha = best_fit$point.est[1:ncol(va)],
        beta = best_fit$point.est[(ncol(va) + 1):(ncol(va) + ncol(vb))],
        convergence = best_fit$convergence,
        step = best_fit$step,
        cv_mean_deviance = cv_mean_deviance,
        fold_deviances = cv_results,
        best_lambda_idx = best_lambda_idx,
        model_history = model_history,
        lambda_grid = lambda_grid
    )
}
