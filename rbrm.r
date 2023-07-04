# va and vb should not include the 1 columns as intercept
# the default currently assumes there are no intercepts given in va and vb
rbrm <- function(va, vb, y, x, alpha.start = NULL, beta.start = NULL,
                 max.step = 3000, thres = 1e-04, lambda = 0,
                 lr.alpha = 0.06, lr.beta = 0.02,
                 intercept = TRUE, early_stopping_rounds = 10) {
    soft_thres <- function(x, lambda) {
        sx <- abs(x) - lambda
        sx[sx < 0] <- 0
        return(sx * sign(x))
    }

    if (is.null(vb)) {
        vb <- va
    }
    pa <- dim(va)[2]
    pb <- dim(vb)[2]

    # sanity check for the intercept term in va, vb
    if (all(va[, 1] == 1) & all(vb[, 1] == 1)) {
        intercept <- TRUE
    }

    ## starting values for parameter optimization
    if (is.null(alpha.start)) {
        alpha.start <- c(rep(0, pa))
    }
    if (is.null(beta.start)) {
        beta.start <- c(rep(0.01, pb))
    }

    nllh.alpha <- function(alpha) {
        logrr <- (va %*% alpha)
        logop <- (vb %*% beta)

        p0 <- brm::getProbRR(logrr, logop)[, 1]
        p1 <- brm::getProbRR(logrr, logop)[, 2]

        p0[p0 == 1] <- 1 - 1e-15
        p0[p0 <= 0] <- 1e-15

        p1[p1 == 1] <- 1 - 1e-15
        p1[p1 <= 0] <- 1e-15


        if (any(is.nan(log1p(-p1[x == 1])))) {
            p1[which(is.nan(log1p(-p1[x == 1])))] <- 1 - 1e-15
        }
        if (any(is.nan(log1p(-p0[x == 0])))) {
            p1[which(is.nan(log1p(-p0[x == 0])))] <- 1 - 1e-15
        }

        return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) +
            (y[x == 0]) * log(p0[x == 0])) -
            sum((1 - y[x == 1]) * log(1 - p1[x == 1]) +
                (y[x == 1]) * log(p1[x == 1])))
    }

    nllh.beta <- function(beta) {
        logrr <- (va %*% alpha)
        logop <- (vb %*% beta)

        p0 <- brm::getProbRR(logrr, logop)[, 1]
        p1 <- brm::getProbRR(logrr, logop)[, 2]

        p0[p0 == 1] <- 1 - 1e-15
        p0[p0 <= 0] <- 1e-15

        p1[p1 == 1] <- 1 - 1e-15
        p1[p1 <= 0] <- 1e-15

        if (any(is.nan(log1p(-p1[x == 1])))) {
            p1[which(is.nan(log1p(-p1[x == 1])))] <- 1 - 1e-15
        }
        if (any(is.nan(log1p(-p0[x == 0])))) {
            p1[which(is.nan(log1p(-p0[x == 0])))] <- 1 - 1e-15
        }

        return(-sum((1 - y[x == 0]) * log1p(-p0[x == 0]) +
            (y[x == 0]) * log(p0[x == 0])) -
            sum((1 - y[x == 1]) * log1p(-p1[x == 1]) +
                (y[x == 1]) * log(p1[x == 1])))
    }

    ## FISTA proximal gradient descent functions
    proximal.gd.alpha.fista <- function(alpha, step_size, lambda, t_old, last_alpha) {
        gradient <- numDeriv::grad(nllh.alpha, alpha, method = "simple")
        gradient[is.na(gradient)] <- 0
        input <- alpha - step_size * gradient
        alpha_new <- soft_thres(input, lambda * step_size)
        if (intercept == TRUE) {
            alpha_new[1] <- input[1]
        }
        t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
        y_alpha_new <- alpha_new + (t_old - 1) / t_new * (alpha_new - last_alpha)
        return(list(alpha_new = alpha_new, t_alpha = t_new, y_alpha = y_alpha_new))
    }

    proximal.gd.beta.fista <- function(beta, step_size, lambda, t_old, last_beta) {
        gradient <- numDeriv::grad(nllh.beta, beta, method = "simple")
        gradient[is.na(gradient)] <- 0
        input <- beta - step_size * gradient
        beta_new <- soft_thres(input, lambda * step_size)
        if (intercept == TRUE) {
            beta_new[1] <- input[1]
        }
        t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
        y_beta_new <- beta_new + (t_old - 1) / t_new * (beta_new - last_beta)
        return(list(beta_new = beta_new, t_beta = t_new, y_beta = y_beta_new))
    }

    penalized.neg.log.likelihood <- function(pars) {
        alpha <- pars[1:pa]
        beta <- pars[(pa + 1):(pa + pb)]

        logrr <- (va %*% alpha)
        logop <- (vb %*% beta)

        p0 <- brm::getProbRR(logrr, logop)[, 1]
        p1 <- brm::getProbRR(logrr, logop)[, 2]

        p0[p0 == 1] <- 1 - 1e-15
        p0[p0 <= 0] <- 1e-15

        p1[p1 == 1] <- 1 - 1e-15
        p1[p1 <= 0] <- 1e-15

        if (any(is.nan(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) +
            (y[x == 0]) * log(p0[x == 0])) -
            sum((1 - y[x == 1]) * log(1 - p1[x == 1]) +
                (y[x == 1]) * log(p1[x == 1]))))) {
            stop("NaN values encountered in unpenalized.nllh. Please check the input matrices.")
        }
        unpenalized.nllh <- (-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) +
            (y[x == 0]) * log(p0[x == 0])) -
            sum((1 - y[x == 1]) * log(1 - p1[x == 1]) +
                (y[x == 1]) * log(p1[x == 1])))
        if (intercept == TRUE) {
            penalty <- lambda * (sum(abs(alpha[-1])) + sum(abs(beta[-1]))) # penalty does not apply to the intercept coefficient
        } else {
            penalty <- lambda * (sum(abs(alpha)) + sum(abs(beta)))
        }
        return(unpenalized.nllh + penalty)
    }

    ## Optimization
    alpha <- alpha.start
    beta <- beta.start
    last_alpha <- last_beta <- alpha
    t_alpha <- t_beta <- 1
    y_alpha <- y_beta <- alpha.start
    step_size_alpha <- lr.alpha
    step_size_beta <- lr.beta
    step <- 0
    for (iter in 1:max.step) {
        # print step size every 100 steps
        # if (step %% 100 == 0) {
        #     print(paste0("step ", step))
        #     # print("alpha step size")
        #     # print(t_alpha)
        #     # print("beta step size")
        #     # print(t_beta)
        # }
        step <- step + 1
        # FISTA update for alpha
        last_alpha <- alpha
        res_alpha <- proximal.gd.alpha.fista(y_alpha, step_size_alpha, lambda, t_alpha, last_alpha)
        alpha_new <- res_alpha$alpha_new
        t_alpha <- res_alpha$t_alpha
        y_alpha <- res_alpha$y_alpha

        # FISTA update for beta
        last_beta <- beta
        res_beta <- proximal.gd.beta.fista(y_beta, step_size_beta, lambda, t_beta, last_beta)
        beta_new <- res_beta$beta_new
        t_beta <- res_beta$t_beta
        y_beta <- res_beta$y_beta

        grad_alpha <- numDeriv::grad(nllh.alpha, y_alpha, method = "simple")
        grad_beta <- numDeriv::grad(nllh.beta, y_beta, method = "simple")
        if ((max(abs(grad_alpha)) < thres) & (max(abs(grad_beta)) < thres)) {
            break
        }

        # Update alpha and beta for the next iteration
        alpha <- alpha_new
        beta <- beta_new
    }
    opt <- list(
        point.est = c(alpha, beta), convergence = (step < max.step),
        value = penalized.neg.log.likelihood(c(alpha, beta)),
        step = step
    )

    return(opt)
}
