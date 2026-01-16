#' Predict method for rbrm objects
#'
#' @param object A fitted `rbrm` object.
#' @param new_va Optional matrix of new predictors for alpha.
#' @param new_vb Optional matrix of new predictors for beta.
#' @param type Type of prediction: "response" (probabilities), "link" (linear predictors), "coefficients".
#' @param ... Additional arguments (ignored).
#' @return Matrix of predictions or list of coefficients.
#' @export
predict.rbrm <- function(object, new_va = NULL, new_vb = NULL, type = c("response", "link", "coefficients"), ...) {
    type <- match.arg(type)

    if (type == "coefficients") {
        return(list(alpha = object$alpha, beta = object$beta))
    }

    # Extract Coefficients
    alpha <- object$alpha
    beta <- object$beta

    # Defaults if not provided
    if (is.null(new_va)) {
        cli::cli_abort("Argument `new_va` is missing. Please provide `new_va` (and `new_vb` if applicable) for prediction.")
    }

    if (is.null(new_vb)) new_vb <- new_va

    new_va <- as.matrix(new_va)
    new_vb <- as.matrix(new_vb)

    # Intercept Handling
    # If model has intercept, but input matrix doesn't match coefficient length, add intercept.
    if (object$intercept) {
        # Check alpha
        if (ncol(new_va) == length(alpha) - 1) {
            new_va <- cbind(Intercept = 1, new_va)
        } else if (ncol(new_va) != length(alpha)) {
            cli::cli_abort("Dimension mismatch for `new_va`: Model has {length(alpha)} alpha coefficients (including intercept), but input has {ncol(new_va)} columns.")
        }

        # Check beta
        if (ncol(new_vb) == length(beta) - 1) {
            new_vb <- cbind(Intercept = 1, new_vb)
        } else if (ncol(new_vb) != length(beta)) {
            cli::cli_abort("Dimension mismatch for `new_vb`: Model has {length(beta)} beta coefficients (including intercept), but input has {ncol(new_vb)} columns.")
        }
    } else {
        # No intercept
        if (ncol(new_va) != length(alpha)) {
            cli::cli_abort("Dimension mismatch for `new_va`: Model has {length(alpha)} alpha coefficients, but input has {ncol(new_va)} columns.")
        }
        if (ncol(new_vb) != length(beta)) {
            cli::cli_abort("Dimension mismatch for `new_vb`: Model has {length(beta)} beta coefficients, but input has {ncol(new_vb)} columns.")
        }
    }

    # Linear Predictors
    theta <- new_va %*% alpha
    phi <- new_vb %*% beta

    if (type == "link") {
        return(list(logrr = theta, logop = phi))
    }

    # Response (Probabilities)
    # Defaulting to getProbRR.org as it's the standard for this package
    probs <- getProbRR.org(theta, phi)

    return(probs) # Returns list(p0, p1)
}
