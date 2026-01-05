#' Compute Probabilities for Robinson's Relative Risk Model specification
#'
#' This function computes the probabilities \( p_0 \) and \( p_1 \) for a binary outcome
#' under an alternative specification of the relative risk regression model.
#'
#' @param logrr A numeric vector or matrix representing the log-relative risk. If `logrr` is
#'   a vector of length 2, the first element is interpreted as the log-relative risk
#'   and the second as the log-odds product. If `logrr` is a matrix with two columns,
#'   the first column is treated as the log-relative risk and the second as the log-odds product.
#' @param logop A numeric value or vector representing the log-odds product. If `logop`
#'   is `NA` and `logrr` is a vector of length 2, `logop` is extracted from the second element of `logrr`.
#'
#' @return A matrix with two columns:
#' \describe{
#'   \item{\code{p0}}{The probability of the binary outcome when the treatment \( X = 0 \).}
#'   \item{\code{p1}}{The probability of the binary outcome when the treatment \( X = 1 \).}
#' }
#'
#' @details The probabilities are computed using the following formulas:
#' \deqn{p_0 = \frac{-\left(1 + \exp(\text{logop}) \cdot (1 + \exp(\text{logrr}))\right)
#'   + \sqrt{\left(1 + \exp(\text{logop}) \cdot (1 + \exp(\text{logrr}))\right)^2
#'   - 4 \cdot \exp(\text{logrr} + 2 \cdot \text{logop})}}{2 \cdot \exp(\text{logrr} + \text{logop})}}
#'
#' \deqn{p_1 = \exp(\text{logrr}) \cdot p_0}
#'
#' The function supports flexible input formats and adjusts automatically when the input is
#' a matrix or vector of appropriate dimensions.
#'
#' @examples
#' # Example with log-relative risk and log-odds product as separate inputs
#' getProbRR.alt(logrr = 0.5, logop = 0.3)
#'
#' # Example with log-relative risk and log-odds product as a vector
#' getProbRR.alt(logrr = c(0.5, 0.3), logop = NA)
#'
#' # Example with log-relative risk and log-odds product as a matrix
#' logrr_matrix <- matrix(c(0.5, 0.3, 0.7, 0.4), ncol = 2)
#' getProbRR.alt(logrr = logrr_matrix)
#'
#' @export
getProbRR.org <- function(logrr, logop = NA,
                          clipping = 1e-15) {
  if (is.matrix(logrr) && ncol(logrr) == 2) {
    logop <- logrr[, 2]
    logrr <- logrr[, 1]
  } else if (length(logop) == 1 && is.na(logop) && length(logrr) == 2) {
    logop <- logrr[2]
    logrr <- logrr[1]
  }

  p0 <- ifelse((logop < (-12)) | (logop > 12) | (logrr < (-12)) | (logrr > 12),
    ## on the boundary South edge: large -ve logrr or (large -ve logop and -ve
    ## logrr)
    ifelse((logrr < (-12)) | ((logop < (-12)) & (logrr < 0)),
      brm:::getPrbAux(logop - logrr),
      ifelse((logrr > 12) | ((logop < (-12)) & (logrr > 0)),
        ## West edge: large +ve logrr or (large -ve logop and +ve logrr)
        0,
        pmin(exp(-logrr), 1)
      )
    ),
    ## not on the boundary logop = 0; solving linear equations logop not 0;
    ## solving a quadratic equation
    ifelse(same(logop, 0),
      1 / (1 + exp(logrr)),
      (-(exp(logrr) + 1) * exp(logop) + sqrt(exp(2 * logop) * (exp(logrr) + 1)^2 + 4 * exp(logrr + logop) * (1 - exp(logop)))) / (2 * exp(logrr) * (1 - exp(logop)))
    )
  )

  p1 <- ifelse((logop < (-12)) | (logop > 12) |
    (logrr < (-12)) | (logrr > 12),
  ## on the boundary South edge: large -ve logrr or (large -ve logop and -ve
  ## logrr
  ifelse((logrr < (-12)) | ((logop < (-12)) &
    (logrr < 0)),
  p0 * exp(logrr),
  ifelse((logrr > 12) | ((logop < (-12)) & (logrr > 0)),
    ## West edge: large +ve logrr or (large -ve logop and +ve logrr)
    brm:::getPrbAux(logop + logrr),
    pmin(exp(logrr), 1)
  )
  ),
  ## not on the boundary logop = 0
  exp(logrr) * p0
  )


  clipping <- ifelse(is.logical(clipping) &&
    isTRUE(clipping), 1e-15, clipping)
  p0 <- pmin(pmax(p0, clipping), 1 - clipping)
  p1 <- pmin(pmax(p1, clipping), 1 - clipping)

  return(structure(
    list(p0 = p0, p1 = p1),
    class = "Richardson"
  ))
}


#' Compute Probabilities for Alternative Relative Risk Model
#'
#' This function computes the probabilities \( p_0 \) and \( p_1 \) for a binary outcome
#' under an alternative specification of the relative risk regression model.
#'
#' @param logrr A numeric vector or matrix representing the log-relative risk. If `logrr` is
#'   a vector of length 2, the first element is interpreted as the log-relative risk
#'   and the second as the log-odds product. If `logrr` is a matrix with two columns,
#'   the first column is treated as the log-relative risk and the second as the log-odds product.
#' @param logop A numeric value or vector representing the log-odds product. If `logop`
#'   is `NA` and `logrr` is a vector of length 2, `logop` is extracted from the second element of `logrr`.
#'
#' @return A matrix with two columns:
#' \describe{
#'   \item{\code{p0}}{The probability of the binary outcome when the treatment \( X = 0 \).}
#'   \item{\code{p1}}{The probability of the binary outcome when the treatment \( X = 1 \).}
#' }
#'
#' @details The probabilities are computed using the following formulas:
#' \deqn{p_0 = \frac{-\left(1 + \exp(\text{logop}) \cdot (1 + \exp(\text{logrr}))\right)
#'   + \sqrt{\left(1 + \exp(\text{logop}) \cdot (1 + \exp(\text{logrr}))\right)^2
#'   - 4 \cdot \exp(\text{logrr} + 2 \cdot \text{logop})}}{2 \cdot \exp(\text{logrr} + \text{logop})}}
#'
#' \deqn{p_1 = \exp(\text{logrr}) \cdot p_0}
#'
#' The function supports flexible input formats and adjusts automatically when the input is
#' a matrix or vector of appropriate dimensions.
#'
#' @examples
#' # Example with log-relative risk and log-odds product as separate inputs
#' getProbRR.alt(logrr = 0.5, logop = 0.3)
#'
#' # Example with log-relative risk and log-odds product as a vector
#' getProbRR.alt(logrr = c(0.5, 0.3), logop = NA)
#'
#' # Example with log-relative risk and log-odds product as a matrix
#' logrr_matrix <- matrix(c(0.5, 0.3, 0.7, 0.4), ncol = 2)
#' getProbRR.alt(logrr = logrr_matrix)
#'
#' @export
getProbRR.alt <- function(logrr, logop,
                          clipping = 1e-15) {
  if (is.matrix(logrr) && ncol(logrr) == 2) {
    logop <- logrr[, 2]
    logrr <- logrr[, 1]
  } else if (length(logop) == 1 && is.na(logop) && length(logrr) == 2) {
    logop <- logrr[2]
    logrr <- logrr[1]
  }

  p0 <- (1 + exp(logop) * (1 + exp(logrr)) -
    sqrt((1 + exp(logop) * (1 + exp(logrr)))^2 - 4 * exp(logrr + 2 * logop))) / (2 * exp(logrr + logop))

  p1 <- exp(logrr) * p0

  clipping <- ifelse(is.logical(clipping) &&
    isTRUE(clipping), 1e-15, clipping)
  p0 <- pmin(pmax(p0, clipping), 1 - clipping)
  p1 <- pmin(pmax(p1, clipping), 1 - clipping)

  attr(c(1, 2), "dim")

  return(structure(
    list(p0 = p0, p1 = p1),
    class = "Pozza"
  ))
}
