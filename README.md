
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rbrm: Regularized Binary Regression Models for Relative Risk

<!-- badges: start -->

<!-- badges: end -->

## Overview

`rbrm` provides a comprehensive framework for estimating **relative
risks** (RR) using penalized binary regression models. Unlike
traditional logistic regression (which estimates odds ratios), this
package directly models relative risk through a novel log odds-product
parameterization, avoiding the limitations of Poisson regression that
can produce invalid probabilities.

## Installation

You can install the development version from
[GitHub](https://github.com/JavierMtzRdz/rbrm):

``` r
# Install remotes if needed
# install.packages("remotes")
remotes::install_github("JavierMtzRdz/rbrm")
```

## Quick Start

``` r
library(rbrm)

# Simulate data
set.seed(123)
n <- 200
p <- 5

data <- generate_data(
  n = n,
  p_a = p,
  p_b = p,
  alpha_true = c(0.5, -0.3, 0.2, 0, 0),
  beta_true = c(-0.2, 0.4, 0, 0, 0)
)

# Fit regularized model
fit <- fit.rbrm(
  va = data$va,
  vb = data$vb,
  x = data$x,
  y = data$y,
  lambda = 0.01,
  intercept = TRUE
)

print(fit)

# Cross-validation for lambda selection
cv_fit <- cv_rbrm(
  va = data$va,
  vb = data$vb,
  x = data$x,
  y = data$y,
  nfold = 5,
  nlambda = 50
)

plot(cv_fit)

# Fit regularization path
path <- rbrm_path(
  va = data$va,
  vb = data$vb,
  x = data$x,
  y = data$y,
  nlambda = 50
)

plot(path)
```

## The Model

This package implements the Richardson-Robins-Wang (2017) binary
regression model for relative risk:

- **Outcome Model**: $P(Y=1 \mid X, W, Z) = p_1$ if $X=1$, $p_0$ if
  $X=0$
- **RR Parameterization**: $\log(\text{RR}) = \theta = W^\top \alpha$
  (treatment effect)
- **Baseline Risk**: $p_0$ via log odds-product $\phi = Z^\top \beta$

The model directly estimates relative risk while ensuring valid
probabilities $(0 \leq p_0, p_1 \leq 1)$.

### Regularization

Supports L1 (lasso) and elastic net penalties for high-dimensional data:

- Variable selection when $p > n$
- Separate penalties for $\alpha$ and $\beta$ parameters
- Unpenalized intercept option

## Available Optimizers

1.  **FISTA** (default): Fast proximal gradient with momentum and
    adaptive restart
2.  **Newton-CD**: Coordinate descent with Newton steps and diagonal
    Hessian
3.  **Active Set Newton-CD**: Efficient for sparse solutions via KKT
    conditions
4.  **L-BFGS**: Quasi-Newton method (best for low-dimensional,
    unpenalized problems)

## References

Richardson, T. S., Robins, J. M., & Wang, L. (2017). On modeling and
estimation for the relative risk and risk difference. *Journal of the
American Statistical Association*, 112(519), 1121-1130.
[arXiv:1510.02430](https://arxiv.org/abs/1510.02430)

## Acknowledgments

This package extends the original
[regularized-RR-regression](https://github.com/ChloeYou/regularized-RR-regression)
implementation and builds upon the
[`brm`](https://github.com/mclements/brm) package.

## License

GPL (\>= 2)
