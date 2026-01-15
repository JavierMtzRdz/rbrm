// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(fntl)]]
#include <R_ext/Applic.h>  // For lbfgsb
#include <RcppArmadillo.h> // For Armadillo matrix/vector operations
#include <cmath>
#include <fntl.h>
#include <limits>
#include <string>
#include <vector>

// Internal helper function (not exported to R)
// Used by getProbRR_org_cpp
double getPrbAux_cpp(double x) {
  if (R_IsNaN(x)) { // Explicitly check for NaN input
    return NA_REAL;
  }
  if ((x < 17.0) && (x > -500.0)) {
    return 0.5 * std::exp(x) * (-1.0 + std::sqrt(1.0 + 4.0 * std::exp(-x)));
  } else if (x < 0.0) {
    return 0.0;
  } else { // x >= 17.0 or x <= -500.0 (and x is not NaN)
    return 1.0;
  }
}

//' @export
// [[Rcpp::export]]
Rcpp::List getProbRR_org_cpp(const Rcpp::NumericVector &logrr,
                             const Rcpp::NumericVector &logop,
                             double clipping = 1e-15) {
  int n = logrr.size();
  if (logop.size() != n) {
    Rcpp::stop(
        "In prob_rr_org_worker: logrr and logop must have the same length.");
  }

  Rcpp::NumericVector p0_out(n);
  Rcpp::NumericVector p1_out(n);

  const double boundary_limit_low = -12.0;
  const double boundary_limit_high = 12.0;
  const double zero_epsilon = std::sqrt(std::numeric_limits<double>::epsilon());

  for (int i = 0; i < n; ++i) {
    double lr = logrr[i];
    double lo = logop[i];

    // Explicit NA check for inputs for this iteration
    if (R_IsNaN(lr) || R_IsNaN(lo)) {
      p0_out[i] = NA_REAL;
      p1_out[i] = NA_REAL;
      continue; // Skip to the next iteration
    }

    double current_p0 = NA_REAL;
    double current_p1 = NA_REAL;

    bool on_boundary = (lo < boundary_limit_low) ||
                       (lo > boundary_limit_high) ||
                       (lr < boundary_limit_low) || (lr > boundary_limit_high);

    if (on_boundary) {
      if ((lr < boundary_limit_low) ||
          ((lo < boundary_limit_low) && (lr < 0.0))) {
        current_p0 = getPrbAux_cpp(lo - lr);
        current_p1 = 0.0;
      } else if ((lr > boundary_limit_high) ||
                 ((lo < boundary_limit_low) && (lr > 0.0))) {
        current_p0 = 0.0;
        current_p1 = getPrbAux_cpp(lo + lr);
      } else {
        current_p0 = std::fmin(std::exp(-lr), 1.0);
        current_p1 = std::fmin(std::exp(lr), 1.0);
      }
    } else {
      if (std::fabs(lo) < zero_epsilon) {
        current_p0 = 1.0 / (1.0 + std::exp(lr));
      } else {
        double exp_lr = std::exp(lr);
        double exp_lo = std::exp(lo);
        double term1_coeff = exp_lr + 1.0;
        double expm1_lo = std::expm1(lo); // Use expm1 for precision
        double discriminant = std::exp(2.0 * lo) * term1_coeff * term1_coeff +
                              4.0 * std::exp(lr + lo) * (-expm1_lo);
        double sqrt_discriminant = std::sqrt(discriminant);
        double numerator_val = -term1_coeff * exp_lo + sqrt_discriminant;
        double denominator_p0 = 2.0 * exp_lr * (-expm1_lo);

        if (std::fabs(denominator_p0) < (zero_epsilon * zero_epsilon)) {
          if (R_IsNaN(numerator_val)) {
            current_p0 = NA_REAL;
          } else if (std::fabs(numerator_val) < (zero_epsilon * zero_epsilon)) {
            current_p0 = NA_REAL;
          } else if (numerator_val > 0) {
            current_p0 = R_PosInf;
          } else {
            current_p0 = R_NegInf;
          }
        } else {
          current_p0 = numerator_val / denominator_p0;
        }
      }
      current_p1 = std::exp(lr) * current_p0;
    }

    if (R_IsNaN(current_p0)) {
      // p0 remains NA_REAL
    } else if (!R_finite(current_p0)) {
      current_p0 = (current_p0 > 0) ? (1.0 - clipping) : clipping;
    } else {
      current_p0 = std::fmin(std::fmax(current_p0, clipping), 1.0 - clipping);
    }

    if (R_IsNaN(current_p1)) {
      // p1 remains NA_REAL
    } else if (!R_finite(current_p1)) {
      current_p1 = (current_p1 > 0) ? (1.0 - clipping) : clipping;
    } else {
      current_p1 = std::fmin(std::fmax(current_p1, clipping), 1.0 - clipping);
    }

    p0_out[i] = current_p0;
    p1_out[i] = current_p1;
  }

  return Rcpp::List::create(Rcpp::Named("p0") = p0_out,
                            Rcpp::Named("p1") = p1_out);
}

//' @export
// [[Rcpp::export]]
Rcpp::List getProbRR_alt_cpp(const Rcpp::NumericVector &logrr,
                             const Rcpp::NumericVector &logop,
                             double clipping = 1e-15) {
  int n = logrr.size();
  if (logop.size() != n) {
    Rcpp::stop(
        "In prob_rr_alt_worker: logrr and logop must have the same length.");
  }

  Rcpp::NumericVector p0_out(n);
  Rcpp::NumericVector p1_out(n);

  for (int i = 0; i < n; ++i) {
    double lr = logrr[i];
    double lo = logop[i];

    // Explicit NA check for inputs for this iteration
    if (R_IsNaN(lr) || R_IsNaN(lo)) {
      p0_out[i] = NA_REAL;
      p1_out[i] = NA_REAL;
      continue; // Skip to the next iteration
    }

    double current_p0 =
        NA_REAL; // Initialize just in case, though it should be set.

    double exp_lr = std::exp(lr);
    double exp_lo = std::exp(lo);
    double exp_lr_plus_2lo = std::exp(lr + 2.0 * lo);

    double term_A = 1.0 + exp_lo * (1.0 + exp_lr);
    double term_B_sqrt_arg = term_A * term_A - 4.0 * exp_lr_plus_2lo;

    double term_B = std::sqrt(term_B_sqrt_arg);
    double denominator = 2.0 * std::exp(lr + lo);

    current_p0 = (term_A - term_B) / denominator;

    double current_p1 = exp_lr * current_p0;

    if (R_IsNaN(current_p0)) {
      // p0 remains NA_REAL
    } else if (!R_finite(current_p0)) {
      current_p0 = (current_p0 > 0) ? (1.0 - clipping) : clipping;
    } else {
      current_p0 = std::fmin(std::fmax(current_p0, clipping), 1.0 - clipping);
    }

    if (R_IsNaN(current_p1)) {
      // p1 remains NA_REAL
    } else if (!R_finite(current_p1)) {
      current_p1 = (current_p1 > 0) ? (1.0 - clipping) : clipping;
    } else {
      current_p1 = std::fmin(std::fmax(current_p1, clipping), 1.0 - clipping);
    }

    p0_out[i] = current_p0;
    p1_out[i] = current_p1;
  }

  return Rcpp::List::create(Rcpp::Named("p0") = p0_out,
                            Rcpp::Named("p1") = p1_out);
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector soft_thres_cpp(Rcpp::NumericVector x, double lambda) {
  if (lambda < 0) {
    // In C++ FISTA, you might handle this error differently or ensure lambda is
    // valid.
    Rcpp::warning(
        "lambda in soft_thres_cpp should be non-negative. Using abs(lambda).");
    lambda = std::fabs(lambda);
  }
  int n = x.size();
  Rcpp::NumericVector result(n);

  for (int i = 0; i < n; ++i) {
    double val_x = x[i];
    if (R_IsNaN(val_x)) {
      result[i] = NA_REAL;
    } else if (val_x ==
               0.0) { // Explicitly handle sign(0)*fmax(0-lambda,0) -> 0
      result[i] = 0.0;
    } else {
      double abs_val_x = std::fabs(val_x);
      double sign_val_x =
          (val_x > 0) - (val_x < 0); // Efficient sign: 1, -1, or 0
      result[i] = sign_val_x * std::fmax(0.0, abs_val_x - lambda);
    }
  }
  return result;
}

//' @export
// [[Rcpp::export]]
double nllh_cpp(const arma::vec &alpha, const arma::vec &beta,
                const arma::mat &va, const arma::mat &vb,
                const Rcpp::NumericVector &x_indicator, // Assumed to be 0/1
                const Rcpp::NumericVector &y_outcome,   // Assumed to be 0/1
                int prob_fun, double clipping = 1e-10) {

  int n = y_outcome.size();
  if (n == 0) {
    return 0.0; // Or R_PosInf or an error, consistent with R version if n can
                // be 0
  }
  if (x_indicator.size() != n) {
    Rcpp::stop("In nllh_cpp: x_indicator length must match y_outcome length.");
  }

  arma::vec logrr_arma(n, arma::fill::zeros);
  arma::vec logop_arma(n, arma::fill::zeros);

  if (alpha.n_elem > 0) {
    if (va.n_rows == n && va.n_cols == alpha.n_elem) {
      logrr_arma = va * alpha;
    } else {
      // Handle mismatched dimensions robustly if this function could be called
      // with them
      Rcpp::Rcout << "Warning: Dimension mismatch or empty va for alpha in "
                     "nllh_cpp. va rows: "
                  << va.n_rows << ", va cols: " << va.n_cols
                  << ", alpha elems: " << alpha.n_elem << ", n: " << n
                  << std::endl;
      // Depending on desired behavior, either stop or ensure logrr_arma remains
      // zeros. For now, it remains zeros if condition not met.
      if (va.n_rows != n && alpha.n_elem > 0)
        Rcpp::stop("va rows incorrect for n");
      if (va.n_cols != alpha.n_elem && alpha.n_elem > 0)
        Rcpp::stop("va cols incorrect for alpha");
    }
  }

  if (beta.n_elem > 0) {
    if (vb.n_rows == n && vb.n_cols == beta.n_elem) {
      logop_arma = vb * beta;
    } else {
      Rcpp::Rcout << "Warning: Dimension mismatch or empty vb for beta in "
                     "nllh_cpp. vb rows: "
                  << vb.n_rows << ", vb cols: " << vb.n_cols
                  << ", beta elems: " << beta.n_elem << ", n: " << n
                  << std::endl;
      if (vb.n_rows != n && beta.n_elem > 0)
        Rcpp::stop("vb rows incorrect for n");
      if (vb.n_cols != beta.n_elem && beta.n_elem > 0)
        Rcpp::stop("vb cols incorrect for beta");
    }
  }

  Rcpp::NumericVector logrr_rcpp = Rcpp::wrap(logrr_arma);
  Rcpp::NumericVector logop_rcpp = Rcpp::wrap(logop_arma);

  Rcpp::List ps;
  // Pass 1e-15 to prob fun as base clipping, but nllh clips to 'clipping'
  // (1e-10) Actually, R nllh relies on prob_fun clipping to 'clipping' (1e-10).
  // So we should pass 'clipping' to getProbRR call here.
  if (prob_fun == 0) {
    ps = getProbRR_org_cpp(logrr_rcpp, logop_rcpp, clipping);
  } else if (prob_fun == 1) {
    ps = getProbRR_alt_cpp(logrr_rcpp, logop_rcpp, clipping);
  } else {
    Rcpp::stop("Invalid prob_fun in nllh_cpp");
  }

  Rcpp::NumericVector p0_vec = Rcpp::as<Rcpp::NumericVector>(ps["p0"]);
  Rcpp::NumericVector p1_vec = Rcpp::as<Rcpp::NumericVector>(ps["p1"]);

  double nll_sum = 0.0;

  for (int i = 0; i < n; ++i) {
    double p0 = p0_vec[i];
    double p1 = p1_vec[i];

    // prob_fun already clipped to 'clipping'.
    // R nllh does NOT clip again.
    // So we don't clip again here unless we want safety against NA.

    // Safety against NaN/Inf already handled in getProbRR clipping (if we trust
    // it) But let's be safe against 0 or 1 if clipping didn't catch something.
    // However, if clipping >= 1e-15, it should be fine.

    double yi = y_outcome[i];

    if (x_indicator[i] == 0) {
      if (R_IsNaN(p0) || p0 <= 0 || (1.0 - p0) <= 0) {
        // If p0 is 0 or 1 despite clipping, clamp it.
        // This handles cases where clipping might be 0.0 if user passed 0.
        if (p0 <= 0)
          p0 = 1e-15;
        if (p0 >= 1)
          p0 = 1.0 - 1e-15;
        if (R_IsNaN(p0)) {
          nll_sum = NA_REAL;
          break;
        }
      }
      nll_sum -= (yi * std::log(p0) + (1.0 - yi) * std::log(1.0 - p0));
    } else { // x_indicator[i] == 1
      if (R_IsNaN(p1) || p1 <= 0 || (1.0 - p1) <= 0) {
        if (p1 <= 0)
          p1 = 1e-15;
        if (p1 >= 1)
          p1 = 1.0 - 1e-15;
        if (R_IsNaN(p1)) {
          nll_sum = NA_REAL;
          break;
        }
      }
      nll_sum -= (yi * std::log(p1) + (1.0 - yi) * std::log(1.0 - p1));
    }
    if (R_IsNaN(nll_sum))
      break; // Propagate NaN early
  }

  if (!R_finite(nll_sum)) {
    return R_PosInf;
  }

  // Ensure n is not zero to prevent division by zero if that's possible
  if (n > 0) {
    return nll_sum / static_cast<double>(n);
  } else {
    return 0.0; // Or appropriate value for n=0, R code returns 0/0 = NaN -> Inf
  }
}

//' @export
// [[Rcpp::export]]
double penalized_nllh_cpp(const arma::vec &alpha, const arma::vec &beta,
                          const arma::mat &va, const arma::mat &vb,
                          const Rcpp::NumericVector &x_indicator,
                          const Rcpp::NumericVector &y_outcome, double lambda,
                          bool intercept, int prob_fun,
                          double clipping = 1e-10) {

  double unpenalized_nllh =
      nllh_cpp(alpha, beta, va, vb, x_indicator, y_outcome, prob_fun, clipping);

  // If nllh is Inf (or NaN which becomes Inf), propagate it.
  if (unpenalized_nllh == R_PosInf || !R_finite(unpenalized_nllh)) {
    return R_PosInf;
  }

  double penalty = 0.0;

  if (lambda > 0) { // Only calculate penalty if lambda is non-zero
    double l1_norm_alpha = 0.0;
    if (alpha.n_elem > 0) {
      int start_idx_alpha = (intercept && alpha.n_elem > 0) ? 1 : 0;
      for (arma::uword i = start_idx_alpha; i < alpha.n_elem; ++i) {
        l1_norm_alpha += std::fabs(alpha[i]);
      }
    }

    double l1_norm_beta = 0.0;
    if (beta.n_elem > 0) {
      // Assuming beta might also have an intercept conceptually, or always
      // penalize all of beta R code was sum(abs(beta[if(intercept) -1 else
      // TRUE])) This implies if intercept=TRUE, beta[1] is also an intercept.
      // Adjust if this assumption is wrong. If beta has no intercept,
      // start_idx_beta should always be 0.
      int start_idx_beta = (intercept && beta.n_elem > 0) ? 1 : 0;
      for (arma::uword i = start_idx_beta; i < beta.n_elem; ++i) {
        l1_norm_beta += std::fabs(beta[i]);
      }
    }
    penalty = lambda * (l1_norm_alpha + l1_norm_beta);
  }

  return unpenalized_nllh + penalty;
}

// Helper: check if value is approximately zero
inline bool same_cpp(double x, double y = 0.0, double eps = 1e-8) {
  return std::fabs(x - y) < eps;
}

// Derivative of p0 with respect to theta (Richardson parameterization)
// Port of R function dp0_theta
arma::vec dp0_theta_richardson_cpp(const arma::vec &theta,
                                   const arma::vec &phi) {
  int n = theta.n_elem;
  arma::vec dp0_dtheta(n);

  for (int i = 0; i < n; ++i) {
    double t = theta[i];
    double p = phi[i];
    double x = p - t;
    double expm1_t = std::expm1(t);
    double expm1_p = std::expm1(p);

    // Boundary conditions
    bool is_boundary = (p < -12) || (p > 12) || (t < -12) || (t > 12);
    bool is_south_edge = (t < -12) || ((p < -12) && (t < 0));
    bool is_west_edge = (t > 12) || ((p < -12) && (t > 0));
    bool is_phi_zero = same_cpp(p, 0.0);

    if (is_boundary) {
      if (is_south_edge) {
        if ((x < 17) && (x > -500)) {
          double y_val = 4.0 * std::exp(-x);
          double sqrt_term = std::sqrt(y_val + 1.0);
          double term1 = -y_val / (1.0 + sqrt_term);
          dp0_dtheta[i] =
              (0.5 * term1 * sqrt_term * std::exp(x) + 1.0) / sqrt_term;
        } else {
          dp0_dtheta[i] = 0.0;
        }
      } else if (is_west_edge) {
        dp0_dtheta[i] = 0.0;
      } else {
        double val = -std::exp(-t);
        dp0_dtheta[i] = (val < -1.0) ? 0.0 : val;
      }
    } else {
      // Not on boundary
      if (is_phi_zero) {
        double exp_t = std::exp(t);
        dp0_dtheta[i] = -exp_t / std::pow(1.0 + exp_t, 2);
      } else {
        // Quadratic equation case
        double sqrt_term = std::sqrt(4.0 * std::exp(p + t) +
                                     expm1_t * expm1_t * std::exp(2.0 * p));
        dp0_dtheta[i] = -std::exp(p - t) / (2.0 * expm1_p) +
                        std::exp(p) / (expm1_p * sqrt_term) +
                        (std::exp(-t) * (-expm1_t) * std::exp(2.0 * p)) /
                            (2.0 * expm1_p * sqrt_term);
      }
    }
  }
  return dp0_dtheta;
}

// Derivative of p0 with respect to phi (Richardson parameterization)
// Port of R function dp0_phi
arma::vec dp0_phi_richardson_cpp(const arma::vec &theta, const arma::vec &phi) {
  int n = theta.n_elem;
  arma::vec dp0_dphi(n);

  for (int i = 0; i < n; ++i) {
    double t = theta[i];
    double p = phi[i];
    double x = p - t;
    double expm1_t = std::expm1(t);
    double expm1_p = std::expm1(p);

    // Boundary conditions
    bool is_boundary = (p < -12) || (p > 12) || (t < -12) || (t > 12);
    bool is_south_edge = (t < -12) || ((p < -12) && (t < 0));
    bool is_west_edge = (t > 12) || ((p < -12) && (t > 0));
    (void)is_west_edge; // Silence unused warning
    bool is_phi_zero = same_cpp(p, 0.0);

    if (is_boundary) {
      if (is_south_edge) {
        if ((x < 17) && (x > -500)) {
          double y_val = 4.0 * std::exp(-x);
          double sqrt_term = std::sqrt(y_val + 1.0);
          double term1 = -y_val / (1.0 + sqrt_term);
          dp0_dphi[i] = (0.5 * std::exp(x) * term1 + 1.0) / sqrt_term;
        } else {
          dp0_dphi[i] = 0.0;
        }
      } else {
        dp0_dphi[i] = 0.0;
      }
    } else {
      // Not on boundary
      if (is_phi_zero) {
        double exp_t = std::exp(t);
        dp0_dphi[i] = exp_t / std::pow(exp_t + 1.0, 3);
      } else {
        // Quadratic equation case
        double exp_t = std::exp(t);
        double exp_p = std::exp(p);
        double sqrt_term = std::sqrt(4.0 * std::exp(p + t) +
                                     expm1_t * expm1_t * std::exp(2.0 * p));
        dp0_dphi[i] =
            -((exp_t + 1.0) * exp_p) / (2.0 * exp_t * expm1_p * expm1_p) +
            exp_p / (expm1_p * expm1_p * sqrt_term) +
            (std::exp(-t) * (std::exp(2.0 * t) + 1.0) * std::exp(2.0 * p)) /
                (2.0 * expm1_p * expm1_p * sqrt_term);
      }
    }
  }
  return dp0_dphi;
}

// Derivative of p0 with respect to theta (Pozza/Alt parameterization)
arma::vec dp0_theta_pozza_cpp(const arma::vec &theta, const arma::vec &phi) {
  int n = theta.n_elem;
  arma::vec dp0_dtheta(n);

  for (int i = 0; i < n; ++i) {
    double t = theta[i];
    double p = phi[i];

    // R Formula: sqrt(((exp(theta) + 1) * exp(phi) + 1)^2 - 4 * exp(2 * phi +
    // theta))
    double exp_t = std::exp(t);
    double exp_p = std::exp(p);
    double term_acc = (exp_t + 1.0) * exp_p + 1.0;
    double sqrt_term =
        std::sqrt(term_acc * term_acc - 4.0 * std::exp(2.0 * p + t));

    // Original formula is split:
    // result = term1 - term2

    double part1_num =
        2.0 * term_acc * std::exp(p + t) - 4.0 * std::exp(2.0 * p + t);
    double part1 = part1_num / (2.0 * sqrt_term);

    double num1 = std::exp(p + t) - part1;
    double term1 = (std::exp(-p - t) * num1) / 2.0;

    double term2_inner = -sqrt_term + term_acc;
    double term2 = (std::exp(-p - t) * term2_inner) / 2.0;

    dp0_dtheta[i] = term1 - term2;
  }
  return dp0_dtheta;
}

// Derivative of p0 with respect to phi (Pozza/Alt parameterization)
arma::vec dp0_phi_pozza_cpp(const arma::vec &theta, const arma::vec &phi) {
  int n = theta.n_elem;
  arma::vec dp0_dphi(n);

  for (int i = 0; i < n; ++i) {
    double t = theta[i];
    double p = phi[i];

    double exp_t = std::exp(t);
    double exp_p = std::exp(p);
    double term_acc = (exp_t + 1.0) * exp_p + 1.0;
    double sqrt_term =
        std::sqrt(term_acc * term_acc - 4.0 * std::exp(2.0 * p + t));

    // R Formula:
    double partA =
        2.0 * (exp_t + 1.0) * exp_p * term_acc - 8.0 * std::exp(2.0 * p + t);
    double partB = partA / (2.0 * sqrt_term);
    double numeratorC = (exp_t + 1.0) * exp_p - partB;
    double term1 = (std::exp(-p - t) * numeratorC) / 2.0;

    double term2 = (std::exp(-p - t) * (-sqrt_term + term_acc)) / 2.0;

    dp0_dphi[i] = term1 - term2;
  }
  return dp0_dphi;
}

// Analytical gradient for alpha using correct gradient formula
//' @export
// [[Rcpp::export]]
arma::vec
grad_nll_alpha_analytical_cpp(const arma::vec &alpha, const arma::vec &beta,
                              const arma::mat &va, const arma::mat &vb,
                              const Rcpp::NumericVector &x_indicator,
                              const Rcpp::NumericVector &y_outcome,
                              int prob_fun_selector, double clipping = 1e-10) {

  int n = y_outcome.size();
  int pa = alpha.n_elem;

  if (pa == 0)
    return arma::vec();

  // Compute probabilities
  arma::vec theta = va * alpha;
  arma::vec phi = vb * beta;

  Rcpp::NumericVector theta_rcpp = Rcpp::wrap(theta);
  Rcpp::NumericVector phi_rcpp = Rcpp::wrap(phi);

  Rcpp::List ps;
  if (prob_fun_selector == 0) {
    ps = getProbRR_org_cpp(theta_rcpp, phi_rcpp, clipping);
  } else {
    ps = getProbRR_alt_cpp(theta_rcpp, phi_rcpp, clipping);
  }

  Rcpp::NumericVector p0_vec = Rcpp::as<Rcpp::NumericVector>(ps["p0"]);
  Rcpp::NumericVector p1_vec = Rcpp::as<Rcpp::NumericVector>(ps["p1"]);

  // Compute derivatives of p0 w.r.t. theta
  arma::vec dp0_dtheta;
  if (prob_fun_selector == 0) {
    dp0_dtheta = dp0_theta_richardson_cpp(theta, phi);
  } else {
    dp0_dtheta = dp0_theta_pozza_cpp(theta, phi);
  }

  //  Compute dp1/dtheta = (p0 + dp0/dtheta) * exp(theta)
  arma::vec dp1_dtheta(n);
  for (int i = 0; i < n; ++i) {
    double p0_i = p0_vec[i];
    if (same_cpp(p0_i + dp0_dtheta[i], 0.0)) {
      dp1_dtheta[i] = 0.0;
    } else {
      dp1_dtheta[i] = (p0_i + dp0_dtheta[i]) * std::exp(theta[i]);
    }
  }

  // Compute gradient
  arma::vec inner_alpha(n);
  for (int i = 0; i < n; ++i) {
    double p0_i = std::fmax(clipping, std::fmin(1.0 - clipping, p0_vec[i]));
    double p1_i = std::fmax(clipping, std::fmin(1.0 - clipping, p1_vec[i]));
    double x_i = x_indicator[i];
    double y_i = y_outcome[i];

    double p_i = (x_i == 0) ? p0_i : p1_i;

    // Check if clipped
    bool is_clipped = (p_i <= clipping) || (p_i >= (1.0 - clipping));

    if (is_clipped) {
      inner_alpha[i] = 0.0;
    } else {
      // dllh/dp1 * dp1/dtheta + dllh/dp0 * dp0/dtheta
      double dllh_dp1 = (y_i * x_i) / p1_i - ((1.0 - y_i) * x_i) / (1.0 - p1_i);
      double dllh_dp0 = (y_i * (1.0 - x_i)) / p0_i -
                        ((1.0 - y_i) * (1.0 - x_i)) / (1.0 - p0_i);
      inner_alpha[i] = dllh_dp1 * dp1_dtheta[i] + dllh_dp0 * dp0_dtheta[i];
    }
  }

  arma::vec grad_alpha = -(va.t() * inner_alpha) / static_cast<double>(n);

  return grad_alpha;
}

// Analytical gradient for beta
//' @export
// [[Rcpp::export]]
arma::vec grad_nll_beta_analytical_cpp(const arma::vec &alpha,
                                       const arma::vec &beta,
                                       const arma::mat &va, const arma::mat &vb,
                                       const Rcpp::NumericVector &x_indicator,
                                       const Rcpp::NumericVector &y_outcome,
                                       int prob_fun_selector,
                                       double clipping = 1e-10) {

  int n = y_outcome.size();
  int pb = beta.n_elem;

  if (pb == 0)
    return arma::vec();

  // Compute probabilities
  arma::vec theta = va * alpha;
  arma::vec phi = vb * beta;

  Rcpp::NumericVector theta_rcpp = Rcpp::wrap(theta);
  Rcpp::NumericVector phi_rcpp = Rcpp::wrap(phi);

  Rcpp::List ps;
  if (prob_fun_selector == 0) {
    ps = getProbRR_org_cpp(theta_rcpp, phi_rcpp, clipping);
  } else {
    ps = getProbRR_alt_cpp(theta_rcpp, phi_rcpp, clipping);
  }

  Rcpp::NumericVector p0_vec = Rcpp::as<Rcpp::NumericVector>(ps["p0"]);
  Rcpp::NumericVector p1_vec = Rcpp::as<Rcpp::NumericVector>(ps["p1"]);

  // Compute derivatives of p0 w.r.t. phi
  arma::vec dp0_dphi;
  if (prob_fun_selector == 0) {
    dp0_dphi = dp0_phi_richardson_cpp(theta, phi);
  } else {
    dp0_dphi = dp0_phi_pozza_cpp(theta, phi);
  }

  // Compute dp1/dphi = dp0/dphi * exp(theta)
  arma::vec dp1_dphi(n);
  for (int i = 0; i < n; ++i) {
    if (same_cpp(dp0_dphi[i], 0.0)) {
      dp1_dphi[i] = 0.0;
    } else {
      dp1_dphi[i] = dp0_dphi[i] * std::exp(theta[i]);
    }
  }

  // Compute gradient
  arma::vec inner_beta(n);
  for (int i = 0; i < n; ++i) {
    double p0_i = std::fmax(clipping, std::fmin(1.0 - clipping, p0_vec[i]));
    double p1_i = std::fmax(clipping, std::fmin(1.0 - clipping, p1_vec[i]));
    double x_i = x_indicator[i];
    double y_i = y_outcome[i];

    double p_i = (x_i == 0) ? p0_i : p1_i;

    // Check if clipped
    bool is_clipped = (p_i <= clipping) || (p_i >= (1.0 - clipping));

    if (is_clipped) {
      inner_beta[i] = 0.0;
    } else {
      // dllh/dp1 * dp1/dphi + dllh/dp0 * dp0/dphi
      double dllh_dp1 = (y_i * x_i) / p1_i - ((1.0 - y_i) * x_i) / (1.0 - p1_i);
      double dllh_dp0 = (y_i * (1.0 - x_i)) / p0_i -
                        ((1.0 - y_i) * (1.0 - x_i)) / (1.0 - p0_i);
      inner_beta[i] = dllh_dp1 * dp1_dphi[i] + dllh_dp0 * dp0_dphi[i];
    }
  }

  arma::vec grad_beta = -(vb.t() * inner_beta) / static_cast<double>(n);

  return grad_beta;
}

// --- Newton-CD Optimizer in C++ ---
//' @export
// [[Rcpp::export]]
Rcpp::List newton_cd_cpp(Rcpp::NumericVector alpha_start_rcpp,
                         Rcpp::NumericVector beta_start_rcpp, double lambda,
                         bool intercept, int max_iter,
                         Rcpp::NumericMatrix va_rcpp,
                         Rcpp::NumericMatrix vb_rcpp,
                         Rcpp::NumericVector x_indicator,
                         Rcpp::NumericVector y_outcome, int prob_fun_selector,
                         double lambda_beta = -1.0, double tol = 1e-5,
                         double clipping = 1e-10, bool save_history = false) {

  // Convert to Armadillo
  arma::vec alpha = Rcpp::as<arma::vec>(alpha_start_rcpp);
  arma::vec beta = Rcpp::as<arma::vec>(beta_start_rcpp);
  arma::mat va = Rcpp::as<arma::mat>(va_rcpp);
  arma::mat vb = Rcpp::as<arma::mat>(vb_rcpp);

  int pa = alpha.n_elem;
  int pb = beta.n_elem;
  int n = y_outcome.size();

  if (lambda_beta < 0)
    lambda_beta = lambda;

  int actual_iter = 0;
  bool converged = false;

  arma::mat alphas_hist;
  arma::mat betas_hist;
  arma::vec nll_hist;

  if (save_history) {
    alphas_hist.zeros(max_iter, pa);
    betas_hist.zeros(max_iter, pb);
    nll_hist.zeros(max_iter);
  }

  // Main Newton Loop
  for (int iter = 0; iter < max_iter; ++iter) {
    actual_iter = iter + 1;
    Rcpp::checkUserInterrupt();

    double obj_prev =
        penalized_nllh_cpp(alpha, beta, va, vb, x_indicator, y_outcome, lambda,
                           intercept, prob_fun_selector, clipping);

    // 1. Compute Analytical Gradient
    arma::vec g_alpha =
        grad_nll_alpha_analytical_cpp(alpha, beta, va, vb, x_indicator,
                                      y_outcome, prob_fun_selector, clipping);
    arma::vec g_beta =
        grad_nll_beta_analytical_cpp(alpha, beta, va, vb, x_indicator,
                                     y_outcome, prob_fun_selector, clipping);

    // Clean gradients
    g_alpha.elem(arma::find_nonfinite(g_alpha)).zeros();
    g_beta.elem(arma::find_nonfinite(g_beta)).zeros();

    // 2. Compute Weights for Diagonal Hessian
    arma::vec theta = va * alpha;
    arma::vec phi = vb * beta;

    Rcpp::NumericVector theta_rcpp = Rcpp::wrap(theta);
    Rcpp::NumericVector phi_rcpp = Rcpp::wrap(phi);

    Rcpp::List ps;
    if (prob_fun_selector == 0) {
      ps = getProbRR_org_cpp(theta_rcpp, phi_rcpp, clipping);
    } else {
      ps = getProbRR_alt_cpp(theta_rcpp, phi_rcpp, clipping);
    }

    Rcpp::NumericVector p0_vec = Rcpp::as<Rcpp::NumericVector>(ps["p0"]);
    Rcpp::NumericVector p1_vec = Rcpp::as<Rcpp::NumericVector>(ps["p1"]);

    // Compute weights
    arma::vec weights(n);
    for (int i = 0; i < n; ++i) {
      double p_i = (x_indicator[i] == 0) ? p0_vec[i] : p1_vec[i];
      p_i = std::fmax(clipping, std::fmin(1.0 - clipping, p_i));
      weights[i] = std::fmax(p_i * (1.0 - p_i), 1e-4);
    }

    // Diagonal Hessian
    arma::vec h_alpha(pa);
    arma::vec h_beta(pb);

    for (int j = 0; j < pa; ++j) {
      h_alpha[j] =
          arma::dot(va.col(j) % va.col(j), weights) / static_cast<double>(n) +
          1e-6;
    }
    for (int j = 0; j < pb; ++j) {
      h_beta[j] =
          arma::dot(vb.col(j) % vb.col(j), weights) / static_cast<double>(n) +
          1e-6;
    }

    // 3. Coordinate Descent Update
    arma::vec z_alpha = alpha % h_alpha - g_alpha;
    arma::vec z_beta = beta % h_beta - g_beta;

    // Soft thresholding
    arma::vec alpha_new(pa);
    arma::vec beta_new(pb);

    for (int j = 0; j < pa; ++j) {
      double lam_j = (intercept && j == 0) ? 0.0 : lambda;
      Rcpp::NumericVector z_j = Rcpp::NumericVector::create(z_alpha[j]);
      Rcpp::NumericVector thresh_j = soft_thres_cpp(z_j, lam_j);
      alpha_new[j] = thresh_j[0] / h_alpha[j];
    }

    for (int j = 0; j < pb; ++j) {
      double lam_j = (intercept && j == 0) ? 0.0 : lambda_beta;
      Rcpp::NumericVector z_j = Rcpp::NumericVector::create(z_beta[j]);
      Rcpp::NumericVector thresh_j = soft_thres_cpp(z_j, lam_j);
      beta_new[j] = thresh_j[0] / h_beta[j];
    }

    // 4. Line Search (Backtracking)
    arma::vec d_alpha = alpha_new - alpha;
    arma::vec d_beta = beta_new - beta;

    // Line search

    double step_ls = 1.0;
    bool accepted = false;
    for (int ls = 0; ls < 30; ++ls) {
      arma::vec a_cand = alpha + step_ls * d_alpha;
      arma::vec b_cand = beta + step_ls * d_beta;
      double obj_cand =
          penalized_nllh_cpp(a_cand, b_cand, va, vb, x_indicator, y_outcome,
                             lambda, intercept, prob_fun_selector, clipping);

      if (obj_cand <= obj_prev && std::isfinite(obj_cand)) {
        alpha = a_cand;
        beta = b_cand;
        accepted = true;
        if (save_history) {
          nll_hist(iter) = obj_cand;
          alphas_hist.row(iter) = alpha.t();
          betas_hist.row(iter) = beta.t();
        }

        // Update obj_prev for next iteration's check!
        obj_prev = obj_cand;

        // Check convergence on the ACTUAL step taken
        double step_norm_a =
            (d_alpha.n_elem > 0) ? arma::abs(d_alpha).max() * step_ls : 0.0;
        double step_norm_b =
            (d_beta.n_elem > 0) ? arma::abs(d_beta).max() * step_ls : 0.0;
        if (std::max(step_norm_a, step_norm_b) < tol) {
          converged = true;
        }
        break;
      }
      step_ls *= 0.5;
    }

    if (converged)
      break;

    if (!accepted) {
      if (save_history) {
        alphas_hist.row(iter) = alpha.t();
        betas_hist.row(iter) = beta.t();
        nll_hist(iter) = obj_prev;
      }
      break;
    }

    if (save_history) {
      alphas_hist.row(iter) = alpha.t();
      betas_hist.row(iter) = beta.t();
    }
  }

  double final_nll =
      penalized_nllh_cpp(alpha, beta, va, vb, x_indicator, y_outcome, lambda,
                         intercept, prob_fun_selector, clipping);

  // Return results (compatible with R interface)
  return Rcpp::List::create(
      Rcpp::Named("alpha") = Rcpp::wrap(alpha),
      Rcpp::Named("beta") = Rcpp::wrap(beta),
      Rcpp::Named("convergence") = converged, Rcpp::Named("step") = actual_iter,
      Rcpp::Named("final_nll") = final_nll,
      Rcpp::Named("alphas") =
          (save_history && actual_iter > 0)
              ? Rcpp::wrap(alphas_hist.rows(0, actual_iter - 1))
              : R_NilValue,
      Rcpp::Named("betas") =
          (save_history && actual_iter > 0)
              ? Rcpp::wrap(betas_hist.rows(0, actual_iter - 1))
              : R_NilValue,
      Rcpp::Named("nllh_results") = (save_history && actual_iter > 0)
                                        ? Rcpp::wrap(nll_hist.head(actual_iter))
                                        : R_NilValue);
}

//' @export
// [[Rcpp::export]]
Rcpp::List active_set_newton_cd_cpp(
    Rcpp::NumericVector alpha_start_rcpp, Rcpp::NumericVector beta_start_rcpp,
    double lambda, bool intercept, int max_iter, Rcpp::NumericMatrix va_rcpp,
    Rcpp::NumericMatrix vb_rcpp, Rcpp::NumericVector x_indicator,
    Rcpp::NumericVector y_outcome, int prob_fun_selector,
    double lambda_beta = -1.0, double tol = 1e-5, double clipping = 1e-10,
    int kkt_check_freq = 10, double active_tol = 1e-6,
    bool save_history = false) {

  // Convert to Armadillo
  arma::vec alpha = Rcpp::as<arma::vec>(alpha_start_rcpp);
  arma::vec beta = Rcpp::as<arma::vec>(beta_start_rcpp);
  arma::mat va = Rcpp::as<arma::mat>(va_rcpp);
  arma::mat vb = Rcpp::as<arma::mat>(vb_rcpp);

  int pa = alpha.n_elem;
  int pb = beta.n_elem;
  int n = y_outcome.size();

  if (lambda_beta < 0)
    lambda_beta = lambda;

  // Active sets (bools)
  arma::uvec active_alpha(pa, arma::fill::ones); // All active initially
  arma::uvec active_beta(pb, arma::fill::ones);

  int actual_iter = 0;
  bool converged = false;

  arma::mat alphas_hist;
  arma::mat betas_hist;
  arma::vec nll_hist;

  if (save_history) {
    alphas_hist.zeros(max_iter, pa);
    betas_hist.zeros(max_iter, pb);
    nll_hist.zeros(max_iter);
  }

  for (int iter = 0; iter < max_iter; ++iter) {
    actual_iter = iter + 1;
    Rcpp::checkUserInterrupt();

    double obj_prev =
        penalized_nllh_cpp(alpha, beta, va, vb, x_indicator, y_outcome, lambda,
                           intercept, prob_fun_selector, clipping);

    // KKT Check logic
    bool check_kkt = (iter == 0) || ((iter + 1) % kkt_check_freq == 0);
    // Note: R implementation checks at iter=1 (first step) which means
    // effectively it computes full gradient first. Here iter 0 is first.

    // Calculate Theta and Phi (Efficiently using active set)
    // Inactive alpha/beta are 0. So full va*alpha is same as
    // va_act * alpha_act.
    arma::uvec alpha_idx = arma::find(active_alpha);
    arma::uvec beta_idx = arma::find(active_beta);

    arma::vec theta(n, arma::fill::zeros);
    if (alpha_idx.n_elem > 0) {
      theta = va.cols(alpha_idx) * alpha.elem(alpha_idx);
    }

    arma::vec phi(n, arma::fill::zeros);
    if (beta_idx.n_elem > 0) {
      phi = vb.cols(beta_idx) * beta.elem(beta_idx);
    }

    // Compute Probabilities and Weights (O(n))
    Rcpp::NumericVector theta_rcpp = Rcpp::wrap(theta);
    Rcpp::NumericVector phi_rcpp = Rcpp::wrap(phi);

    Rcpp::List ps;
    if (prob_fun_selector == 0) {
      ps = getProbRR_org_cpp(theta_rcpp, phi_rcpp, clipping);
    } else {
      ps = getProbRR_alt_cpp(theta_rcpp, phi_rcpp, clipping);
    }

    Rcpp::NumericVector p0_vec = Rcpp::as<Rcpp::NumericVector>(ps["p0"]);
    Rcpp::NumericVector p1_vec = Rcpp::as<Rcpp::NumericVector>(ps["p1"]);

    arma::vec weights(n);
    arma::vec inner_alpha(n);
    arma::vec inner_beta(n);

    // This loop computes weights AND inner vector elements for gradient reuse
    // Gradient calculation:
    // grad = -X' * inner
    // inner_alpha_i = dllh/dp * dp/dtheta
    // inner_beta_i = dllh/dp * dp/dphi
    // Reuse logic from analytical gradient functions, but optimized loop?
    // We can't easily reuse the exported functions because they take full
    // alpha/beta and recompute everything. We want to reuse `theta`, `phi`,
    // `ps`.

    // We need dp0_dtheta/dphi.
    arma::vec dp0_dtheta;
    arma::vec dp0_dphi;
    if (prob_fun_selector == 0) {
      dp0_dtheta = dp0_theta_richardson_cpp(theta, phi);
      dp0_dphi = dp0_phi_richardson_cpp(theta, phi);
    } else {
      dp0_dtheta = dp0_theta_pozza_cpp(theta, phi);
      dp0_dphi = dp0_phi_pozza_cpp(theta, phi);
    }

    arma::vec dp1_dtheta(n);
    arma::vec dp1_dphi(n);

    for (int i = 0; i < n; ++i) {

      // Let's copy exact logic from analytical_gradient
      double p0_i = p0_vec[i]; // raw
      if (same_cpp(p0_i + dp0_dtheta[i], 0.0)) {
        dp1_dtheta[i] = 0.0;
      } else {
        dp1_dtheta[i] = (p0_i + dp0_dtheta[i]) * std::exp(theta[i]);
      }

      // dp1/dphi
      if (same_cpp(dp0_dphi[i], 0.0)) {
        dp1_dphi[i] = 0.0;
      } else {
        dp1_dphi[i] = dp0_dphi[i] * std::exp(theta[i]);
      }

      // Weights
      double p_i = (x_indicator[i] == 0) ? p0_i : p1_vec[i];
      p_i = std::fmax(clipping, std::fmin(1.0 - clipping, p_i));
      weights[i] = std::fmax(p_i * (1.0 - p_i), 1e-4);

      // Inner terms
      double p0_c = std::fmax(clipping, std::fmin(1.0 - clipping, p0_vec[i]));
      double p1_c = std::fmax(clipping, std::fmin(1.0 - clipping, p1_vec[i]));
      double yi = y_outcome[i];
      double xi = x_indicator[i];

      bool is_clipped_p = (xi == 0)
                              ? (p0_c <= clipping || p0_c >= 1.0 - clipping)
                              : (p1_c <= clipping || p1_c >= 1.0 - clipping);

      if (is_clipped_p &&
          false) { // R doesn't seem to zero gradient on clipped?
        // R grad_nll: clipping affects P used in dllh/dp terms.
        // My C++ grad uses explicit `is_clipped` check to zero `inner`.
        // Let's stick to C++ grad logic.
        inner_alpha[i] = 0.0;
        inner_beta[i] = 0.0;
      } else {
        double dllh_dp1 = (yi * xi) / p1_c - ((1.0 - yi) * xi) / (1.0 - p1_c);
        double dllh_dp0 =
            (yi * (1.0 - xi)) / p0_c - ((1.0 - yi) * (1.0 - xi)) / (1.0 - p0_c);
        inner_alpha[i] = dllh_dp1 * dp1_dtheta[i] + dllh_dp0 * dp0_dtheta[i];
        inner_beta[i] = dllh_dp1 * dp1_dphi[i] + dllh_dp0 * dp0_dphi[i];
      }
    }

    // Now compute Gradient and Hessian (Subset or Full)
    arma::vec g_alpha;
    arma::vec h_alpha;
    arma::vec g_beta;
    arma::vec h_beta;

    // We need full containers if check_kkt is true, but we can just compute
    // subsets + violations To minimize memory, let's keep g as sparse or full?
    // Armadillo vec is dense.

    // Alpha
    g_alpha = arma::vec(pa, arma::fill::zeros);
    h_alpha = arma::vec(pa, arma::fill::zeros);

    if (check_kkt) {
      // Compute full G and H
      if (pa > 0) {
        g_alpha = -(va.t() * inner_alpha) / static_cast<double>(n);
        // Full H diagonal
        // h[j] = dot(v_j^2, w)
        // This is va.each_col() op?
        // arma::sum(square(va).each_col() % weights) ?
        // Or loop. O(np).
        for (int j = 0; j < pa; ++j) {
          h_alpha[j] = arma::dot(va.col(j) % va.col(j), weights) / n + 1e-6;
        }
      }
    } else {
      // Compute subset
      if (alpha_idx.n_elem > 0) {
        arma::mat va_sub = va.cols(alpha_idx);
        arma::vec g_sub = -(va_sub.t() * inner_alpha) / n;
        g_alpha.elem(alpha_idx) = g_sub;

        for (unsigned int k = 0; k < alpha_idx.n_elem; ++k) {
          int j = alpha_idx[k];
          h_alpha[j] = arma::dot(va.col(j) % va.col(j), weights) / n + 1e-6;
        }
      }
    }

    // Beta
    g_beta = arma::vec(pb, arma::fill::zeros);
    h_beta = arma::vec(pb, arma::fill::zeros);

    if (check_kkt) {
      if (pb > 0) {
        g_beta = -(vb.t() * inner_beta) / static_cast<double>(n);
        for (int j = 0; j < pb; ++j) {
          h_beta[j] = arma::dot(vb.col(j) % vb.col(j), weights) / n + 1e-6;
        }
      }
    } else {
      if (beta_idx.n_elem > 0) {
        arma::mat vb_sub = vb.cols(beta_idx);
        arma::vec g_sub = -(vb_sub.t() * inner_beta) / n;
        g_beta.elem(beta_idx) = g_sub;

        for (unsigned int k = 0; k < beta_idx.n_elem; ++k) {
          int j = beta_idx[k];
          h_beta[j] = arma::dot(vb.col(j) % vb.col(j), weights) / n + 1e-6;
        }
      }
    }
    // Clean gradients
    g_alpha.elem(arma::find_nonfinite(g_alpha)).zeros();
    g_beta.elem(arma::find_nonfinite(g_beta)).zeros();

    // KKT Check (if check_kkt)
    if (check_kkt) {
      if (lambda > 0) {
        arma::uvec inactive_idx = arma::find(active_alpha == 0);
        if (inactive_idx.n_elem > 0) {
          // z = alpha * h - g. But alpha=0. so z = -g.
          // violated: abs(-g) > lambda + tol
          for (arma::uword k = 0; k < inactive_idx.n_elem; ++k) {
            int j = inactive_idx[k];
            // Skip intercept KKT check (intercept always active usually, but if
            // not?)
            if (intercept && j == 0)
              continue;
            if (std::abs(-g_alpha[j]) > lambda + active_tol) {
              active_alpha[j] = 1;
              // Update H for this new active var (needed for update step)
              h_alpha[j] = arma::dot(va.col(j) % va.col(j), weights) / n + 1e-6;
            }
          }
        }
      }
      if (lambda_beta > 0) {
        arma::uvec inactive_idx = arma::find(active_beta == 0);
        if (inactive_idx.n_elem > 0) {
          for (arma::uword k = 0; k < inactive_idx.n_elem; ++k) {
            int j = inactive_idx[k];
            if (intercept && j == 0)
              continue;
            if (std::abs(-g_beta[j]) > lambda_beta + active_tol) {
              active_beta[j] = 1;
              h_beta[j] = arma::dot(vb.col(j) % vb.col(j), weights) / n + 1e-6;
            }
          }
        }
      }
    }

    // Coordinate Descent Update (Only Active)
    arma::vec alpha_new = alpha;
    arma::vec beta_new = beta;

    // Refresh active indices after KKT
    alpha_idx = arma::find(active_alpha);
    beta_idx = arma::find(active_beta);

    // Alpha Update
    for (arma::uword k = 0; k < alpha_idx.n_elem; ++k) {
      int j = alpha_idx[k];
      double lam_j = (intercept && j == 0) ? 0.0 : lambda;
      // z = alpha*h - g
      double z_val = alpha[j] * h_alpha[j] - g_alpha[j];

      // Soft thres
      double thresh_val = 0.0;
      double abs_z = std::fabs(z_val);
      if (lam_j >= 0) {
        if (z_val > 0)
          thresh_val = std::fmax(0.0, abs_z - lam_j) * 1.0;
        else if (z_val < 0)
          thresh_val = std::fmax(0.0, abs_z - lam_j) * -1.0;
      } else { // Should not happen
        thresh_val = z_val;
      }

      alpha_new[j] = thresh_val / h_alpha[j];
    }

    // Beta Update
    for (arma::uword k = 0; k < beta_idx.n_elem; ++k) {
      int j = beta_idx[k];
      double lam_j = (intercept && j == 0) ? 0.0 : lambda_beta;
      double z_val = beta[j] * h_beta[j] - g_beta[j];

      double thresh_val = 0.0;
      double abs_z = std::fabs(z_val);
      if (lam_j >= 0) {
        if (z_val > 0)
          thresh_val = std::fmax(0.0, abs_z - lam_j) * 1.0;
        else if (z_val < 0)
          thresh_val = std::fmax(0.0, abs_z - lam_j) * -1.0;
      }
      beta_new[j] = thresh_val / h_beta[j];
    }

    // Shrink inactive to 0 (already 0, but safe)
    // Actually active vars might become 0.

    // Line Search
    arma::vec d_alpha = alpha_new - alpha;
    arma::vec d_beta = beta_new - beta;

    double step_ls = 1.0;
    bool accepted = false;
    for (int ls = 0; ls < 30; ++ls) {
      arma::vec a_cand = alpha + step_ls * d_alpha;
      arma::vec b_cand = beta + step_ls * d_beta;
      double obj_cand =
          penalized_nllh_cpp(a_cand, b_cand, va, vb, x_indicator, y_outcome,
                             lambda, intercept, prob_fun_selector, clipping);

      if (obj_cand <= obj_prev && std::isfinite(obj_cand)) {
        alpha = a_cand;
        beta = b_cand;
        accepted = true;
        if (save_history) {
          nll_hist(iter) = obj_cand;
          alphas_hist.row(iter) = alpha.t();
          betas_hist.row(iter) = beta.t();
        }

        // Update obj_prev for next iteration
        obj_prev = obj_cand;

        // Check convergence on the ACTUAL step taken
        double step_norm_a =
            (d_alpha.n_elem > 0) ? arma::abs(d_alpha).max() * step_ls : 0.0;
        double step_norm_b =
            (d_beta.n_elem > 0) ? arma::abs(d_beta).max() * step_ls : 0.0;
        if (std::max(step_norm_a, step_norm_b) < tol) {
          converged = true;
        }
        break;
      }
      step_ls *= 0.5;
    }

    if (converged)
      break;

    if (!accepted) {
      if (save_history) {
        alphas_hist.row(iter) = alpha.t();
        betas_hist.row(iter) = beta.t();
        nll_hist(iter) = obj_prev;
      }
      break;
    }

    // Update Active Sets (Shrink)
    // Remove if abs <= active_tol
    // Intercept always active
    if (lambda > 0) {
      for (int j = 0; j < pa; ++j) {
        if (active_alpha[j]) {
          if (std::abs(alpha[j]) <= active_tol && (!intercept || j != 0)) {
            active_alpha[j] = 0;
          }
        }
      }
    }
    if (lambda_beta > 0) {
      for (int j = 0; j < pb; ++j) {
        if (active_beta[j]) {
          if (std::abs(beta[j]) <= active_tol && (!intercept || j != 0)) {
            active_beta[j] = 0;
          }
        }
      }
    }
  }

  double final_nll =
      penalized_nllh_cpp(alpha, beta, va, vb, x_indicator, y_outcome, lambda,
                         intercept, prob_fun_selector, clipping);

  return Rcpp::List::create(
      Rcpp::Named("alpha") = Rcpp::wrap(alpha),
      Rcpp::Named("beta") = Rcpp::wrap(beta),
      Rcpp::Named("convergence") = converged, Rcpp::Named("step") = actual_iter,
      Rcpp::Named("final_nll") = final_nll,
      Rcpp::Named("alphas") =
          (save_history && actual_iter > 0)
              ? Rcpp::wrap(alphas_hist.rows(0, actual_iter - 1))
              : R_NilValue,
      Rcpp::Named("betas") =
          (save_history && actual_iter > 0)
              ? Rcpp::wrap(betas_hist.rows(0, actual_iter - 1))
              : R_NilValue,
      Rcpp::Named("nllh_results") = (save_history && actual_iter > 0)
                                        ? Rcpp::wrap(nll_hist.head(actual_iter))
                                        : R_NilValue);
}

// --- FISTA Optimizer in C++ ---

// ' @export
// [[Rcpp::export]]
Rcpp::List fista_cpp(Rcpp::NumericVector alpha_start_rcpp,
                     Rcpp::NumericVector beta_start_rcpp, double lambda,
                     bool intercept, int max_iter, Rcpp::NumericMatrix va_rcpp,
                     Rcpp::NumericMatrix vb_rcpp,
                     Rcpp::NumericVector x_indicator,
                     Rcpp::NumericVector y_outcome, int prob_fun_selector,
                     double lambda_beta = -1.0, double tol = 1e-5,
                     double step_size_init = 0.5, double armijo_c = 1e-4,
                     double shrink_factor = 0.5, double clipping = 1e-10,
                     bool save_history = false) {

  // Convert Rcpp objects to Armadillo
  arma::vec alpha = Rcpp::as<arma::vec>(alpha_start_rcpp);
  arma::vec beta = Rcpp::as<arma::vec>(beta_start_rcpp);
  arma::mat va = Rcpp::as<arma::mat>(va_rcpp);
  arma::mat vb = Rcpp::as<arma::mat>(vb_rcpp);

  // Handle lambda_beta
  if (lambda_beta < 0)
    lambda_beta = lambda;

  // Helper for full objective (handling lambda_beta correctly)
  auto calc_obj = [&](const arma::vec &a, const arma::vec &b) {
    double val = nllh_cpp(a, b, va, vb, x_indicator, y_outcome,
                          prob_fun_selector, clipping);
    double pen = 0.0;
    int s_a = (intercept && a.n_elem > 0) ? 1 : 0;
    for (arma::uword i = s_a; i < a.n_elem; ++i)
      pen += std::abs(a[i]) * lambda;
    int s_b = (intercept && b.n_elem > 0) ? 1 : 0;
    for (arma::uword i = s_b; i < b.n_elem; ++i)
      pen += std::abs(b[i]) * lambda_beta;
    return val + pen;
  };

  // Initialize FISTA variables
  arma::vec alpha_y = alpha;
  arma::vec beta_y = beta;
  arma::vec alpha_old = alpha; // x_k
  arma::vec beta_old = beta;

  double t_k = 1.0;
  double t_k_next;

  // Initial step size
  double step_size = step_size_init;

  int final_iter = 0;
  bool converged = false;

  arma::mat alphas_hist;
  arma::mat betas_hist;
  arma::vec nll_hist;

  if (save_history) {
    alphas_hist.set_size(max_iter, alpha.n_elem);
    betas_hist.set_size(max_iter, beta.n_elem);
    nll_hist.set_size(max_iter);
  }

  double nll_current = calc_obj(alpha, beta);

  for (int k = 0; k < max_iter; ++k) {
    // 1. Compute Gradient at y_k
    arma::vec grad_alpha =
        grad_nll_alpha_analytical_cpp(alpha_y, beta_y, va, vb, x_indicator,
                                      y_outcome, prob_fun_selector, clipping);
    arma::vec grad_beta =
        grad_nll_beta_analytical_cpp(alpha_y, beta_y, va, vb, x_indicator,
                                     y_outcome, prob_fun_selector, clipping);

    // 2. Backtracking Line Search
    double nll_y = nllh_cpp(alpha_y, beta_y, va, vb, x_indicator, y_outcome,
                            prob_fun_selector, clipping);

    // Fallback if nll_y is not finite (should rarely happen if previous steps
    // were safe)
    if (!R_finite(nll_y)) {
      // Can't trust gradient if point is invalid.
      // Usually means step size too big previously or momentum pushed out.
      // Reset to x_k (no momentum) and retry?
      alpha_y = alpha;
      beta_y = beta;
      t_k = 1.0;
      // Recompute gradient at x_k
      grad_alpha =
          grad_nll_alpha_analytical_cpp(alpha_y, beta_y, va, vb, x_indicator,
                                        y_outcome, prob_fun_selector, clipping);
      grad_beta =
          grad_nll_beta_analytical_cpp(alpha_y, beta_y, va, vb, x_indicator,
                                       y_outcome, prob_fun_selector, clipping);
      nll_y = nllh_cpp(alpha_y, beta_y, va, vb, x_indicator, y_outcome,
                       prob_fun_selector, clipping);
      nll_y = nllh_cpp(alpha_y, beta_y, va, vb, x_indicator, y_outcome,
                       prob_fun_selector, clipping);
    }

    // Monotone FISTA: Ensure F(x_{k+1}) <= F(x_k)
    // We generated candidate x_{k+1} (alpha_new) from y_k.
    // We will check its NLL later.

    arma::vec alpha_new, beta_new;
    double current_step = step_size;

    // Try to increase step size slightly from previous success (adaptive)?
    // Standard approach: start with previous step size.
    // Some implementations try current_step * 1.2 first. Let's stick to safe
    // start.

    double nll_new = nll_y; // Default if loop doesn't run or first iter check
                            // fails? No, loop updates it.
    for (int j = 0; j < 50; ++j) { // Limit backtracks
      // Gradient Step
      arma::vec z_alpha = alpha_y - current_step * grad_alpha;
      arma::vec z_beta = beta_y - current_step * grad_beta;

      // Proximal Operator (Soft Thresholding)
      alpha_new = z_alpha; // Copies structure
      int start_idx_a = (intercept && alpha.n_elem > 0) ? 1 : 0;
      for (arma::uword i = start_idx_a; i < alpha.n_elem; ++i) {
        double val = z_alpha[i];
        double thresh = lambda * current_step;
        if (val > thresh)
          alpha_new[i] = val - thresh;
        else if (val < -thresh)
          alpha_new[i] = val + thresh;
        else
          alpha_new[i] = 0.0;
      }

      beta_new = z_beta;
      int start_idx_b = (intercept && beta.n_elem > 0) ? 1 : 0;
      for (arma::uword i = start_idx_b; i < beta.n_elem; ++i) {
        double val = z_beta[i];
        double thresh = lambda_beta * current_step;
        if (val > thresh)
          beta_new[i] = val - thresh;
        else if (val < -thresh)
          beta_new[i] = val + thresh;
        else
          beta_new[i] = 0.0;
      }

      // Check Armijo Condition on NLL (Smooth part)
      nll_new = nllh_cpp(alpha_new, beta_new, va, vb, x_indicator, y_outcome,
                         prob_fun_selector, clipping);

      if (!R_finite(nll_new)) {
        current_step *= shrink_factor;
        continue;
      }

      arma::vec diff_a = alpha_new - alpha_y;
      arma::vec diff_b = beta_new - beta_y;
      double sq_norm_diff =
          arma::dot(diff_a, diff_a) + arma::dot(diff_b, diff_b);
      double dot_grad =
          arma::dot(grad_alpha, diff_a) + arma::dot(grad_beta, diff_b);

      // Q(x, y) = f(y) + <g, x-y> + 1/2t ||x-y||^2
      double rhs =
          nll_y + dot_grad + (1.0 / (2.0 * current_step)) * sq_norm_diff;

      if (nll_new <= rhs + 1e-10) {
        break;
      }

      current_step *= shrink_factor;
    }

    step_size = current_step; // Save for next iter

    // Check convergence on x_k (alpha vs alpha_new) or x_{k+1} - x_k?
    // Usually ||x_{k+1} - x_k|| / max(1, ||x_k||) or similar.
    // Using absolute change as consistent with Newton-CD check.

    // alpha_new is x_{k+1}. alpha is x_k.
    arma::vec change_a = alpha_new - alpha;
    arma::vec change_b = beta_new - beta;
    double max_change =
        std::max(arma::abs(change_a).max(), arma::abs(change_b).max());

    if (max_change < tol) {
      beta = beta_new;
      converged = true;
      final_iter = k + 1;
      if (save_history) {
        alphas_hist.row(k) = alpha.t();
        betas_hist.row(k) = beta.t();
        nll_hist(k) = nll_new;
      }
      break;
    }

    // --- Monotonicity Check (Restart) ---
    if (nll_new > nll_current) {
      // Reject momentum step. Restart to x_k.
      alpha_y = alpha;
      beta_y = beta;
      t_k = 1.0;

      // Recompute gradients at x_k
      arma::vec g_a_k =
          grad_nll_alpha_analytical_cpp(alpha, beta, va, vb, x_indicator,
                                        y_outcome, prob_fun_selector, clipping);
      arma::vec g_b_k =
          grad_nll_beta_analytical_cpp(alpha, beta, va, vb, x_indicator,
                                       y_outcome, prob_fun_selector, clipping);

      double nll_x = nll_current;
      double cur_step = step_size;

      // Simple backtracking GD from x_k
      bool step_found = false;
      arma::vec a_gd, b_gd;
      double nll_gd = nll_new;

      for (int r = 0; r < 20; ++r) {
        arma::vec z_a = alpha - cur_step * g_a_k;
        arma::vec z_b = beta - cur_step * g_b_k;

        a_gd = z_a;
        int s_a = (intercept && alpha.n_elem > 0) ? 1 : 0;
        for (arma::uword i = s_a; i < alpha.n_elem; ++i) {
          double v = z_a[i];
          double th = lambda * cur_step;
          if (v > th)
            a_gd[i] = v - th;
          else if (v < -th)
            a_gd[i] = v + th;
          else
            a_gd[i] = 0.0;
        }

        b_gd = z_b;
        int s_b = (intercept && beta.n_elem > 0) ? 1 : 0;
        for (arma::uword i = s_b; i < beta.n_elem; ++i) {
          double v = z_b[i];
          double th = lambda_beta * cur_step;
          if (v > th)
            b_gd[i] = v - th;
          else if (v < -th)
            b_gd[i] = v + th;
          else
            b_gd[i] = 0.0;
        }

        nll_gd =
            penalized_nllh_cpp(a_gd, b_gd, va, vb, x_indicator, y_outcome,
                               lambda, intercept, prob_fun_selector, clipping);

        // Check sufficient decrease or simple monotonicity
        if (nll_gd <= nll_x) {
          step_found = true;
          break;
        }
        cur_step *= shrink_factor;
      }

      if (step_found) {
        alpha_new = a_gd;
        beta_new = b_gd;
        nll_new = nll_gd;
        step_size = cur_step;
      } else {
        alpha_new = alpha;
        beta_new = beta;
        nll_new = nll_current;
      }
    }

    nll_current = nll_new;

    // Update Momentum
    t_k_next = (1.0 + std::sqrt(1.0 + 4.0 * t_k * t_k)) / 2.0;
    double momentum = (t_k - 1.0) / t_k_next;

    // Restart logic? If gradient direction opposes momentum.
    // Optional. For now standard FISTA.

    alpha_y = alpha_new + momentum * (alpha_new - alpha);
    beta_y = beta_new + momentum * (beta_new - beta);

    // Update x_k
    alpha = alpha_new;
    beta = beta_new;

    t_k = t_k_next;
    final_iter = k + 1;

    if (save_history) {
      alphas_hist.row(k) = alpha.t();
      betas_hist.row(k) = beta.t();
      // nll_hist(k) = nll_new; // Need to access nll_new, but it is local to
      // loop. nll_new is computed in line search. I need to surface it or
      // recompute. Recomputing is redundant. Assuming line search accepted,
      // nll_new is valid. But nll_new is strictly local. I should bubble it
      // out. For now, recomputing or storing last valid nll_y.
      nll_hist(k) = nll_current;
    }
  }

  double final_nll =
      penalized_nllh_cpp(alpha, beta, va, vb, x_indicator, y_outcome, lambda,
                         intercept, prob_fun_selector, clipping);

  return Rcpp::List::create(
      Rcpp::Named("alpha") = Rcpp::wrap(alpha),
      Rcpp::Named("beta") = Rcpp::wrap(beta),
      Rcpp::Named("convergence") = converged, Rcpp::Named("step") = final_iter,
      Rcpp::Named("final_nll") = final_nll,
      Rcpp::Named("alphas") =
          (save_history && final_iter > 0)
              ? Rcpp::wrap(alphas_hist.rows(0, final_iter - 1))
              : R_NilValue,
      Rcpp::Named("betas") =
          (save_history && final_iter > 0)
              ? Rcpp::wrap(betas_hist.rows(0, final_iter - 1))
              : R_NilValue,
      Rcpp::Named("nllh_results") = (save_history && final_iter > 0)
                                        ? Rcpp::wrap(nll_hist.head(final_iter))
                                        : R_NilValue);
}

// --- L-BFGS Optimizer (Wrapper for R API) ---

struct LBFGS_Data {
  const arma::mat &va;
  const arma::mat &vb;
  const Rcpp::NumericVector &x_indicator;
  const Rcpp::NumericVector &y_outcome;
  double lambda;
  double lambda_beta;
  bool intercept;
  int prob_fun_selector;
  double clipping;
  int pa;
  int pb;
};

// Objective function wrapper
double lbfgs_obj(int n, double *par, void *ex) {
  LBFGS_Data *data = (LBFGS_Data *)ex;

  // Map existing memory to arma vectors (no copy)
  arma::vec alpha(par, data->pa, false, true);
  arma::vec beta(par + data->pa, data->pb, false, true);

  return penalized_nllh_cpp(alpha, beta, data->va, data->vb, data->x_indicator,
                            data->y_outcome, data->lambda, data->intercept,
                            data->prob_fun_selector, data->clipping);
}

// Gradient function wrapper
void lbfgs_grad(int n, double *par, double *gr, void *ex) {
  LBFGS_Data *data = (LBFGS_Data *)ex;

  arma::vec alpha(par, data->pa, false, true);
  arma::vec beta(par + data->pa, data->pb, false, true);

  arma::vec ga = grad_nll_alpha_analytical_cpp(
      alpha, beta, data->va, data->vb, data->x_indicator, data->y_outcome,
      data->prob_fun_selector, data->clipping);

  arma::vec gb = grad_nll_beta_analytical_cpp(
      alpha, beta, data->va, data->vb, data->x_indicator, data->y_outcome,
      data->prob_fun_selector, data->clipping);

  // Add L1 penalty subgradient (approximated or omitted? Standard optim adds it
  // if fn has it) penalized_nllh includes penalty. Gradients should too. L1
  // derivative is sign(x) * lambda. Note: L-BFGS on non-smooth is risky.

  // Copy to gr output
  // Alpha part
  for (int i = 0; i < data->pa; ++i) {
    double penalty = 0.0;
    // Don't penalize intercept if it exists and is first index
    if (!(data->intercept && i == 0)) {
      if (alpha[i] > 0)
        penalty = data->lambda;
      else if (alpha[i] < 0)
        penalty = -data->lambda;
    }
    gr[i] = ga[i] + penalty;
  }

  // Beta part
  for (int i = 0; i < data->pb; ++i) {
    double penalty = 0.0;
    if (!(data->intercept && i == 0)) {
      if (beta[i] > 0)
        penalty = data->lambda_beta;
      else if (beta[i] < 0)
        penalty = -data->lambda_beta;
    }
    gr[data->pa + i] = gb[i] + penalty;
  }
}

// ' @export
// [[Rcpp::export]]
Rcpp::List lbfgs_cpp(Rcpp::NumericVector alpha_start_rcpp,
                     Rcpp::NumericVector beta_start_rcpp, double lambda,
                     bool intercept, int max_iter, Rcpp::NumericMatrix va_rcpp,
                     Rcpp::NumericMatrix vb_rcpp,
                     Rcpp::NumericVector x_indicator,
                     Rcpp::NumericVector y_outcome, int prob_fun_selector,
                     double lambda_beta = -1.0, double tol = 1e-5,
                     double clipping = 1e-10, bool save_history = false) {

  arma::mat va = Rcpp::as<arma::mat>(va_rcpp);
  arma::mat vb = Rcpp::as<arma::mat>(vb_rcpp);

  int pa = alpha_start_rcpp.length();
  int pb = beta_start_rcpp.length();
  int n_params = pa + pb;

  if (lambda_beta < 0)
    lambda_beta = lambda;

  // Initial parameters
  std::vector<double> par(n_params);
  for (int i = 0; i < pa; ++i)
    par[i] = alpha_start_rcpp[i];
  for (int i = 0; i < pb; ++i)
    par[pa + i] = beta_start_rcpp[i];

  // Bounds (Unbounded for now, but lbfgsb requires them)
  // We can set +/- Inf
  std::vector<double> lower(n_params, -std::numeric_limits<double>::infinity());
  std::vector<double> upper(n_params, std::numeric_limits<double>::infinity());
  std::vector<int> nbd(n_params, 0); // 0 means unbounded

  // L-BFGS-B parameters
  int m = 5; // Memory
  int fail = 0;
  double factr = 1e7; // Default factr ~ 1e-5 relative tolerance?
  // factr is tolerance on reduction factor. 1e7 is low precision (approx 1e-8
  // machine epsilon?) R optim default is 1e7.
  double pgtol = tol; // Gradient tolerance
  int fncount = 0;
  int grcount = 0;
  int trace = 0;
  int nREPORT = 10;
  char msg[60];
  double Fmin = 0.0;

  // Data struct
  LBFGS_Data ex_data = {va,       vb,          x_indicator, y_outcome,
                        lambda,   lambda_beta, intercept,   prob_fun_selector,
                        clipping, pa,          pb};

  // Call lbfgsb
  lbfgsb(n_params, m, par.data(), lower.data(), upper.data(), nbd.data(), &Fmin,
         lbfgs_obj, lbfgs_grad, &fail, &ex_data, factr, pgtol, &fncount,
         &grcount, max_iter, msg, trace, nREPORT);

  // Unpack results
  arma::vec alpha_out(pa);
  arma::vec beta_out(pb);
  for (int i = 0; i < pa; ++i)
    alpha_out[i] = par[i];
  for (int i = 0; i < pb; ++i)
    beta_out[i] = par[pa + i];

  bool converged =
      (fail == 0); // fail=0 usually means convergence (or no error)

  return Rcpp::List::create(
      Rcpp::Named("alpha") = alpha_out, Rcpp::Named("beta") = beta_out,
      Rcpp::Named("convergence") = converged, // fail code?
      Rcpp::Named("step") =
          fncount, // or generic iterations? fncount roughly resembles
      Rcpp::Named("final_nll") = Fmin,
      Rcpp::Named("message") = std::string(msg));
}
