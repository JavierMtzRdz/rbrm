// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(fntl)]]
#include <RcppArmadillo.h>   // For Armadillo matrix/vector operations
#include <fntl.h>
#include <cmath>
#include <vector>
#include <string>
#include <limits> 

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
Rcpp::List getProbRR_org_cpp(const Rcpp::NumericVector& logrr,
                              const Rcpp::NumericVector& logop,
                              bool clipping = true) {
   int n = logrr.size();
   if (logop.size() != n) {
     Rcpp::stop("In prob_rr_org_worker: logrr and logop must have the same length.");
   }
   
   Rcpp::NumericVector p0_out(n);
   Rcpp::NumericVector p1_out(n);
   
   const double boundary_limit_low = -12.0;
   const double boundary_limit_high = 12.0;
   const double zero_epsilon = 1e-9;
   
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
     
     bool on_boundary = (lo < boundary_limit_low) || (lo > boundary_limit_high) ||
       (lr < boundary_limit_low) || (lr > boundary_limit_high);
     
     if (on_boundary) {
       if ((lr < boundary_limit_low) || ((lo < boundary_limit_low) && (lr < 0.0))) {
         current_p0 = getPrbAux_cpp(lo - lr);
         current_p1 = 0.0;
       } else if ((lr > boundary_limit_high) || ((lo < boundary_limit_low) && (lr > 0.0))) {
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
         double discriminant = exp_lo * exp_lo * term1_coeff * term1_coeff +
           4.0 * exp_lr * exp_lo * (1.0 - exp_lo);
         double sqrt_discriminant = std::sqrt(discriminant);
         double numerator_val = -term1_coeff * exp_lo + sqrt_discriminant;
         double denominator_p0 = 2.0 * exp_lr * (1.0 - exp_lo);
         
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
     
     if (clipping) {
       if (R_IsNaN(current_p0)) {
         // p0 remains NA_REAL
       } else if (!R_finite(current_p0)) {
         current_p0 = (current_p0 > 0) ? (1.0 - 1e-15) : 1e-15;
       } else {
         current_p0 = std::fmin(std::fmax(current_p0, 1e-15), 1.0 - 1e-15);
       }
       
       if (R_IsNaN(current_p1)) {
         // p1 remains NA_REAL
       } else if (!R_finite(current_p1)) {
         current_p1 = (current_p1 > 0) ? (1.0 - 1e-15) : 1e-15;
       } else {
         current_p1 = std::fmin(std::fmax(current_p1, 1e-15), 1.0 - 1e-15);
       }
     }
     p0_out[i] = current_p0;
     p1_out[i] = current_p1;
   }
   
   return Rcpp::List::create(Rcpp::Named("p0") = p0_out,
                             Rcpp::Named("p1") = p1_out);
 }


//' @export
// [[Rcpp::export]]
Rcpp::List getProbRR_alt_cpp(const Rcpp::NumericVector& logrr,
                              const Rcpp::NumericVector& logop,
                              bool clipping = true) {
   int n = logrr.size();
   if (logop.size() != n) {
     Rcpp::stop("In prob_rr_alt_worker: logrr and logop must have the same length.");
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
     
     double current_p0 = NA_REAL; // Initialize just in case, though it should be set.
     
     double exp_lr = std::exp(lr);
     double exp_lo = std::exp(lo);
     double exp_lr_plus_2lo = std::exp(lr + 2.0 * lo);
     
     double term_A = 1.0 + exp_lo * (1.0 + exp_lr);
     double term_B_sqrt_arg = term_A * term_A - 4.0 * exp_lr_plus_2lo;
     
     double term_B = std::sqrt(term_B_sqrt_arg);
     double denominator = 2.0 * std::exp(lr + lo);
     
     current_p0 = (term_A - term_B) / denominator;
     
     double current_p1 = exp_lr * current_p0;
     
     if (clipping) {
       if (R_IsNaN(current_p0)) {
         // p0 remains NA_REAL
       } else if (!R_finite(current_p0)) {
         current_p0 = (current_p0 > 0) ? (1.0 - 1e-15) : 1e-15;
       } else {
         current_p0 = std::fmin(std::fmax(current_p0, 1e-15), 1.0 - 1e-15);
       }
       
       if (R_IsNaN(current_p1)) {
         // p1 remains NA_REAL
       } else if (!R_finite(current_p1)) {
         current_p1 = (current_p1 > 0) ? (1.0 - 1e-15) : 1e-15;
       } else {
         current_p1 = std::fmin(std::fmax(current_p1, 1e-15), 1.0 - 1e-15);
       }
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
    // In C++ FISTA, you might handle this error differently or ensure lambda is valid.
    Rcpp::warning("lambda in soft_thres_cpp should be non-negative. Using abs(lambda).");
    lambda = std::fabs(lambda);
  }
  int n = x.size();
  Rcpp::NumericVector result(n);
  
  for (int i = 0; i < n; ++i) {
    double val_x = x[i];
    if (R_IsNaN(val_x)) {
      result[i] = NA_REAL;
    } else if (val_x == 0.0) { // Explicitly handle sign(0)*fmax(0-lambda,0) -> 0
      result[i] = 0.0;
    } else {
      double abs_val_x = std::fabs(val_x);
      double sign_val_x = (val_x > 0) - (val_x < 0); // Efficient sign: 1, -1, or 0
      result[i] = sign_val_x * std::fmax(0.0, abs_val_x - lambda);
    }
  }
  return result;
}

//' @export
// [[Rcpp::export]]
double nllh_cpp(const arma::vec& alpha,
                const arma::vec& beta,
                const arma::mat& va,
                const arma::mat& vb,
                const Rcpp::NumericVector& x_indicator, // Assumed to be 0/1
                const Rcpp::NumericVector& y_outcome,   // Assumed to be 0/1
                int prob_fun) {
  
  int n = y_outcome.size();
  if (n == 0) {
    return 0.0; // Or R_PosInf or an error, consistent with R version if n can be 0
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
      // Handle mismatched dimensions robustly if this function could be called with them
      Rcpp::Rcout << "Warning: Dimension mismatch or empty va for alpha in nllh_cpp. va rows: " << va.n_rows << ", va cols: " << va.n_cols << ", alpha elems: " << alpha.n_elem << ", n: " << n << std::endl;
      // Depending on desired behavior, either stop or ensure logrr_arma remains zeros.
      // For now, it remains zeros if condition not met.
      if(va.n_rows != n && alpha.n_elem > 0) Rcpp::stop("va rows incorrect for n");
      if(va.n_cols != alpha.n_elem && alpha.n_elem > 0) Rcpp::stop("va cols incorrect for alpha");
    }
  }
  
  if (beta.n_elem > 0) {
    if (vb.n_rows == n && vb.n_cols == beta.n_elem) {
      logop_arma = vb * beta;
    } else {
      Rcpp::Rcout << "Warning: Dimension mismatch or empty vb for beta in nllh_cpp. vb rows: " << vb.n_rows << ", vb cols: " << vb.n_cols << ", beta elems: " << beta.n_elem << ", n: " << n << std::endl;
      if(vb.n_rows != n && beta.n_elem > 0) Rcpp::stop("vb rows incorrect for n");
      if(vb.n_cols != beta.n_elem && beta.n_elem > 0) Rcpp::stop("vb cols incorrect for beta");
    }
  }
  
  Rcpp::NumericVector logrr_rcpp = Rcpp::wrap(logrr_arma);
  Rcpp::NumericVector logop_rcpp = Rcpp::wrap(logop_arma);
  
  Rcpp::List ps;
  if (prob_fun == 0) {
    ps = getProbRR_org_cpp(logrr_rcpp, logop_rcpp);
  } else if (prob_fun == 1) {
    ps = getProbRR_alt_cpp(logrr_rcpp, logop_rcpp);
  } else {
    Rcpp::stop("Invalid prob_fun in nllh_cpp");
  }
  
  Rcpp::NumericVector p0_vec = Rcpp::as<Rcpp::NumericVector>(ps["p0"]);
  Rcpp::NumericVector p1_vec = Rcpp::as<Rcpp::NumericVector>(ps["p1"]);
  
  double eps = 1e-15;
  double nll_sum = 0.0;
  
  for (int i = 0; i < n; ++i) {
    double p0 = p0_vec[i];
    double p1 = p1_vec[i];
    
    // If prob_fun_workers returned NA, p0/p1 will be NA_REAL.
    // log(NA_REAL) is NaN. Sum will become NaN.
    // The R_finite check at the end will catch this.
    
    // Local clipping for p0, p1 before log, as in R's nllh
    // This clipping also handles cases where p0/p1 might be exactly 0 or 1
    // from a non-NA calculation in prob_fun_worker if clipping was false there.
    if (!R_IsNaN(p0)) p0 = std::fmax(eps, std::fmin(1.0 - eps, p0));
    if (!R_IsNaN(p1)) p1 = std::fmax(eps, std::fmin(1.0 - eps, p1));
    
    double yi = y_outcome[i];
    
    if (x_indicator[i] == 0) {
      if (R_IsNaN(p0) || p0 <= 0 || (1.0-p0) <=0) { // Ensure log arguments are valid
        nll_sum = NA_REAL; break;
      }
      nll_sum -= (yi * std::log(p0) + (1.0 - yi) * std::log(1.0 - p0));
    } else { // x_indicator[i] == 1
      if (R_IsNaN(p1) || p1 <= 0 || (1.0-p1) <=0) { // Ensure log arguments are valid
        nll_sum = NA_REAL; break;
      }
      nll_sum -= (yi * std::log(p1) + (1.0 - yi) * std::log(1.0 - p1));
    }
    if (R_IsNaN(nll_sum)) break; // Propagate NaN early
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
double penalized_nllh_cpp(const arma::vec& alpha,
                          const arma::vec& beta,
                          const arma::mat& va,
                          const arma::mat& vb,
                          const Rcpp::NumericVector& x_indicator,
                          const Rcpp::NumericVector& y_outcome,
                          double lambda,
                          bool intercept,
                          int prob_fun) {
  
  double unpenalized_nllh = nllh_cpp(alpha, beta, va, vb, x_indicator, y_outcome,
                                     prob_fun);
  
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
      // Assuming beta might also have an intercept conceptually, or always penalize all of beta
      // R code was sum(abs(beta[if(intercept) -1 else TRUE]))
      // This implies if intercept=TRUE, beta[1] is also an intercept. Adjust if this assumption is wrong.
      // If beta has no intercept, start_idx_beta should always be 0.
      int start_idx_beta = (intercept && beta.n_elem > 0) ? 1 : 0; 
      for (arma::uword i = start_idx_beta; i < beta.n_elem; ++i) {
        l1_norm_beta += std::fabs(beta[i]);
      }
    }
    penalty = lambda * (l1_norm_alpha + l1_norm_beta);
  }
  
  return unpenalized_nllh + penalty;
}






//' @export
 // [[Rcpp::export]]
 arma::vec grad_nll_alpha_cpp(
     const arma::vec& alpha_eval,
     const arma::vec& beta_current,
     const arma::mat& va,
     const arma::mat& vb,
     const Rcpp::NumericVector& x_indicator,
     const Rcpp::NumericVector& y_outcome,
     int prob_fun_selector, // Changed name to avoid conflict with R's 'prob_fun' argument if it were a function
     bool clipping_for_prob_fun_passed_to_nllh = true) { // Added for completeness, ensure nllh_cpp uses it
   
   if (alpha_eval.n_elem == 0) {
     return arma::vec();
   }
   Rcpp::NumericVector alpha_eval_rcpp = Rcpp::wrap(alpha_eval);
   
   auto f_alpha = [&](Rcpp::NumericVector current_alpha_rcpp_lambda_arg) -> double {
     arma::vec current_alpha_arma = Rcpp::as<arma::vec>(current_alpha_rcpp_lambda_arg);
     // Ensure nllh_cpp is called with the correct number of arguments
     // If your nllh_cpp takes 7 args (no clipping_for_prob_fun), remove it from this call:
     return nllh_cpp(current_alpha_arma, beta_current, va, vb, 
                     x_indicator, y_outcome, prob_fun_selector
                       // If your nllh_cpp has the 8th arg for clipping:
                       // , clipping_for_prob_fun_passed_to_nllh
     );
   };
   
   Rcpp::List fntl_output_list = Rcpp::wrap(fntl::gradient(f_alpha, alpha_eval_rcpp));
   Rcpp::NumericVector grad_rcpp;
   
   if (fntl_output_list.size() > 0 && Rcpp::is<Rcpp::NumericVector>(fntl_output_list[0])) {
     grad_rcpp = Rcpp::as<Rcpp::NumericVector>(fntl_output_list[0]);
   } else {
     // This Rcout and Rf_PrintValue can help you debug the list structure if needed
     Rcpp::Rcout << "Unexpected structure for fntl::gradient output list in grad_nll_alpha_cpp:" << std::endl;
     Rf_PrintValue(fntl_output_list); // Prints the R list structure to the R console
     Rcpp::stop("Could not extract gradient vector from fntl output (expected NumericVector as first element).");
   }
   
   return Rcpp::as<arma::vec>(grad_rcpp);
 }

//' @export
// [[Rcpp::export]]
arma::vec grad_nll_beta_cpp(
     const arma::vec& alpha_current,
     const arma::vec& beta_eval,
     const arma::mat& va,
     const arma::mat& vb,
     const Rcpp::NumericVector& x_indicator,
     const Rcpp::NumericVector& y_outcome,
     int prob_fun_selector,
     bool clipping_for_prob_fun_passed_to_nllh = true) { // Added for completeness
   
   if (beta_eval.n_elem == 0) {
     return arma::vec();
   }
   Rcpp::NumericVector beta_eval_rcpp = Rcpp::wrap(beta_eval);
   
   auto f_beta = [&](Rcpp::NumericVector current_beta_rcpp_lambda_arg) -> double {
     arma::vec current_beta_arma = Rcpp::as<arma::vec>(current_beta_rcpp_lambda_arg);
     // Ensure nllh_cpp is called with the correct number of arguments
     return nllh_cpp(alpha_current, current_beta_arma, va, vb, 
                     x_indicator, y_outcome, prob_fun_selector
                       // If your nllh_cpp has the 8th arg for clipping:
                       // , clipping_for_prob_fun_passed_to_nllh
     );
   };
   
   Rcpp::List fntl_output_list = Rcpp::wrap(fntl::gradient(f_beta, beta_eval_rcpp));
   Rcpp::NumericVector grad_rcpp;
   
   if (fntl_output_list.size() > 0 && Rcpp::is<Rcpp::NumericVector>(fntl_output_list[0])) {
     grad_rcpp = Rcpp::as<Rcpp::NumericVector>(fntl_output_list[0]);
   } else {
     Rcpp::Rcout << "Unexpected structure for fntl::gradient output list in grad_nll_beta_cpp:" << std::endl;
     Rf_PrintValue(fntl_output_list);
     Rcpp::stop("Could not extract gradient vector from fntl output (expected NumericVector as first element).");
   }
   
   return Rcpp::as<arma::vec>(grad_rcpp);
 }






struct StepFistaResult {
  arma::vec value_new; // x_{k+1}
  double t_new;        // t_{k+1}
  // arma::vec y_eval_point; // y_{k+1} where gradient was evaluated (optional return)
  // double step_size_used; // If using line search, return the step size found
};

// --- step_fista_cpp ---
// Performs one FISTA step for either alpha or beta
StepFistaResult step_fista_cpp(
    const arma::vec& current_param_val,      // Current iterate (x_k)
    const arma::vec& other_param_val,        // The other parameter (beta or alpha)
    const arma::vec& prev_param_val,         // Previous iterate (x_{k-1})
    const std::string& opt_target,           // "alpha" or "beta"
    double step_size,                        // step size 's'
    double lambda,                           // penalty strength
    double t_current,                        // t_k
    bool intercept,                          // Penalize intercept?
    const arma::mat& va,
    const arma::mat& vb,
    const Rcpp::NumericVector& x_indicator,
    const Rcpp::NumericVector& y_outcome,
    int prob_fun_selector) {
  
  StepFistaResult result;
  
  // Calculate next t value: t_{k+1} = (1 + sqrt(1 + 4*t_k^2)) / 2
  result.t_new = (1.0 + std::sqrt(1.0 + 4.0 * t_current * t_current)) / 2.0;
  double momentum_coeff = (t_current - 1.0) / result.t_new;
  
  // Calculate momentum point: y_{k+1} = x_k + momentum_coeff * (x_k - x_{k-1})
  arma::vec y_eval_point;
  if (current_param_val.n_elem > 0 && prev_param_val.n_elem == current_param_val.n_elem) {
    y_eval_point = current_param_val + momentum_coeff * (current_param_val - prev_param_val);
  } else {
    y_eval_point = current_param_val; // No momentum (e.g., first iteration or empty params)
  }
  // result.y_eval_point = y_eval_point; // Store if needed outside
  
  // --- Gradient Calculation ---
  // Gradient of unpenalized nllh evaluated at the momentum point y_{k+1}
  arma::vec gradient_val;
  if (opt_target == "alpha") {
    gradient_val = grad_nll_alpha_cpp(y_eval_point, other_param_val,
                                      va, vb, x_indicator, y_outcome,
                                      prob_fun_selector);
  } else if (opt_target == "beta") {
    gradient_val = grad_nll_beta_cpp(other_param_val, y_eval_point,
                                     va, vb, x_indicator, y_outcome,
                                     prob_fun_selector);
  } else {
    Rcpp::stop("Invalid opt_target in step_fista_cpp. Must be 'alpha' or 'beta'.");
  }
  
  // Clean NA/NaN/Inf gradients (R code replaced NA with 0)
  gradient_val.elem(arma::find_nonfinite(gradient_val)).zeros();
  
  // --- Proximal Update (Gradient Descent + Soft Thresholding) ---
  // Input to soft-thresholding: y_{k+1} - step_size * gradient(y_{k+1})
  arma::vec prox_input_arma = y_eval_point - step_size * gradient_val;
  
  // soft_thres_cpp expects Rcpp::NumericVector, call it
  Rcpp::NumericVector prox_input_rcpp = Rcpp::wrap(prox_input_arma);
  Rcpp::NumericVector value_new_rcpp = soft_thres_cpp(prox_input_rcpp, lambda * step_size);
  
  // Convert result back to Armadillo vec
  result.value_new = Rcpp::as<arma::vec>(value_new_rcpp);
  
  // Apply intercept logic (don't penalize first element)
  if (intercept && result.value_new.n_elem > 0) {
    // R code: value_new[1] <- input[1]
    // C++: result.value_new[0] should not be thresholded -> take value from prox_input_arma
    result.value_new[0] = prox_input_arma[0]; // 0-based index
  }
  
  // result.step_size_used = step_size; // Add if using line search later
  return result;
}

// --- stop_crit_cpp ---
// Basic stopping criteria based on parameter change and optionally gradient norm
bool stop_crit_cpp(const arma::vec& grad_alpha, // Gradient at current alpha, beta
                   const arma::vec& grad_beta,
                   const arma::vec& alpha,     // Current alpha (x_k)
                   const arma::vec& beta,      // Current beta (x_k)
                   const arma::vec& last_alpha,// Previous alpha (x_{k-1})
                   const arma::vec& last_beta, // Previous beta (x_{k-1})
                   double tol_param_change,    // Tolerance for parameter change (e.g., max abs diff)
                   double tol_grad_norm,       // Tolerance for gradient norm (e.g., max abs component)
                   bool check_grad) {          // Whether to check gradient norm
  
  // Check parameter convergence (relative or absolute change can be used)
  // Using max absolute difference here for simplicity
  double alpha_diff = (alpha.n_elem > 0 && last_alpha.n_elem == alpha.n_elem) ?
  arma::norm(alpha - last_alpha, "fro") : 0.0;
  double beta_diff = (beta.n_elem > 0 && last_beta.n_elem == beta.n_elem) ?
  arma::norm(beta - last_beta, "fro") : 0.0;
  
  // R code's stop_crit function isn't shown, but often requires *both* params to converge
  bool params_converged = (alpha_diff < tol_param_change) && (beta_diff < tol_param_change);
  
  if (params_converged) {
    return true;
  }
  
  // Optional: Check gradient norm convergence (using gradient at current x_k)
  if (check_grad) {
    double grad_alpha_norm = (grad_alpha.n_elem > 0 && arma::is_finite(grad_alpha)) ?
    arma::norm(grad_alpha, "fro") : std::numeric_limits<double>::infinity();
    double grad_beta_norm = (grad_beta.n_elem > 0 && arma::is_finite(grad_beta)) ?
    arma::norm(grad_beta, "fro") : std::numeric_limits<double>::infinity();
    
    if (grad_alpha_norm < tol_grad_norm && grad_beta_norm < tol_grad_norm) {
      return true;
    }
  }
  
  return false; // Continue iterations
}


// --- fista_opt2_cpp ---
//' @export
// [[Rcpp::export]]
Rcpp::List fista_opt2_cpp(
    Rcpp::NumericVector alpha_start_rcpp,
    Rcpp::NumericVector beta_start_rcpp,
    double step_size_alpha,           // Initial step sizes
    double step_size_beta,
    double lambda,                    // Penalty parameter
    bool intercept,                   // Penalize intercept?
    int max_iter,                     // Max iterations
    Rcpp::NumericMatrix va_rcpp,
    Rcpp::NumericMatrix vb_rcpp,
    Rcpp::NumericVector x_indicator,
    Rcpp::NumericVector y_outcome,
    int prob_fun_selector,            // 0 for org, 1 for alt
    bool clipping_for_prob_fun = true,// Clipping for internal prob funcs
    bool eval_grad_for_output_and_stop_crit = true, // Whether to calc grad for output/stopping
    double tol_param_change = 1e-6,   // Stopping tolerance for parameters
    double tol_grad_norm = 1e-6) {    // Stopping tolerance for gradient norm
  
  // --- Input Conversions and Initializations ---
  arma::vec alpha = Rcpp::as<arma::vec>(alpha_start_rcpp);
  arma::vec beta = Rcpp::as<arma::vec>(beta_start_rcpp);
  arma::mat va = Rcpp::as<arma::mat>(va_rcpp);
  arma::mat vb = Rcpp::as<arma::mat>(vb_rcpp);
  
  arma::vec alpha_prev = alpha; // x_{k-1} for alpha
  arma::vec beta_prev = beta;   // x_{k-1} for beta
  
  double t_alpha = 1.0;
  double t_beta = 1.0;
  
  // History storage (using std::vector for flexibility)
  std::vector<arma::vec> alphas_hist;
  std::vector<arma::vec> betas_hist;
  std::vector<arma::vec> g_alphas_hist; // Gradient at end of iteration
  std::vector<arma::vec> g_betas_hist;  // Gradient at end of iteration
  std::vector<double> nllh_results_hist; // Penalized NLLH
  
  int actual_iter = 0;
  bool converged = false;
  
  // --- Main FISTA Loop ---
  for (int iter = 0; iter < max_iter; ++iter) {
    actual_iter = iter + 1; // R-like iteration count (starts at 1)
    Rcpp::checkUserInterrupt(); // Allow user to interrupt long computation
    
    arma::vec alpha_at_iter_start = alpha; // Store x_k before updates
    arma::vec beta_at_iter_start = beta;
    
    // --- Alpha Update ---
    StepFistaResult res_alpha = step_fista_cpp(
      alpha, beta, alpha_prev, "alpha", // Pass current alpha(x_k), current beta, prev alpha(x_k-1)
      step_size_alpha, lambda, t_alpha, intercept,
      va, vb, x_indicator, y_outcome,
      prob_fun_selector
    );
    alpha_prev = alpha_at_iter_start; // Update previous alpha for next iter's momentum
    alpha = res_alpha.value_new;      // Update current alpha (now x_{k+1} conceptually)
    t_alpha = res_alpha.t_new;
    
    // --- Beta Update ---
    // Note: Uses the *updated* alpha from this iteration
    StepFistaResult res_beta = step_fista_cpp(
      beta, alpha, beta_prev, "beta",  // Pass current beta(x_k), updated alpha, prev beta(x_k-1)
      step_size_beta, lambda, t_beta, intercept,
      va, vb, x_indicator, y_outcome,
      prob_fun_selector
    );
    beta_prev = beta_at_iter_start; // Update previous beta for next iter's momentum
    beta = res_beta.value_new;      // Update current beta (now x_{k+1} conceptually)
    t_beta = res_beta.t_new;
    
    // --- Store History and Check Convergence ---
    arma::vec grad_alpha_final_iter; // Gradient at the *updated* alpha, beta
    arma::vec grad_beta_final_iter;
    bool check_grad_stop = false; // Only check grad if evaluated
    
    if (eval_grad_for_output_and_stop_crit) {
      check_grad_stop = true;
      // Calculate gradient at the newly updated point (alpha, beta)
      grad_alpha_final_iter = grad_nll_alpha_cpp(
        alpha, beta, va, vb, x_indicator, y_outcome,
        prob_fun_selector);
      grad_alpha_final_iter.elem(arma::find_nonfinite(grad_alpha_final_iter)).zeros(); // Clean
      
      grad_beta_final_iter = grad_nll_beta_cpp(
        alpha, beta, va, vb, x_indicator, y_outcome,
        prob_fun_selector);
      grad_beta_final_iter.elem(arma::find_nonfinite(grad_beta_final_iter)).zeros(); // Clean
      
      g_alphas_hist.push_back(grad_alpha_final_iter);
      g_betas_hist.push_back(grad_beta_final_iter);
    } else {
      // Store empty vectors or vectors of NA if grads not evaluated
      g_alphas_hist.push_back(arma::vec());
      g_betas_hist.push_back(arma::vec());
    }
    
    // Calculate penalized NLLH at end of iteration
    double nllh_iter = penalized_nllh_cpp(
      alpha, beta, va, vb, x_indicator, y_outcome,
      lambda, intercept, prob_fun_selector
    );
    
    // Store other history items
    alphas_hist.push_back(alpha);
    betas_hist.push_back(beta);
    nllh_results_hist.push_back(nllh_iter);
    
    // Check stopping criteria using the parameters from *start* of this iter (x_k)
    // and the end of the previous iter (x_{k-1}) which were saved as alpha/beta_at_iter_start
    // R code compared end-of-iter alpha/beta with previous end-of-iter (last_alpha/beta)
    // This corresponds to comparing alpha/beta with alpha_at_iter_start/beta_at_iter_start
    if (iter > 0) { // Start checking after first iteration
      if (stop_crit_cpp(grad_alpha_final_iter, grad_beta_final_iter, // grad at end of current iter
                        alpha, beta,                  // alpha/beta at end of current iter
                        alpha_at_iter_start, beta_at_iter_start, // alpha/beta at start of current iter
                        tol_param_change, tol_grad_norm,
                        check_grad_stop)) {
        converged = true;
        break;
      }
    }
  } // End main loop
  
  if (!converged && actual_iter == max_iter) {
    Rcpp::warning("FISTA did not converge after %d iterations.", max_iter);
  }
  
  // --- Prepare Output ---
  // Convert history vectors of vectors/doubles to Rcpp types for output
  // For simplicity, just returning final values and basic info.
  // Returning full history requires converting vector<arma::vec> to matrix.
  
  return Rcpp::List::create(
    Rcpp::Named("alpha") = Rcpp::wrap(alpha),
    Rcpp::Named("beta") = Rcpp::wrap(beta),
    Rcpp::Named("iterations") = actual_iter,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("nllh_final") = nllh_results_hist.empty() ? R_NaN : nllh_results_hist.back()
    // Add history if implemented:
    // Rcpp::Named("nllh_history") = Rcpp::wrap(nllh_results_hist)
    // Rcpp::Named("alphas_history") = wrap_history_to_matrix(alphas_hist),
    // Rcpp::Named("betas_history") = wrap_history_to_matrix(betas_hist),
    // Rcpp::Named("grad_alpha_history") = wrap_history_to_matrix(g_alphas_hist),
    // Rcpp::Named("grad_beta_history") = wrap_history_to_matrix(g_betas_hist)
  );
}

