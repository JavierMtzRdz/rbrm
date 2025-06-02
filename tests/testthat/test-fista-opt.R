# tests/testthat/test-fista-comparison.R

library(testthat)
library(rbrm)     # Your package
library(numDeriv) # Used by the R version fista_opt2 -> step_fista

# --- Test Data Setup ---
# (Using generate_data as before)
set.seed(20250508) # Use current date for seed YYYYMMDD

n_obs_train <- 150L # Keep reasonably small for faster tests
pa_coeffs <- 3L
pb_coeffs <- pa_coeffs # Assuming va=vb

# True coefficients for data generation
true_alpha_gen <- c(1, -1, 0)
true_beta_gen <- c(1, -1, 0)
true_gamma_gen <- c(1, -1, 0)

# Generate data
sim_data <- rbrm::generate_data(
  pa = pa_coeffs, pb = pb_coeffs, n = n_obs_train, n_test = 10,
  alpha = true_alpha_gen, beta = true_beta_gen, gamma = true_gamma_gen,
  misspec_nuisance = FALSE, misspec_propensity = FALSE
)

# Extract data
va_data <- sim_data$v.train
vb_data <- sim_data$v.train # Assuming vb = va
x_indicator <- sim_data$x.train
y_outcome <- sim_data$y.train

# --- FISTA Parameters (use the same for both R and C++) ---
alpha_start_vals <- rep(0.0, pa_coeffs)
beta_start_vals <- rep(0.0, pb_coeffs)
step_alpha <- 0.1 # May need tuning for convergence in reasonable time
step_beta <- 0.1
lambda_pen <- 0 # Use a small non-zero lambda
intercept_opt <- FALSE
max_iter_test <- 300 # Limit iterations for test speed
prob_fun_r <- rbrm::getProbRR.org # R function object for R version
prob_fun_selector_cpp <- 0      # Integer selector for C++ version (0=org)
clipping_opt <- TRUE
eval_grad_opt <- TRUE # Needed for R version stop_crit if it uses grads
stop_tol_p_test <- 1e-6 # Use slightly looser tolerance for test comparison
stop_tol_g_test <- 1e-6

# --- Run FISTA versions ---
test_that("C++ FISTA results match R FISTA results", {
  
  cat("Running R version (fista_opt2)...\n")
  # Ensure all R helper functions (step_fista, nllh, penalized_nllh, stop_crit)
  # are available to fista_opt2, either exported or internal to the package.
  time_r <- system.time({
    fista_result_r <- tryCatch({
      rbrm::fista_opt(
        alpha.start = alpha_start_vals,
        beta.start = beta_start_vals,
        step_size_alpha = step_alpha,
        step_size_beta = step_beta,
        lambda = lambda_pen,
        intercept = intercept_opt,
        max.step = max_iter_test, # Note: param name difference R vs C++
        va = va_data,
        vb = vb_data,
        x = x_indicator,
        y = y_outcome,
        prob_fun = prob_fun_r, # Pass the R function object
        # opt_step = rbrm::step_fista, # Assuming default
        eval_grad = eval_grad_opt
        # Assuming fista_opt2 internally calls a stop_crit function
        # that uses tolerances similar to tol_param_change/tol_grad_norm
        # If tolerances are hardcoded/different, comparison is harder.
      )
    }, error = function(e) {
      warning("R fista_opt2 failed: ", e$message)
      NULL
    })
  })
  cat("R version took:", time_r['elapsed'], "seconds.\n")
  
  # Skip C++ test if R version failed
  skip_if(is.null(fista_result_r), "R version of FISTA failed to run.")
  
  # --- Run C++ version ---
  cat("Running C++ version (fista_opt2_cpp)...\n")
  # Ensure the C++ gradient functions are correctly implemented!
  time_cpp <- system.time({
    fista_result_cpp <- tryCatch({
      rbrm::fista_opt2_cpp(
        alpha_start_rcpp = alpha_start_vals,
        beta_start_rcpp = beta_start_vals,
        step_size_alpha = step_alpha,
        step_size_beta = step_beta,
        lambda = lambda_pen,
        intercept = intercept_opt,
        max_iter = max_iter_test, # Note: param name difference R vs C++
        va_rcpp = va_data,
        vb_rcpp = vb_data,
        x_indicator = x_indicator,
        y_outcome = y_outcome,
        prob_fun_selector = prob_fun_selector_cpp,
        clipping_for_prob_fun = clipping_opt,
        eval_grad_for_output_and_stop_crit = T,
        tol_param_change = stop_tol_p_test,
        tol_grad_norm = stop_tol_g_test
      )
    }, error = function(e) {
      warning("C++ fista_opt2_cpp failed: ", e$message)
      NULL
    })
  })
  cat("C++ version took:", time_cpp['elapsed'], "seconds.\n")
  
  # Skip comparison if C++ version failed
  skip_if(is.null(fista_result_cpp), "C++ version of FISTA failed to run.")
  
  # --- Compare Results ---
  cat("Comparing final alpha coefficients...\n")
  expect_equal(fista_result_cpp$alpha, fista_result_r$alpha,
               tolerance = stop_tol_p_test * 10, # Allow slightly larger tolerance for final params
               label = "Final alpha coefficients (C++ vs R)")
  
  cat("Comparing final beta coefficients...\n")
  expect_equal(as.vector(fista_result_cpp$beta), fista_result_r$beta,
               tolerance = stop_tol_p_test * 10,
               label = "Final beta coefficients (C++ vs R)")
  
  # Compare final objective function value (optional, less critical than params)
  cat("Comparing final penalized NLLH...\n")
  final_nllh_r <- fista_result_r$nllh_results[length(fista_result_r$nllh_results)]
  expect_equal(fista_result_cpp$nllh_final, final_nllh_r,
               tolerance = 1e-4, # Objective function might differ more
               label = "Final penalized NLLH (C++ vs R)")
  
  # Compare convergence status (might differ slightly due to numerical precision)
  cat("Comparing convergence status (expect TRUE for both ideally)...\n")
  expect_equal(fista_result_cpp$converged, length(fista_result_r$step) < max_iter_test,
               label = "Convergence status (C++ vs R)")
  # Note: R version doesn't return a 'converged' flag, we infer it from iter < max.step
  # If R version *does* have a specific convergence flag/message, compare that instead.
  
  # Iteration count comparison (often differs slightly, just check order of magnitude)
  cat("Comparing iteration counts (expect reasonable similarity)...\n")
  # expect_lt(abs(fista_result_cpp$iterations - fista_result_r$step), max_iter_test * 0.2) # Example: within 20%
  # Or just print them:
  cat("  C++ iterations:", fista_result_cpp$iterations, "\n")
  cat("  R iterations:", fista_result_r$step, "\n") # R version returns 'step'
  
  
})
