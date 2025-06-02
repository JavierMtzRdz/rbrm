# tests/testthat/test-gradients.R

library(testthat)
library(rbrm)     # Your package, exports C++ functions and R functions like generate_data, nllh, etc.
library(numDeriv) # For reference gradient calculation


# --- Test Data Setup using generate_data ---
set.seed(12345) # for reproducibility

n_obs_train <- 50L # Number of training observations
n_obs_test <- 10L # Not used in these nllh gradient tests, but generate_data produces it
pa_coeffs <- 3L   # Number of covariates for v, and thus length of alpha and beta
pb_coeffs <- pa_coeffs # Assuming vb uses the same covariates as va for simplicity here,
# so beta will also have pa_coeffs length. Adjust if vb is different.

# True coefficients for data generation
true_alpha_gen <- rnorm(pa_coeffs, sd = 0.5)
true_beta_gen <- rnorm(pa_coeffs, sd = 0.3) # Must match number of columns in vb (here, pa_coeffs)
true_gamma_gen <- rnorm(pa_coeffs, sd = 0.4) # For propensity score model, length pa_coeffs

sim_data_list <- rbrm::generate_data(
  pa = pa_coeffs,
  pb = pb_coeffs, # Used by generate_data for v.train dimensions if different, but here v.train is n x pa
  n = n_obs_train,
  n_test = n_obs_test,
  alpha = true_alpha_gen,
  beta = true_beta_gen,
  gamma = true_gamma_gen,
  misspec_nuisance = FALSE,   # Use correctly specified models for simpler gradient checking
  misspec_propensity = FALSE
)

# Data for nllh function (using training set from sim_data_list)
# In your generate_data, v.train is used for both alpha and beta terms.
# So va_mat and vb_mat will both be sim_data_list$v.train
va_test_mat <- sim_data_list$v.train
vb_test_mat <- sim_data_list$v.train # Assuming vb is the same as va
x_indicator_test_vec <- sim_data_list$x.train
y_outcome_test_vec <- sim_data_list$y.train

# Parameter values at which to evaluate the gradient
# These can be different from true_alpha_gen/true_beta_gen
alpha_eval_params <- rnorm(pa_coeffs, mean = 0.1, sd = 0.2)
beta_eval_params <- rnorm(pb_coeffs, mean = -0.1, sd = 0.2) # Should match length of beta_true_gen

clipping_prob_test <- TRUE # Clipping setting for the probability functions

# --- Test for Alpha Gradient (using getProbRR.org logic) ---
test_that("Alpha gradient (org prob_fun) from C++ matches numDeriv from R", {
  prob_fun_selector_test <- 0 # 0 for getProbRR.org
  
  # R wrapper for nllh to use with numDeriv::grad (for alpha)
  nllh_r_for_alpha_grad <- function(alpha_val) {
    rbrm::nllh( # Call the R version of nllh from your package
      alpha = alpha_val,
      beta = beta_eval_params, # Fixed beta at the evaluation point
      va = va_test_mat,
      vb = vb_test_mat,
      x = x_indicator_test_vec,
      y = y_outcome_test_vec
    )
  }
  
  # Calculate gradient using numDeriv
  grad_alpha_r_numderiv <- numDeriv::grad(
    func = nllh_r_for_alpha_grad,
    x = alpha_eval_params, # Point at which to evaluate the gradient
    method = "Richardson"
  )
  
  # Calculate gradient using your C++ function
  grad_alpha_cpp <- rbrm::grad_nll_alpha_cpp(
    alpha = alpha_eval_params,
    beta = beta_eval_params, # Fixed beta
    va = va_test_mat,
    vb = vb_test_mat,
    x_indicator = x_indicator_test_vec,
    y_outcome = y_outcome_test_vec,
    prob_fun = prob_fun_selector_test
  )
  
  # bench::mark(num = numDeriv::grad(
  #   func = nllh_r_for_alpha_grad,
  #   x = alpha_eval_params, # Point at which to evaluate the gradient
  #   method = "simple"
  # ),
  # rcpp = rbrm::grad_nll_alpha_cpp(
  #   alpha = alpha_eval_params,
  #   beta = beta_eval_params, # Fixed beta
  #   va = va_test_mat,
  #   vb = vb_test_mat,
  #   x_indicator = x_indicator_test_vec,
  #   y_outcome = y_outcome_test_vec,
  #   prob_fun = prob_fun_selector_test
  # ),
  # check = F)
  
  expect_equal(as.vector(grad_alpha_cpp), grad_alpha_r_numderiv, tolerance = 1e-3,
               label = "C++ alpha gradient (fntl, org)",
               expected.label = "R numDeriv alpha gradient (org)")
})

# --- Test for Beta Gradient (using getProbRR.org logic) ---
test_that("Beta gradient (org prob_fun) from C++ matches numDeriv from R", {
  prob_fun_selector_test <- 0 # 0 for getProbRR.org
  
  # R wrapper for nllh to use with numDeriv::grad (for beta)
  nllh_r_for_beta_grad <- function(beta_val) {
    rbrm::nllh( # Call the R version of nllh
      alpha = alpha_eval_params, # Fixed alpha
      beta = beta_val,
      va = va_test_mat,
      vb = vb_test_mat,
      x = x_indicator_test_vec,
      y = y_outcome_test_vec,
      prob_fun = rbrm::getProbRR.org
    )
  }
  
  # Calculate gradient using numDeriv
  grad_beta_r_numderiv <- numDeriv::grad(
    func = nllh_r_for_beta_grad,
    x = beta_eval_params, # Point at which to evaluate the gradient
    method = "Richardson"
  )
  
  # Calculate gradient using your C++ function
  grad_beta_cpp <- rbrm::grad_nll_beta_cpp(
    alpha = alpha_eval_params, # Fixed alpha
    beta = beta_eval_params,
    va = va_test_mat,
    vb = vb_test_mat,
    x_indicator = x_indicator_test_vec,
    y_outcome = y_outcome_test_vec,
    prob_fun = prob_fun_selector_test,
  )
  
  expect_equal(as.vector(grad_beta_cpp), grad_beta_r_numderiv, tolerance = 1e-2,
               label = "C++ beta gradient (fntl, org)",
               expected.label = "R numDeriv beta gradient (org)")
})


# You can add similar tests for prob_fun_selector = 1 (getProbRR.alt)
# by changing the prob_fun_selector_test and the prob_fun passed to nllh.R.
# For example:
test_that("Alpha gradient (alt prob_fun) from C++ matches numDeriv from R", {
  prob_fun_selector_test <- 1 # 1 for getProbRR.alt
  
  nllh_r_for_alpha_grad_alt <- function(alpha_val) {
    rbrm::nllh(
      alpha = alpha_val, beta = beta_eval_params,
      va = va_test_mat, vb = vb_test_mat,
      x = x_indicator_test_vec, y = y_outcome_test_vec,
      prob_fun = rbrm::getProbRR.alt # Use R version of getProbRR.alt
    )
  }
  
  grad_alpha_r_numderiv_alt <- numDeriv::grad(func = nllh_r_for_alpha_grad_alt, x = alpha_eval_params)
  
  grad_alpha_cpp_alt <- rbrm::grad_nll_alpha_cpp(
    alpha = alpha_eval_params, beta = beta_eval_params,
    va = va_test_mat, vb = vb_test_mat,
    x_indicator = x_indicator_test_vec, y_outcome = y_outcome_test_vec,
    prob_fun = prob_fun_selector_test
  )
  
  expect_equal(as.vector(grad_alpha_cpp_alt), grad_alpha_r_numderiv_alt, tolerance = 1e-3)
})

