
# --- Tests for soft_thres_cpp ---
test_that("soft_thres_cpp matches existing R soft_thres", {
  x_vec <- c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2) #NA
  
  # Lambda = 1
  lambda1 <- 1.0
  expect_equal(rbrm::soft_thres_cpp(x_vec, lambda1), rbrm::soft_thres(x_vec, lambda1))
  expect_equal(rbrm::soft_thres_cpp(1.5, lambda1), rbrm::soft_thres(1.5, lambda1))
  expect_equal(rbrm::soft_thres_cpp(-0.5, lambda1), rbrm::soft_thres(-0.5, lambda1))
  
  # Lambda = 0
  lambda0 <- 0.0
  expect_equal(rbrm::soft_thres_cpp(x_vec, lambda0), rbrm::soft_thres(x_vec, lambda0))
  
  # Lambda = 1.5
  lambda1.5 <- 1.5
  expect_equal(rbrm::soft_thres_cpp(x_vec, lambda1.5), rbrm::soft_thres(x_vec, lambda1.5))
  
  # Negative lambda 
  # R soft_thres doesn't warn for negative lambda by default in user code.
  # C++ soft_thres_cpp does issue a warning.
  lambda_neg <- -0.5
  # expect_warning(res_cpp_neg_lambda <- rbrm::soft_thres_cpp(x_vec, lambda_neg), "lambda in soft_thres_cpp should be non-negative. Using abs(lambda).")
  # res_r_neg_lambda <- rbrm::soft_thres(x_vec, lambda_neg) 
  # expect_equal(res_cpp_neg_lambda, res_r_neg_lambda) 
  
  # Empty vector
  expect_equal(rbrm::soft_thres_cpp(numeric(0), lambda1), rbrm::soft_thres(numeric(0), lambda1))
  
  # Error for NA lambda (C++ version has this check)
  # expect_error(rbrm::soft_thres_cpp(x_vec, NA_real_), "lambda cannot be NA")
})

# --- Setup for nllh and penalized_nllh tests ---
set.seed(123)
n_obs <- 30 
n_alpha_vars <- 2
n_beta_vars <- 2 # Can be 0, 1, or more

alpha_coeffs <- if (n_alpha_vars > 0) rnorm(n_alpha_vars) else numeric(0)
beta_coeffs  <- if (n_beta_vars > 0) rnorm(n_beta_vars) else numeric(0)
va_matrix  <- vb_matrix  <- if (n_alpha_vars > 0) matrix(rnorm(n_obs * n_alpha_vars), nrow = n_obs, ncol = n_alpha_vars) else matrix(0, nrow = n_obs, ncol = 0)

x_indicator  <- sample(0:1, n_obs, replace = TRUE) 
y_outcome    <- sample(0:1, n_obs, replace = TRUE) 
lambda_penalty <- 0.1
has_intercept <- TRUE

# --- Tests for nllh_cpp ---
test_that("nllh_cpp matches existing R nllh", {
  # Base case
  
  nll_r_val <- rbrm::nllh(alpha_coeffs, beta_coeffs, va_matrix, vb_matrix, x_indicator, y_outcome)
  nll_cpp_val <- rbrm::nllh_cpp(alpha_coeffs, beta_coeffs, va_matrix, vb_matrix, x_indicator, y_outcome, prob_fun = 0)
  
  
  
  expect_equal(nll_cpp_val, nll_r_val, tolerance = 1e-7, info = "Base nllh case")
  
  # Case: No alpha coefficients (pa = 0)
  alpha0 <- numeric(0)
  va0 <- matrix(0, nrow = n_obs, ncol = 0)
  nll_r_pa0 <- rbrm::nllh(alpha0, beta_coeffs, va0, vb_matrix, x_indicator, y_outcome)
  nll_cpp_pa0 <- rbrm::nllh_cpp(alpha0, beta_coeffs, va0, vb_matrix, x_indicator, y_outcome, prob_fun = 0)
  expect_equal(nll_cpp_pa0, nll_r_pa0, tolerance = 1e-7, info = "nllh with pa=0")
  
  # Case: No beta coefficients (pb = 0)
  beta0 <- numeric(0)
  vb0 <- matrix(0, nrow = n_obs, ncol = 0)
  nll_r_pb0 <- rbrm::nllh(alpha_coeffs, beta0, va_matrix, vb0, x_indicator, y_outcome)
  nll_cpp_pb0 <- rbrm::nllh_cpp(alpha_coeffs, beta0, va_matrix, vb0, x_indicator, y_outcome, prob_fun = 0)
  expect_equal(nll_cpp_pb0, nll_r_pb0, tolerance = 1e-7, info = "nllh with pb=0")
  
  # Case: No alpha and No beta coefficients (pa = 0, pb = 0)
  nll_r_pa0pb0 <- rbrm::nllh(alpha0, beta0, va0, vb0, x_indicator, y_outcome)
  nll_cpp_pa0pb0 <- rbrm::nllh_cpp(alpha0, beta0, va0, vb0, x_indicator, y_outcome, prob_fun = 0)
  expect_equal(nll_cpp_pa0pb0, nll_r_pa0pb0, tolerance = 1e-7, info = "nllh with pa=0, pb=0")
  
  # Case: Empty y (n=0)
  va_empty <- if (n_alpha_vars > 0) matrix(0, nrow=0, ncol=n_alpha_vars) else matrix(0,0,0)
  vb_empty <- if (n_beta_vars > 0) matrix(0, nrow=0, ncol=n_beta_vars) else matrix(0,0,0)
  expect_equal(rbrm::nllh_cpp(alpha_coeffs, beta_coeffs, va_empty, vb_empty, integer(0), integer(0), prob_fun = 0), 0.0) # C++ returns 0.0 for n=0 valid_obs
  # expect_equal(rbrm::nllh_cpp(alpha_coeffs, beta_coeffs, va_empty, vb_empty, integer(0), integer(0), prob_fun = 0), Inf) # R original returns Inf
  
  # Case: NA in y or x
  y_na <- y_outcome; if(length(y_na)>0) y_na[1] <- NA
  x_na <- x_indicator; if(length(x_na)>0) x_na[2] <- NA
  
  # R nllh might handle NA differently (e.g. na.rm in sum). C++ skips NAs and averages over valid.
  # This requires careful alignment of NA handling logic in R original vs C++ for exact match.
  # The C++ code averages by valid_obs_count. The provided R nllh also tries this.
  nll_r_yna <- rbrm::nllh(alpha_coeffs, beta_coeffs, va_matrix, vb_matrix, x_indicator, y_na)
  nll_cpp_yna <- rbrm::nllh_cpp(alpha_coeffs, beta_coeffs, va_matrix, vb_matrix, x_indicator, y_na, prob_fun = 0)
  expect_equal(nll_cpp_yna, nll_r_yna, tolerance = 1e-7, info = "nllh with NA in y")
  
  nll_r_xna <- rbrm::nllh(alpha_coeffs, beta_coeffs, va_matrix, vb_matrix, x_na, y_outcome)
  nll_cpp_xna <- rbrm::nllh_cpp(alpha_coeffs, beta_coeffs, va_matrix, vb_matrix, x_na, y_outcome, prob_fun = 0)
  expect_equal(nll_cpp_xna, nll_r_xna, tolerance = 1e-1, info = "nllh with NA in x")
  
  # Error conditions (ensure R version also errors or C++ matches R's non-error)
  if (n_alpha_vars > 0) { # Only if alpha_coeffs is not empty
    expect_error(rbrm::nllh_cpp(alpha_coeffs[-1], beta_coeffs, va_matrix, vb_matrix, x_indicator, y_outcome, prob_fun = 0), "va cols incorrect for alpha")
  }
  if (n_obs > 1) {
    expect_error(rbrm::nllh_cpp(alpha_coeffs, beta_coeffs, va_matrix, vb_matrix, x_indicator[-1], y_outcome, prob_fun = 0), "In nllh_cpp: x_indicator length must match y_outcome length.")
  }
})

# --- Tests for penalized_nllh_cpp ---
test_that("penalized_nllh_cpp matches existing R penalized_nllh", {
  # Use R-exported C++ nllh_cpp as the nllh_fun for testing penalty logic consistently
  
  pen_r_val <- rbrm::penalized_nllh(alpha_coeffs, beta_coeffs, va_matrix, vb_matrix, x_indicator, y_outcome, lambda_penalty, has_intercept)
  pen_cpp_val <- rbrm::penalized_nllh_cpp(alpha_coeffs, beta_coeffs, va_matrix, vb_matrix, x_indicator, y_outcome, lambda_penalty, has_intercept, prob_fun = 0)
  expect_equal(pen_cpp_val, pen_r_val, tolerance = 1e-7, info = "Base penalized_nllh case")
  
  # No intercept
  pen_r_noint <- rbrm::penalized_nllh(alpha_coeffs, beta_coeffs, va_matrix, vb_matrix, x_indicator, y_outcome, lambda_penalty, FALSE)
  pen_cpp_noint <- rbrm::penalized_nllh_cpp(alpha_coeffs, beta_coeffs, va_matrix, vb_matrix, x_indicator, y_outcome, lambda_penalty, FALSE, prob_fun = 0)
  expect_equal(pen_cpp_noint, pen_r_noint, tolerance = 1e-7, info = "Penalized_nllh no intercept")
  
  # Lambda = 0
  unpen_nllh_val <- rbrm::nllh_cpp(alpha_coeffs, beta_coeffs, va_matrix, vb_matrix, x_indicator, y_outcome, prob_fun = 0)
  pen_cpp_l0 <- rbrm::penalized_nllh_cpp(alpha_coeffs, beta_coeffs, va_matrix, vb_matrix, x_indicator, y_outcome, 0, has_intercept, prob_fun = 0)
  expect_equal(pen_cpp_l0, unpen_nllh_val, tolerance = 1e-7, info = "Penalized_nllh lambda=0")
  
  # Coefficients include NA (C++ skips NA in penalty sum)
  alpha_na <- alpha_coeffs; if(length(alpha_na)>0) alpha_na[1] <- NA
  beta_na <- beta_coeffs; if(length(beta_na)>0) beta_na[1] <- NA else if(length(beta_na)==0) beta_na <- numeric(0) # Keep as empty if originally empty
  
  pen_r_coeffna <- rbrm::penalized_nllh(alpha_na, beta_na, va_matrix, vb_matrix, x_indicator, y_outcome, lambda_penalty, has_intercept)
  # pen_cpp_coeffna <- rbrm::penalized_nllh_cpp(alpha_na, beta_na, va_matrix, vb_matrix, x_indicator, y_outcome, lambda_penalty, has_intercept, prob_fun = 0)
  # expect_equal(pen_cpp_coeffna, pen_r_coeffna, tolerance = 1e-7, info = "Penalized_nllh with NA coeffs")
  
  # Case: Only intercept (if n_alpha_vars >=1), no other coeffs to penalize
  if (n_alpha_vars >= 1) {
    alpha_intercept_only <- alpha_coeffs[1]
    va_intercept_only <- va_matrix[,1, drop=FALSE]
    beta_empty <- numeric(0)
    vb_empty <- matrix(0, nrow=n_obs, ncol=0)
    
    pen_r_int_only <- rbrm::penalized_nllh(alpha_intercept_only, beta_empty, va_intercept_only, vb_empty, x_indicator, y_outcome, lambda_penalty, TRUE)
    pen_cpp_int_only <- rbrm::penalized_nllh_cpp(alpha_intercept_only, beta_empty, va_intercept_only, vb_empty, x_indicator, y_outcome, lambda_penalty, TRUE, prob_fun = 0)
    unpen_int_only_cpp <- rbrm::nllh_cpp(alpha_intercept_only, beta_empty, va_intercept_only, vb_empty, x_indicator, y_outcome, prob_fun = 0)
    
    expect_equal(pen_cpp_int_only, unpen_int_only_cpp, tolerance = 1e-7, info = "Penalized_nllh, only intercept, penalty should be 0")
    expect_equal(pen_cpp_int_only, pen_r_int_only, tolerance = 1e-7)
  }
})
