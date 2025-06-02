# # Assuming your package is 'rbrm' and it's loaded.
# # If not, use rbrm::getProbRR.org, rbrm::getProbRR_org_cpp etc.
# 
# # Helper for comparing results, especially handling NA/NaN consistently
# expect_results_equal <- function(res_r_wrapper, res_direct_cpp, tolerance = 1e-9) {
#   p0_wrapper <- res_r_wrapper$p0
#   p0_cpp   <- res_direct_cpp$p0
#   p1_wrapper <- res_r_wrapper$p1
#   p1_cpp   <- res_direct_cpp$p1
#   
#   # Convert NaNs to NAs for consistent comparison
#   p0_wrapper[is.nan(p0_wrapper)] <- NA
#   p0_cpp[is.nan(p0_cpp)]     <- NA
#   p1_wrapper[is.nan(p1_wrapper)] <- NA
#   p1_cpp[is.nan(p1_cpp)]     <- NA
#   
#   expect_equal(p0_wrapper, p0_cpp, tolerance = tolerance, label = "p0 values (R wrapper vs C++)")
#   expect_equal(p1_wrapper, p1_cpp, tolerance = tolerance, label = "p1 values (R wrapper vs C++)")
# }
# 
# # --- Test Scenarios (Input Data) ---
# logrr_scalar <- 0.5
# logop_scalar <- 0.3
# logrr_vec <- c(0.5, -0.2, 1.0, 0.0, NA, -13, 13, 0.5, -700, 0.2, 13, 0)
# logop_vec <- c(0.3, 0.1, -0.5, 0.0, 0.4, 0.2, 0.2, 1e-10, 0.1, -13, -13, 100)
# mat_in <- matrix(c(0.5, 0.3, -0.2, 0.1, 1.0, -0.5, NA, 0.4, -13, -13, -1000, -1000, 0, 100), ncol = 2, byrow = TRUE) # Added one more row for neg_disc case
# logrr_for_na_logop <- c(0.7, 0.8)
# empty_vec <- numeric(0)
# 
# # Test getProbRR.org R wrapper against its C++ backend
# test_that("R wrapper getProbRR.org matches its C++ backend getProbRR_org_cpp", {
#   # Scalar
#   expect_results_equal(
#     rbrm::getProbRR.org(logrr_scalar, logop_scalar), 
#     rbrm::getProbRR_org_cpp(logrr_scalar, logop_scalar)
#   )
#   expect_results_equal(
#     rbrm::getProbRR.org(logrr_scalar, logop_scalar, clipping = FALSE), 
#     rbrm::getProbRR_org_cpp(logrr_scalar, logop_scalar, clipping = FALSE)
#   )
#   
#   # Vector
#   # expect_results_equal(
#   #   rbrm::getProbRR.org(logrr_vec, logop_vec), 
#   #   rbrm::getProbRR_org_cpp(logrr_vec, logop_vec)
#   # )
#   
#   # Scalar NA inputs
#   # expect_results_equal(rbrm::getProbRR.org(NA, logop_scalar), rbrm::getProbRR_org_cpp(NA_real_, logop_scalar))
#   # expect_results_equal(rbrm::getProbRR.org(logrr_scalar, NA), rbrm::getProbRR_org_cpp(logrr_scalar, NA_real_))
#   # expect_results_equal(rbrm::getProbRR.org(NA, NA), rbrm::getProbRR_org_cpp(NA_real_, NA_real_))
#   
#   # Recycling
#   # expect_results_equal(rbrm::getProbRR.org(logrr_scalar, logop_vec), rbrm::getProbRR_org_cpp(logrr_scalar, logop_vec))
#   
#   
#   # expect_results_equal(rbrm::getProbRR.org(logrr_vec, logop_scalar), rbrm::getProbRR_org_cpp(logrr_vec, logop_scalar))
#   
#   # Recycling with scalar NA
#   
#   # expect_results_equal(rbrm::getProbRR.org(rep(NA, length(logrr_vec)), NA), rbrm::getProbRR_org_cpp(logrr_vec, NA_real_))
#   # expect_results_equal(rbrm::getProbRR.org(NA, logop_vec), rbrm::getProbRR_org_cpp(NA_real_, logop_vec))
#   # 
#   # Empty vectors
#   expect_results_equal(rbrm::getProbRR.org(empty_vec, empty_vec), rbrm::getProbRR_org_cpp(empty_vec, empty_vec))
#   # expect_results_equal(rbrm::getProbRR.org(empty_vec, logop_scalar), rbrm::getProbRR_org_cpp(empty_vec, logop_scalar))
#   # expect_results_equal(rbrm::getProbRR.org(logrr_scalar, empty_vec), rbrm::getProbRR_org_cpp(logrr_scalar, empty_vec))
# })
# 
# 
# # Test getProbRR.alt R wrapper against its C++ backend
# test_that("R wrapper getProbRR.alt matches its C++ backend getProbRR_alt_cpp", {
#   # Scalar
#   expect_results_equal(
#     rbrm::getProbRR.alt(logrr_scalar, logop_scalar), 
#     rbrm::getProbRR_alt_cpp(logrr_scalar, logop_scalar)
#   )
#   expect_results_equal(
#     rbrm::getProbRR.alt(logrr_scalar, logop_scalar, clipping = FALSE), 
#     rbrm::getProbRR_alt_cpp(logrr_scalar, logop_scalar, clipping = FALSE)
#   )
#   
#   # Vector
#   expect_results_equal(
#   rbrm::getProbRR.alt(logrr_vec, logop_vec),
#   rbrm::getProbRR_alt_cpp(logrr_vec, logop_vec)
#   )
#   
#   # logrr 2-element vector, logop = NA
#   # expect_results_equal(
#   #   rbrm::getProbRR.alt(logrr_for_na_logop, NA), 
#   #   rbrm::getProbRR_alt_cpp(logrr_for_na_logop, NA_real_)
#   # )
# 
#   # Scalar NA inputs
#   # expect_results_equal(rbrm::getProbRR.alt(NA, logop_scalar), rbrm::getProbRR_alt_cpp(NA_real_, logop_scalar))
#   # expect_results_equal(rbrm::getProbRR.alt(logrr_scalar, NA), rbrm::getProbRR_alt_cpp(logrr_scalar, NA_real_))
#   # expect_results_equal(rbrm::getProbRR.alt(NA,NA), rbrm::getProbRR_alt_cpp(NA_real_,NA_real_))
#   
#   # Recycling
#   # expect_results_equal(rbrm::getProbRR.alt(logrr_scalar, logop_vec), rbrm::getProbRR_alt_cpp(logrr_scalar, logop_vec))
#   # expect_results_equal(rbrm::getProbRR.alt(logrr_vec, logop_scalar), rbrm::getProbRR_alt_cpp(logrr_vec, logop_scalar))
#   
#   # Recycling with scalar NA
#   # expect_results_equal(rbrm::getProbRR.alt(logrr_vec, NA), rbrm::getProbRR_alt_cpp(logrr_vec, NA_real_))
#   # expect_results_equal(rbrm::getProbRR.alt(NA, logop_vec), rbrm::getProbRR_alt_cpp(NA_real_, logop_vec))
#   
#   # Empty vectors
#   expect_results_equal(rbrm::getProbRR.alt(empty_vec, empty_vec), rbrm::getProbRR_alt_cpp(empty_vec, empty_vec))
#   # expect_results_equal(rbrm::getProbRR.alt(empty_vec, logop_scalar), rbrm::getProbRR_alt_cpp(empty_vec, logop_scalar))
#   # expect_results_equal(rbrm::getProbRR.alt(logrr_scalar, empty_vec), rbrm::getProbRR_alt_cpp(logrr_scalar, empty_vec))
#   
#   # Negative discriminant case for .alt (C++ should give NA and Rcout warning)
#   # The R wrapper will just reflect this NA.
#   # suppressWarnings might be needed if the Rcpp::Rcout printouts are treated as warnings by some test runners
#   res_r_wrapper_neg_disc <- rbrm::getProbRR.alt(0, 100) 
#   res_cpp_direct_neg_disc <- rbrm::getProbRR_alt_cpp(0, 100)
#   expect_results_equal(res_r_wrapper_neg_disc, res_cpp_direct_neg_disc)
# })
# 
# # Error handling tests: ensure R wrappers propagate C++ errors or have their own.
# test_that("Error handling for mismatched lengths (non-recyclable) in R wrappers", {
#   # These errors come from the C++ code if wrappers are thin.
#   expect_error(rbrm::getProbRR_org_cpp(c(1,2,3), c(1,2)),
#                "In prob_rr_org_worker: logrr and logop must have the same length.") 
#   expect_error(rbrm::getProbRR_alt_cpp(c(1,2,3), c(1,2)),
#                "In prob_rr_alt_worker: logrr and logop must have the same length.") 
# })

