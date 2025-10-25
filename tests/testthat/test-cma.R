# Load data
data(pupil)
library(mgcv)
pupil_fpca <- SCoRES::prepare_pupil_fpca(pupil, k_mean = 5, k_fpca = 5)
fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=5, bs="cr") +
                        s(seconds, by = use, k=5, bs = "cr") +
                        s(id, by = Phi1, bs="re"),
                      method = "fREML", data = pupil_fpca, discrete = TRUE)

test_that("Function works well", {

  # Call the function
  results <- cma(pupil_fpca, fosr_mod, outcome = "percent_change", domain = "seconds",
                                   subset = c("use = 1"), id = "id")

  expect_type(results, "list")

  expected_components <- c("mu_hat", "domain", "se_hat", "scb_low", "scb_up", "type")
  expect_named(results, expected_components)

  expect_type(results$mu_hat, "double")
  expect_type(results$domain, "double")
  expect_type(results$se_hat, "double")
  expect_type(results$scb_low, "double")
  expect_type(results$scb_up, "double")
  expect_type(results$type, "character")

  n <- length(results$mu_hat)
  expect_length(results$domain, n)
  expect_length(results$se_hat, n)
  expect_length(results$scb_low, n)
  expect_length(results$scb_up, n)

})
