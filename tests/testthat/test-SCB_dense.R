test_that("Function works well", {
  set.seed(123)
  library(SCoRES)
  A <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)

  results <- SCoRES:::SCB_dense(A)
  expect_type(results, "list")

  expected_components <- c("mu_hat", "se_hat", "scb_low", "scb_up", "type")
  expect_named(results, expected_components)

  expect_type(results$mu_hat, "double")
  expect_type(results$se_hat, "double")
  expect_type(results$scb_low, "double")
  expect_type(results$scb_up, "double")
  expect_type(results$type, "character")

  n <- length(results$mu_hat)
  expect_length(results$se_hat, n)
  expect_length(results$scb_low, n)
  expect_length(results$scb_up, n)
})
