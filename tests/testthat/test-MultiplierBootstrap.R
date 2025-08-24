test_that("Function works well", {
  set.seed(123)
  library(SCoRES)
  R <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)

  result <- SCoRES:::MultiplierBootstrap(R)
  expect_type(result, "list")

  expected_components <- c("z", "q", "samples")
  expect_named(result, expected_components)

  expect_type(result$z, "double")
  expect_type(result$q, "double")
  expect_type(result$samples, "double")
})
