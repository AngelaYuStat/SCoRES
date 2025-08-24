# tests/testthat/test-SCB_regression_coef.R

set.seed(123)

testthat::test_that("Input validation: df_fit and model", {
  expect_error(SCB_regression_coef(df_fit = NULL, model = y ~ x),
               "`df_fit` must be provided\\.")
  expect_error(SCB_regression_coef(df_fit = data.frame(x = 1:3, y = 1:3), model = NULL),
               "`model` must be provided\\.")
  m <- matrix(rnorm(40), nrow = 20, ncol = 2)
  colnames(m) <- c("x", "y")
  expect_silent(
    SCB_regression_coef(df_fit = m, model = y ~ x)
  ) # can be coverted to data.frame
  df <- data.frame(y = rnorm(20), x1 = rnorm(20))
  expect_error(
    SCB_regression_coef(df_fit = df, model = y ~ x1 + x2, n_boot = 10),
    "`df_fit` is missing variables: x2"
  )
  df <- data.frame(y = rnorm(30), x1 = rnorm(30), x2 = rnorm(30))
  expect_silent(SCB_regression_coef(df_fit = df, model = y ~ ., n_boot = 10))
  expect_error(SCB_regression_coef(df_fit = data.frame(x = 1:3, y = 1:3), model = 1),
               "`model` must be a formula or string\\.")
})

testthat::test_that("Input validation: n_boot, type and alpha", {
  df <- data.frame(y = rnorm(30), x = rnorm(30))

  expect_error(SCB_regression_coef(df_fit = df, model = y ~ x, n_boot = 0),
               "`n_boot` must be a positive integer\\.")
  expect_error(SCB_regression_coef(df_fit = df, model = y ~ x, n_boot = 10.5),
               "`n_boot` must be a positive integer\\.")

  expect_error(SCB_regression_coef(df_fit = df, model = y ~ x, alpha = 0),
               "`alpha` must be in \\(0, 1\\)\\.")
  expect_error(SCB_regression_coef(df_fit = df, model = y ~ x, alpha = 1),
               "`alpha` must be in \\(0, 1\\)\\.")
  expect_error(SCB_regression_coef(df_fit = df, model = y ~ x, alpha = -0.1),
               "`alpha` must be in \\(0, 1\\)\\.")

  expect_error(SCB_regression_coef(df_fit = df, model = y ~ x, type = "poisson"),
               "`type` must be either 'linear' or 'logistic'\\.")
})

testthat::test_that("Function works well: linear branch.", {
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + 2*x1 - 3*x2 + rnorm(n, sd = 0.8)
  df <- data.frame(y, x1, x2)

  res <- SCB_regression_coef(df_fit = df, model = y ~ x1 + x2, n_boot = 150, alpha = 0.1, type = "linear")

  expect_s3_class(res, "data.frame")
  expect_true(all(c("scb_low", "Mean", "scb_up") %in% names(res)))

  fit <- lm(y ~ x1 + x2, df)
  expect_identical(rownames(res), names(coef(fit)))
  expect_true(all(res$scb_low <= res$scb_up))
  expect_true(all(res$Mean >= res$scb_low & res$Mean <= res$scb_up))

  df <- data.frame(y = rnorm(80))
  # y ~ ., only intercept
  res <- SCB_regression_coef(df_fit = df, model = y ~ ., n_boot = 80, type = "linear")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_identical(rownames(res), "(Intercept)")
  expect_true(all(res$scb_low <= res$scb_up))
})

testthat::test_that("Function works well: logistic branch.", {
  n <- 300
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  eta <- -0.5 + 1.2*x1 - 0.8*x2
  p <- 1/(1 + exp(-eta))
  y <- rbinom(n, 1, p)
  df <- data.frame(y, x1, x2)

  res <- SCB_regression_coef(df_fit = df, model = y ~ x1 + x2, n_boot = 100, alpha = 0.1, type = "logistic")

  expect_s3_class(res, "data.frame")
  expect_true(all(c("scb_low", "Mean", "scb_up") %in% names(res)))

  fit <- suppressWarnings(glm(y ~ x1 + x2, data = df, family = binomial()))
  expect_identical(rownames(res), names(coef(fit)))

  expect_true(all(res$scb_low <= res$scb_up))
})

