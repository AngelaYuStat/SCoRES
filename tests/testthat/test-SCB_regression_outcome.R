set.seed(262)

# generate simulated data
x1 <- rnorm(100)
epsilon <- rnorm(100,0,sqrt(2))
y <- -1 + x1 + epsilon
df <- data.frame(x1 = x1, y = y)
grid <- data.frame(x1 = seq(-1, 1, length.out = 100))
# fit the linear regression model and obtain the SCB for y
model <- "y ~ x1"
w <- matrix(c(0, 1), ncol = 2)

test_that("Function works well", {

  results <- SCB_regression_outcome(df_fit = df, model = model,
                                grid_df = grid, n_boot = 50, alpha = 0.1,
                                fitted = FALSE, w = w)

  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), nrow(w))
  expect_true(all(c("scb_low", "Mean", "scb_up") %in% names(results)))

  expect_equal(length(results$scb_low), length(results$Mean))
  expect_equal(length(results$Mean), length(results$scb_up))
  expect_true(all(results$scb_low <= results$scb_up))
})
