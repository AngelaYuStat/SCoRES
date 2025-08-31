set.seed(262)

x1 <- rnorm(100)
x2 <- rnorm(100)
epsilon <- rnorm(100,0,sqrt(2))
y <- -1 + x1 + 0.5 * x1^2 - 1.1 * x1^3 - 0.5 * x2 + 0.8 * x2^2 - 1.1 * x2^3 + epsilon
y <- SCoRES:::expit(y)
df <- data.frame(x1 = x1, x2 = x2, y = y)
grid <- data.frame(x1 = seq(-1, 1, length.out = 100), x2 = seq(-1, 1, length.out = 100))
# fit the logistic regression model and obtain the SCB for y
model <- "y ~ x1 + I(x1^2) + I(x1^3) + x2 + I(x2^2) + I(x2^3)"

test_that("Input validation works", {

  expect_error(SCB_logistic_outcome(NULL, model, grid), "`df_fit` must be provided")
  expect_error(SCB_logistic_outcome(df, NULL, grid), "`model` must be provided")

  expect_error(SCB_logistic_outcome(df, model, grid, n_boot = 0),
               "`n_boot` must be a positive integer")
  expect_error(SCB_logistic_outcome(df, model, grid, n_boot = 10.5),
               "`n_boot` must be a positive integer")
  expect_error(SCB_logistic_outcome(df, model, grid, alpha = 1),
               "`alpha` must be in \\(0, 1\\)")
  expect_error(SCB_logistic_outcome(df, model, grid, alpha = 0),
               "`alpha` must be in \\(0, 1\\)")

  expect_error(SCB_logistic_outcome(df, 123, grid),
               "`model` must be a formula or string")
})

test_that("Function works well: Output has correct shape and names when grid_df is provided", {

  results <- SCB_logistic_outcome(df_fit = df, model = model,
                                grid_df = grid, n_boot = 50, alpha = 0.1)

  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), nrow(grid))
  expect_true(all(c("scb_low", "Mean", "scb_up") %in% names(results)))

  expect_equal(length(results$scb_low), length(results$Mean))
  expect_equal(length(results$Mean), length(results$scb_up))
  expect_true(all(results$scb_low <= results$scb_up))
})

test_that("Function works well: Character predictors are coerced to factors and predictions succeed", {
  df2 <- df
  df2$grp <- ifelse(df2$x1 > 0, "A", "B")          # character
  grid2 <- data.frame(x1 = seq(-1, 1, length.out = 100),
                      x2 = seq(-1, 1, length.out = 100),
                      grp = factor(c("A", "B"))[1])
  model2 <- "y ~ x1 + x2 + grp"

  results <- SCB_logistic_outcome(df_fit = df2, model = model2,
                                grid_df = grid2, n_boot = 50, alpha = 0.1)
  expect_equal(nrow(results), nrow(grid2))

  expect_equal(length(results$scb_low), length(results$Mean))
  expect_equal(length(results$Mean), length(results$scb_up))
  expect_true(all(results$scb_low <= results$scb_up))
})

test_that("Function works well: y ~ 1 with no input for grid_df", {
  model3 <- "y ~ 1"
  results <- SCB_logistic_outcome(df_fit = df, model = model3,
                                grid_df = NULL, n_boot = 50, alpha = 0.1)

  expect_s3_class(results, "data.frame")
  expect_true(all(c("scb_low","Mean","scb_up") %in% names(results)))

  expect_equal(length(results$scb_low), length(results$Mean))
  expect_equal(length(results$Mean), length(results$scb_up))
  expect_true(all(results$scb_low <= results$scb_up))
})

test_that("Function works well: y ~ .", {
  model4 <- "y ~ ."
  expect_silent(SCB_logistic_outcome(df_fit = df, model = model4,
                                   grid_df = NULL, n_boot = 50, alpha = 0.1))
})

test_that("Function works well: y ~ x1 + x2, only provide values for x1 in grid_df", {
  grid4 <- data.frame(x1 = seq(-1, 1, length.out = 5))
  expect_silent(
    SCB_logistic_outcome(df_fit = df, model = model,
                       grid_df = grid4,
                       n_boot = 50, alpha = 0.1)
  )
})

test_that("Function works well: y ~ x1 + x2, no input for grid_df", {
  expect_silent(
    SCB_logistic_outcome(df_fit = df, model = model,
                       n_boot = 50, alpha = 0.1)
  )
})

test_that("Function works well: y ~ ., no input for grid_df, only y is contained in df_fit", {
  expect_silent(
    SCB_logistic_outcome(df_fit = data.frame(y = as.vector(df$y)), model = "y~.",
                       n_boot = 50, alpha = 0.1)
  )
})
