testthat::test_that("Input validation: scb_low and scb_up", {
  expect_error(scb_to_cs(scb_low = 1:3, scb_up = 1:4, levels = c(0)),
               "Dimensions of `scb_up` and `scb_low` must match\\.")
  scb_chr <- list(scb_low = as.character(1:3), scb_up = as.character(1:3))
  expect_error(scb_to_cs(scb_low = as.character(1:3), scb_up = as.character(1:3), levels = c(0)),
               "Values of `scb_up` and `scb_low` must be numeric\\.")
})

testthat::test_that("Input validation: level", {
  expect_error(scb_to_cs(scb_low = 1:5, scb_up = 2:6, levels = NULL),
               "Must provide input for `levels`\\." ) # must provide levels
  expect_error(scb_to_cs(scb_low = 1:5, scb_up = 2:6, levels = letters[1:3]),
               "Values of `levels` must be numeric\\.") # levels not numeric
  expect_error(scb_to_cs(scb_low = 1:5, scb_up = 2:6, levels = matrix(1:4, 2), type = "lower"),
               "`levels` should be a vector if `type` = upper or lower\\."
  )
  expect_error(scb_to_cs(scb_low = 1:5, scb_up = 2:6, levels = matrix(1:4, 2), type = "interval"),
               "`levels` should be a list if `type` = interval") # levels not list with up and low
  expect_error(scb_to_cs(scb_low = 1:5, scb_up = 2:6, levels = list(up = 1:2, down = -1:0), type = "interval"),
               "`levels` must have elements named 'low' and 'up'\\.") # levels not list with up and low
})

testthat::test_that("Input validation: true_mean", {
  expect_error(scb_to_cs(scb_low = 1:5, scb_up = 2:6, levels = c(0), true_mean = letters[1:5]),
               "Values of `true_mean` must be numeric.")
  expect_error(scb_to_cs(scb_low = 1:5, scb_up = 2:6, levels = c(0), true_mean = 1:4),
               "Dimensions of `scb_up`, `scb_low` and `true_mean` must match.")
  nx <- 12
  ny <- 10
  base <- outer(seq(-1, 1, length.out = nx), seq(-1, 1, length.out = ny),
                function(a,b) a^2 + b^2)
  scb_low = base - 0.1
  scb_up = base + 0.1
  true_mean <- base

  expect_error(scb_to_cs(scb_low = scb_low, scb_up = scb_up, levels = c(0.2, 0.6), type = "upper",
                         true_mean = true_mean[-1, -1]),
               "Dimensions of `scb_up`, `scb_low` and `true_mean` must match\\.")
})

test_that("scb_to_cs returns correct structure.", {
  # Create test data
  scb_up <- matrix(2, nrow = 3, ncol = 3)
  scb_low <- matrix(-2, nrow = 3, ncol = 3)
  true_mean <- matrix(0, nrow = 3, ncol = 3)
  levels <- c(-1, 1)

  result <- scb_to_cs(scb_up, scb_low, levels, type = "upper")

  # Check basic structure
  expect_type(result, "list")
  expect_named(result, c("levels", "U_in", "U_out", "plot_cs"))

  # Check levels
  expect_equal(result$levels, levels)

  # Check U_in and U_out
  expect_type(result$U_in, "list")
  expect_type(result$U_out, "list")
  expect_length(result$U_in, length(levels))
  expect_length(result$U_out, length(levels))

  # Check all elements are logical matrices
  expect_true(all(sapply(result$U_in, is.matrix)))
  expect_true(all(sapply(result$U_out, is.matrix)))
  expect_true(all(sapply(result$U_in, is.logical)))
  expect_true(all(sapply(result$U_out, is.logical)))

  # Check dimensions match input
  expect_equal(dim(result$U_in[[1]]), dim(scb_up))
  expect_equal(dim(result$U_out[[1]]), dim(scb_up))

  # Check contain_* values (should be NULL when true_mean not provided)
  expect_null(result$contain_individual)
  expect_null(result$contain_all)

  # Check plot_cs (should be NULL when return_plot = FALSE)
  expect_null(result$plot_cs)

  result <- scb_to_cs(scb_up, scb_low, levels, type = "upper", true_mean = true_mean)

  # Should have containment results
  expect_type(result$contain_individual, "logical")
  expect_type(result$contain_all, "logical")
  expect_length(result$contain_individual, length(levels))

  result <- scb_to_cs(scb_up, scb_low, levels, type = "upper", true_mean = true_mean, return_contain_only = T)
  expect_null(result$U_in)
  expect_null(result$U_out)

  result <- scb_to_cs(scb_up, scb_low, levels, type = "upper", return_contain_only = T)
  expect_null(result$contain_individual)
  expect_null(result$contain_all)

  result <- scb_to_cs(scb_up, scb_low, levels, type = "two-sided")
  expect_type(result$L_out, "list")
  expect_length(result$L_out, length(levels))

  expect_true(all(sapply(result$L_out, is.matrix)))
  expect_true(all(sapply(result$L_out, is.logical)))

  results <- scb_to_cs(scb_up, scb_low, levels, type = "upper", return_plot = T, x1 = 1:3, x2 = 1:3, true_mean = true_mean)
  testthat::expect_type(results$plot_cs, "list")
})


