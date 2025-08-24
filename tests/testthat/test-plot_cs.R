# tests/testthat/test-plot_cs.R

testthat::test_that("Input validation: SCB", {
  expect_error(plot_cs(SCB = matrix(1,1,1), levels = 0, x = 1:3, mu_hat = 1:3),
    "`SCB` should be a list\\.") # SCB not a list
  expect_error(plot_cs(SCB = list(low = 1:3, up = 1:3), levels = 0, x = 1:3, mu_hat = 1:3),
    "`SCB` must have elements named 'scb_low' and 'scb_up'\\.") # SCB missing names
  expect_error(plot_cs(SCB = list(scb_low = 1:3, scb_up = 1:4), levels = 0, x = 1:3, mu_hat = 1:3),
    "Dimensions of `SCB\\$scb_up` and `SCB\\$scb_low` must match\\."
  ) # SCB dimension mismatch
  scb_chr <- list(scb_low = as.character(1:3), scb_up = as.character(1:3))
  expect_error(plot_cs(SCB = scb_chr, levels = 0, x = 1:3, mu_hat = 1:3),
    "Values of `SCB\\$scb_up` and `SCB\\$scb_low` must be numeric\\."
  ) # SCB non-numeric
})

testthat::test_that("Input validation: level", {
  scb <- list(scb_low = 1:5, scb_up = 2:6)
  expect_error(plot_cs(SCB = scb, levels = NULL, x = 1:5, mu_hat = 1:5),
    "Must provide input for `levels`\\." ) # must provide levels
  expect_error(plot_cs(SCB = scb, levels = letters[1:3], type = "upper", x = 1:5, mu_hat = 1:5),
    "Values of `levels` must be numeric\\.") # levels not numeric
  expect_error(plot_cs(SCB = scb, levels = matrix(1:4, 2), type = "lower", x = 1:5, mu_hat = 1:5),
    "`levels` should be a vector if `type` = upper or lower\\."
  ) # levels not vector

  # Unsupported types
  expect_error(plot_cs(SCB = scb, levels = 0, type = "interval", x = 1:5, mu_hat = 1:5),
    "'interval' is not avaliable for plotting")
  expect_error(plot_cs(SCB = scb, levels = 0, type = "two-sided", x = 1:5, mu_hat = 1:5),
    "'two-sided' is not avaliable for plotting")
  expect_error(plot_cs(SCB = scb, levels = 0, type = "oops", x = 1:5, mu_hat = 1:5),
    "`type` must be chosen between 'upper' and 'lower'\\.")
})

testthat::test_that("Input validation: x and y", {
  scb <- list(scb_low = 1:5, scb_up = 2:6)
  expect_error(plot_cs(SCB = scb, levels = 0, x = NULL, mu_hat = 1:5),
    "Must provide input for `x`\\.")

  n <- 50
  x <- seq(0, 1, length.out = n)
  scb <- list(
    scb_low = sin(seq(0, 2*pi, length.out = n)) - 0.5,
    scb_up  = sin(seq(0, 2*pi, length.out = n)) + 0.5
  )
  mu_hat <- sin(seq(0, 2*pi, length.out = n))
  expect_error(plot_cs(SCB = scb, levels = c(-0.5, 0, 0.5), x = as.list(x), mu_hat = mu_hat),
    "For 1D, `x` must be a vector\\.")# Non-vector x
  expect_error(plot_cs(SCB = scb, levels = 0, x = x[-1], mu_hat = mu_hat[-1]),
    "For 1D, `length\\(x\\)` must match length of `SCB\\$scb_up/scb_low`\\.") # x length mismatch
  expect_error(plot_cs(SCB = scb, levels = c(-0.5, 0, 0.5), x = as.factor(x), mu_hat = mu_hat),
               "`x` must be numeric/character\\.")# x is not numeric/character
})

testthat::test_that("Input validation: mu_hat / mu_true", {
  scb <- list(scb_low = 1:5, scb_up = 2:6)
  expect_error(plot_cs(SCB = scb, levels = 0, x = 1:5, mu_hat = letters[1:5]),
    "Input values of `mu_hat` must be numeric\\.")
  expect_error(plot_cs(SCB = scb, levels = 0, x = 1:5, mu_hat = 1:4),
    "Dimensions of `SCB\\$scb_up`, `SCB\\$scb_low` and `mu_hat` must match\\.")
  expect_error(plot_cs(SCB = scb, levels = 0, x = 1:5, mu_true = as.character(1:5)),
    "Input values of `mu_true` must be numeric\\.")
  expect_error(plot_cs(SCB = scb, levels = 0, x = 1:5, mu_true = 1:4),
    "Dimensions of `SCB\\$scb_up`, `SCB\\$scb_low` and `mu_true` must match\\.")
  expect_error(plot_cs(SCB = scb, levels = 0, x = 1:5),
    "An input must be provided for either `mu_hat` or `mu_true`\\.")
  expect_silent(plot_cs(SCB = scb, levels = 0, x = 1:5, mu_hat = 1:5, mu_true = 1:5))
})

testthat::test_that("palette / color label must be single strings", {
  scb <- list(scb_low = 1:5, scb_up = 2:6)
  expect_error(plot_cs(SCB = scb, levels = 0, x = 1:5, mu_hat = 1:5, palette = c("a","b")),
    "`palette` must be a single character string\\.")
  expect_error(plot_cs(SCB = scb, levels = 0, x = 1:5, mu_hat = 1:5, color_level_label = c("black","white")),
    "`color_level_label` must be a single character string\\.")
})

testthat::test_that("1D: continuous x requires numeric x (branch check)", {
  n <- 40
  x <- seq(0, 1, length.out = n)
  scb <- list(
    scb_low = rep(-0.2, n),
    scb_up  = rep( 0.8, n)
  )
  mu_hat <- seq(0, 1, length.out = n)

  # Valid numeric x -> ggplot
  p <- plot_cs(SCB = scb, levels = c(0.2, 0.5), type = "upper",
               x = x, mu_hat = mu_hat, xlab = "t", ylab = "f")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("2D branch: x/y lengths and together flag", {
  nx <- 12
  ny <- 10
  base <- outer(seq(-1, 1, length.out = nx), seq(-1, 1, length.out = ny),
                function(a,b) a^2 + b^2)
  scb <- list(scb_low = base - 0.1, scb_up = base + 0.1)
  mu_hat <- base
  x <- seq(0, 1, length.out = nx)
  y <- seq(0, 1, length.out = ny)

  expect_error(plot_cs(SCB = scb, levels = c(0.2, 0.6), type = "upper",
                       x = x[-1], y = y, mu_hat = mu_hat),
    "For 2D, `length\\(x\\)` must equal `nrow\\(SCB\\$scb_up\\)`\\.")# x length mismatch
  expect_error(plot_cs(SCB = scb, levels = c(0.2, 0.6), type = "upper",
                       x = x, y = y[-1], mu_hat = mu_hat),
    "For 2D, `length\\(y\\)` must equal `ncol\\(SCB\\$scb_up\\)`\\.")# y length mismatch
  expect_error(plot_cs(SCB = scb, levels = c(0.2, 0.6), type = "upper",
                       x = x, y = y, mu_hat = mu_hat[-1, -1]),
               "Dimensions of `SCB\\$scb_up`, `SCB\\$scb_low` and `mu_hat` must match\\.")
  expect_error(plot_cs(SCB = scb, levels = c(0.2, 0.6), type = "upper",
                       x = x, y = y, mu_hat = mu_hat, together = c(TRUE, FALSE)),
    "`together` must be a single logical value\\.") # together not logical

  # Valid call returns a patchwork object
  p <- plot_cs(SCB = scb, levels = c(0.2, 0.6, 1.0), type = "upper",
               x = x, y = y, mu_hat = mu_hat, together = TRUE)
  testthat::expect_s3_class(p, "patchwork")

  p2 <- plot_cs(SCB = scb, levels = c(0.2, 0.6, 1.0), type = "upper",
                x = x, y = y, mu_hat = mu_hat, together = FALSE)
  testthat::expect_s3_class(p2, "patchwork")
})

testthat::test_that("dims > 2 should error", {
  arr <- array(runif(2*3*4), dim = c(2,3,4))
  scb <- list(scb_low = arr - 0.1, scb_up = arr + 0.1)
  expect_error(plot_cs(SCB = scb, levels = 0, x = 1:2, mu_hat = arr[,,1]),
    "The dimension of `SCB\\$scb_up` and `SCB\\$scb_low` exceed 2\\.")
})

data(pupil)
testthat::test_that("Function works well", {
  library(mgcv)
  pupil_fpca <- prepare_pupil_fpca(pupil)
  fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
       s(seconds, by = use, k=30, bs = "cr") +
       s(seconds, by = age, k = 30, bs = "cr") +
       s(seconds, by = gender, k = 30, bs = "cr") +
       s(id, by = Phi1, bs="re") +
       s(id, by = Phi2, bs="re")+
       s(id, by = Phi3, bs="re") +
       s(id, by = Phi4, bs="re"),
       method = "fREML", data = pupil_fpca, discrete = TRUE)
  pupil_cma <- SCB_functional_outcome(data = pupil_fpca, object = fosr_mod, method = "cma",
                                      outcome = "percent_change",
                                      domain = "seconds", subset= c("use = 1"),
                                      id = "id")
  pupil_cma <- tibble::as_tibble(pupil_cma)
  p <- plot_cs(pupil_cma,levels = c(-18, -20, -22, -24), x = pupil_cma$domain,
          mu_hat = pupil_cma$mu_hat, xlab = "", ylab = "",
          level_label = T, min.size = 40, palette = "Spectral",
           color_level_label = "black")
  testthat::expect_s3_class(p, "ggplot")
})
