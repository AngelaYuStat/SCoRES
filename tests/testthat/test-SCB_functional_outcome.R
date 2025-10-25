# Load data
library(mgcv)
data(pupil)
pupil_fpca <- SCoRES::prepare_pupil_fpca(pupil)
fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
                        s(seconds, by = use, k=30, bs = "cr") +
                        s(id, by = Phi1, bs="re") +
                        s(id, by = Phi2, bs="re") +
                        s(id, by = Phi3, bs="re") +
                        s(id, by = Phi4, bs="re"),
                      method = "fREML", data = pupil_fpca, discrete = TRUE)
s_pred <- unique(pupil_fpca$seconds)

test_that("Input validation: basic checks", {
  expect_error(SCB_functional_outcome(NULL, method = "multiplier",
                                      outcome = "percent_change", domain = "seconds",
                                      id = "id"),
               "`data_df` must be provided\\.")

  # alpha
  expect_error(SCB_functional_outcome(pupil_fpca, method = "multiplier",
                                      alpha = 0, outcome = "percent_change",
                                      domain = "seconds", id = "id"),
               "`alpha` must be in \\(0, 1\\)\\.")

  # nboot
  expect_error(SCB_functional_outcome(pupil_fpca, method = "multiplier",
                                      nboot = -1, outcome = "percent_change",
                                      domain = "seconds", id = "id"),
               "`nboot` must be a positive integer\\.")
  expect_error(SCB_functional_outcome(pupil_fpca, method = "multiplier",
                                      nboot = 2.5, outcome = "percent_change",
                                      domain = "seconds", id = "id"),
               "`nboot` must be a positive integer\\.")

  # subset
  expect_error(SCB_functional_outcome(pupil_fpca, method = "multiplier",
                                      subset = 123, outcome = "percent_change",
                                      domain = "seconds", id = "id"),
               "`subset` should be character\\.")
})

test_that("method must be cma or multiplier", {
  expect_error(SCB_functional_outcome(pupil_fpca, method = "bootstrap",
                                      outcome = "percent_change",
                                      domain = "seconds", id = "id"),
               "`method` must be either 'cma' or 'multiplier'\\.")
})

test_that("Column existence checks for multiplier", {
  # missing domain
  df1 <- pupil_fpca
  names(df1)[names(df1) == "seconds"] <- "time"
  expect_error(
    SCB_functional_outcome(df1, method = "multiplier",
                           outcome = "percent_change", domain = "seconds", id = "id"),
    "`domain` should be a variable in `data_df`"
  )

  # missing id
  df2 <- pupil_fpca
  names(df2)[names(df2) == "id"] <- "subject"
  expect_error(
    SCB_functional_outcome(df2, method = "multiplier",
                           outcome = "percent_change", domain = "seconds", id = "id"),
    "`id` should be a variable in `data_df`"
  )

  # missing outcome
  df3 <- pupil_fpca
  names(df3)[names(df3) == "percent_change"] <- "y"
  expect_error(
    SCB_functional_outcome(df3, method = "multiplier",
                           outcome = "percent_change", domain = "seconds", id = "id"),
    "`outcome` should be a variable in `data_df`"
  )
})

test_that("CMA path errors when object missing", {
  expect_error(
    SCB_functional_outcome(
      data_df = pupil_fpca, object = NULL, method = "cma",
      fitted = TRUE, alpha = 0.05, outcome = "percent_change",
      domain = "seconds", id = "id"
    ),
    "`object` must be provided when `method = 'cma'`\\."
  )
})

# CMA Branch
test_that("CMA path works well", {
  skip_on_cran()
  results <- SCB_functional_outcome(
    data_df = pupil_fpca, object = fosr_mod, method = "cma",
    fitted = TRUE, alpha = 0.05, outcome = "percent_change",
    domain = "seconds", subset = c("use = 1"), id = "id")

  expect_type(results, "list")
  expect_named(results, c("mu_hat","domain","se_hat","scb_low","scb_up","type"))
  expect_identical(results$type, "CMA Confidence Band")

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

# Multiplier branch
test_that("Multiplier path works well", {
  skip_on_cran()
  results <- SCB_functional_outcome(
      data_df = pupil_fpca, object = fosr_mod, method = "multiplier",
      fitted = TRUE, alpha = 0.1, outcome = "percent_change",
      domain = "seconds", subset = c("use = 1"), id = "id")

  expect_type(results, "list")
  expect_named(results, c("mu_hat","domain","se_hat","scb_low","scb_up","type"))
  expect_identical(results$type, "Multiplier Bootstrap Confidence Interval")

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

test_that("Handles NAs with fpca.face imputation", {
  skip_on_cran()
  df_na <- pupil_fpca
  df_na$percent_change[df_na$seconds == 1 & df_na$id == "003-011"] <- NA

  results <- SCB_functional_outcome(
    data_df = df_na, object = NULL, method = "multiplier",
    fitted = TRUE, alpha = 0.05, outcome = "percent_change",
    domain = "seconds", id = "id")

  expect_type(results, "list")
  expect_named(results, c("mu_hat","se_hat","scb_low","scb_up","type","domain"))
  expect_false(any(is.na(results$mu_hat)))

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

test_that("Multiplier path without object: computes overall SCB and errors on all-zero rows", {
  skip_on_cran()

  library(dplyr)

  pupil_fpca2 <- pupil_fpca %>%
    group_by(id) %>%
    mutate(
      is_min = seconds == min(seconds),
      seconds = ifelse(is_min, 0, seconds),
      percent_change = ifelse(is_min, 0, percent_change)
    ) %>%
    ungroup() %>%
    select(-is_min)

  fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
                          s(seconds, by = use, k=30, bs = "cr") +
                          s(id, by = Phi1, bs="re") +
                          s(id, by = Phi2, bs="re") +
                          s(id, by = Phi3, bs="re") +
                          s(id, by = Phi4, bs="re"),
                        method = "fREML", data = pupil_fpca2, discrete = TRUE)
  s_pred <- unique(pupil_fpca2$seconds)

  results <- SCB_functional_outcome(
    data_df = pupil_fpca2, object = NULL, method = "multiplier",
    fitted = TRUE, alpha = 0.05, outcome = "percent_change",
    domain = "seconds", id = "id")
  expect_type(results, "list")
  expect_named(results, c("mu_hat","se_hat","scb_low","scb_up","type","domain"))
  expect_false(any(is.na(results$mu_hat)))

  expect_type(results$mu_hat, "double")
  expect_type(results$domain, "double")
  expect_type(results$se_hat, "double")
  expect_type(results$scb_low, "double")
  expect_type(results$scb_up, "double")
  expect_type(results$type, "character")

  expect_length(results$domain, length(s_pred) - 1) # remove all zeros entry on domain index zero

  n <- length(results$mu_hat)
  expect_length(results$domain, n)
  expect_length(results$se_hat, n)
  expect_length(results$scb_low, n)
  expect_length(results$scb_up, n)

  df_bad <- pupil_fpca
  df_bad$percent_change[df_bad$seconds == 2.674] <- 0

  expect_error(
    SCB_functional_outcome(
      data_df = df_bad, object = NULL, method = "multiplier",
      fitted = TRUE, alpha = 0.05, outcome = "percent_change",
      domain = "seconds", id = "id", nboot = 100
    ),
    "Expected no rows where all entries are zero in `data_df`\\."
  )
})
