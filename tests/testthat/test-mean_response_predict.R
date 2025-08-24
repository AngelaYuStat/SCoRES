# tests/testthat/test-mean_response_predict.R

set.seed(123)

library(tibble)
library(mgcv)
library(dplyr)

data(pupil)

testthat::test_that("Function works well.", {
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

  results <- mean_response_predict(pupil_fpca, fosr_mod, fitted = TRUE,
             outcome = "percent_change", domain = "seconds", subset = c("use = 1"), id = "id")

  testthat::expect_true(is.list(results))
  testthat::expect_true(all(c("s_pred","pred_df","lpmat","mod_coef","mod_cov") %in% names(results)))

  testthat::expect_equal(nrow(results$pred_df), length(results$s_pred))
  testthat::expect_true(all(c("mean","se") %in% names(results$pred_df)))

  testthat::expect_equal(nrow(results$lpmat), length(results$s_pred))
  testthat::expect_equal(length(results$mod_coef), ncol(results$lpmat))
  testthat::expect_equal(dim(results$mod_cov), c(ncol(results$lpmat), ncol(results$lpmat)))
})

testthat::test_that("Input validation: subset", {
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

  testthat::expect_error(
    mean_response_predict(pupil_fpca, fosr_mod, fitted = TRUE,
                                     outcome = "percent_change", domain = "seconds", subset = c("height = 170"), id = "id"),
    "Variable 'height' not found in `data_df`\\."
  )

  testthat::expect_error(
    mean_response_predict(pupil_fpca, fosr_mod, fitted = TRUE,
                          outcome = "percent_change", domain = "seconds", subset = c("use = m"), id = "id"),
    "is not numeric"
  )
  testthat::expect_warning(
    testthat::expect_error(
      mean_response_predict(pupil_fpca, fosr_mod, fitted = TRUE,
                            outcome = "percent_change", domain = "seconds", subset = c("use = 1", "badpattern"),
        id = "id"
      )
    ),
    "did not match the required <name> = <value> pattern"
  )

  pupil_fpca <- pupil_fpca %>%
    mutate(use = as.character(use))

  testthat::expect_error(
    mean_response_predict(pupil_fpca, fosr_mod, fitted = TRUE,
                          outcome = "percent_change", domain = "seconds", subset = c("use = 1"), id = "id"),
    "is of type character"
  )
  pupil_fpca <- pupil_fpca %>%
    mutate(use = as.factor(use))
  testthat::expect_error(
    mean_response_predict(pupil_fpca, fosr_mod, fitted = TRUE,
                          outcome = "percent_change", domain = "seconds", subset = c("use = 1"), id = "id"),
    "is of type factor"
  )
})
