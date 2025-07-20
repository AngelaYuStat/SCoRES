test_that("mean_response_predict returns expected structure and values", {
  # Load ccds data
  data(ccds)
  ccds_fpca <- prepare_ccds_fpca(ccds)
  fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
                          s(seconds, by = use, k=30, bs = "cr") +
                          s(subject, by = Phi1, bs="re") +
                          s(subject, by = Phi2, bs="re")+
                          s(subject, by = Phi3, bs="re") +
                          s(subject, by = Phi4, bs="re"),
                        method = "fREML", data = ccds_fpca, discrete = TRUE)

  # Call the fuction
  results <- mean_response_predict(ccds_fpca, fosr_mod, time = "seconds",
                                   group_name = "use", group_value = 1, subject = "subject")

  # Check return value
  expect_type(results, "list")
  expect_named(results, c("s_pred", "pred_df", "lpmat", "mod_coef", "mod_cov"))

  # Check dimension
  expect_equal(nrow(results$pred_df), length(results$s_pred))
  expect_equal(ncol(results$lpmat), length(results$mod_coef))
  expect_equal(ncol(results$mod_cov), nrow(results$mod_cov))

  # Check the data type and dimension of mean and se
  expect_type(results$pred_df$mean, "double")
  expect_type(results$pred_df$se, "double")
  expect_equal(length(results$pred_df$mean), length(results$s_pred))

})

test_that("cma errors when required inputs are missing or invalid", {
  # Load ccds data
  data(ccds)
  ccds_fpca <- prepare_ccds_fpca(ccds)
  fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
                          s(seconds, by = use, k=30, bs = "cr") +
                          s(subject, by = Phi1, bs="re") +
                          s(subject, by = Phi2, bs="re")+
                          s(subject, by = Phi3, bs="re") +
                          s(subject, by = Phi4, bs="re"),
                        method = "fREML", data = ccds_fpca, discrete = TRUE)

  # Missing data
  expect_error(cma(NULL, fosr_mod, time = "seconds", group_name = "use", group_value = 1),
               "Must provide the origin data")

  # Missing object
  expect_error(cma(ccds_fpca, NULL, time = "seconds", group_name = "use", group_value = 1),
               "Must provide a fitted functional regression model")

  # Missing time
  expect_error(cma(ccds_fpca, fosr_mod, time = NULL, group_name = "use", group_value = 1),
               "Must provide the time variable name")

  # Dismatch group_name and group_value
  expect_error(cma(ccds_fpca, fosr_mod, time = "seconds", group_name = "use", group_value = 2),
               "not found in numeric variable")
})

test_that("cma returns correct structure and dimension", {
  set.seed(1)
  df <- data.frame(
    y = rnorm(50),
    time = rep(seq(0, 1, length.out = 10), 5),
    group = factor(rep(c(1, 0), each = 25)),
    id = rep(1:5, each = 10)
  )
  mod <- mgcv::bam(y ~ s(time, k = 5) + group, data = df)

  result <- cma(df, mod, time = "time", group_name = "group", group_value = 1, nboot = 100)

  expect_type(result, "list")
  expect_named(result, c("yhat", "time", "se_hat", "scb_low", "scb_up", "type"))

  # Check all returned vectors have same length
  n <- length(result$time)
  expect_equal(length(result$yhat), n)
  expect_equal(length(result$se_hat), n)
  expect_equal(length(result$scb_low), n)
  expect_equal(length(result$scb_up), n)

  # Check type field
  expect_match(result$type, "Confidence Interval")
})
