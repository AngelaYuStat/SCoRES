#' Construct Simultaneous Confidence Bands for a Linear Model
#'
#' This function fits a linear model and constructs simultaneous confidence bands (SCB)
#' using a non-parametric bootstrap method for the mean outcome of regression on a fixed test set design matrix
#'
#' @param df_fit A data frame containing the training design matrix used to fit the linear model.
#' @param model A character string representing the formula for the linear model (e.g., \code{"y ~ x1 + x2"}).
#' @param grid_df A data frame containing the test set design matrix at which to evaluate the predictions and construct confidence bands.
#' @param n_boot Number of bootstrap samples used in the non-parametric bootstrap procedure to generate the empirical distribution. Default is 1000.
#' @param alpha Significance level for the confidence band (e.g., 0.05 for 95% confidence). Default is 0.05.
#' @param grid_df_boot An optional data frame specifying the input grid at which predictions are evaluated during bootstrap resampling.
#'                     This allows SCBs to be constructed on a denser or alternative set of covariate values if desired.. If NULL, uses \code{grid_df}.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{LowerBound}{Lower bound of the simultaneous confidence band.}
#'   \item{Mean}{Predicted mean response from the fitted model.}
#'   \item{UpperBound}{Upper bound of the simultaneous confidence band.}
#'   \item{...}{All columns from \code{grid_df}, representing the prediction grid.}
#' }
#'
#' @import stats
#' @export
#'
#' @examples
#' set.seed(262)
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' y <- -1 + x1 + 0.5 * x1^2 - 1.1 * x1^3 - 0.5 * x2 + 0.8 * x2^2 - 1.1 * x2^3 - rnorm(100,0,sqrt(2))
#' df <- data.frame(x1 = x1, x2 = x2, y = y)
#' grid <- data.frame(x1 = seq(-1, 1, length.out = 100), x2 = seq(-1, 1, length.out = 100))
#' model <- "y ~ x1 + I(x1^2) + I(x1^3) + x2 + I(x2^2) + I(x2^3)"
#' SCB_linear_outcome(df_fit = df, model = model, grid_df = grid)
#'
SCB_linear_outcome = function(df_fit, model, grid_df, n_boot = 1000, alpha = 0.05, grid_df_boot = NULL){
  formula_ = as.formula(model)
  fit = lm(model, df_fit)
  y_hat = predict(fit, grid_df, se.fit = TRUE, level = 1 - alpha) # for constructing the whole
  res_max_v = rep(0,n_boot)
  if(is.null(grid_df_boot)){
    grid_df_boot = grid_df
    y_hat_level_grid = y_hat
  }else{
    y_hat_level_grid = predict(fit, grid_df_boot, se.fit = TRUE, level = 1 - alpha)# bootstrap true mean
  }
  for(i in 1:n_boot){
    df_boot = df_fit[sample(1:dim(df_fit)[1], replace = T),]
    fit_boot = lm(model, df_boot)
    y_hat_boot = predict(fit_boot, grid_df_boot, level = 1 - alpha, se.fit = TRUE	)
    residual = abs(y_hat_boot$fit - y_hat_level_grid$fit)/y_hat_boot$se.fit
    res_max_v[i] = max(residual)
  }
  thres = quantile(res_max_v, probs = 1 - alpha)
  sim_CB = data.frame(LowerBound = y_hat$fit - thres*y_hat$se.fit, Mean = y_hat$fit, UpperBound = y_hat$fit + thres*y_hat$se.fit, grid_df)
}

#' Expit (inverse logit) function
#'
#' Computes the inverse logit transformation.
#'
#' @param x A numeric input.
#' @return Value between 0 and 1.
#' @export
#'
#' @examples
#' expit(0)         # returns 0.5
#' expit(c(-2, 0, 2))
expit = function(x){
  1/(1+exp(-x))
}


#' Construct Simultaneous Confidence Bands for a Logistic Regression Model
#'
#' This function fits a logistic regression model and constructs simultaneous confidence bands (SCB)
#' using a non-parametric bootstrap method for the mean outcome of regression on a fixed test set design matrix
#'
#' @param df_fit A data frame containing the training design matrix used to fit the logistic model.
#' @param model A character string representing the formula for the logistic model (e.g., \code{"y ~ x1 + x2"}).
#' @param grid_df A data frame containing the test set design matrix at which to evaluate the predictions and construct confidence bands.
#' @param n_boot Number of bootstrap samples used in the non-parametric bootstrap procedure to generate the empirical distribution. Default is 1000.
#' @param alpha Significance level for the confidence band (e.g., 0.05 for 95% confidence). Default is 0.05.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{LowerBound}{Lower bound of the simultaneous confidence band.}
#'   \item{Mean}{Predicted mean response from the fitted model.}
#'   \item{UpperBound}{Upper bound of the simultaneous confidence band.}
#'   \item{grid_df}{All columns from \code{grid_df}, representing the prediction grid.}
#' }
#'
#' @import stats
#' @export
#'
#' @examples
#' set.seed(262)
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' y <- -1 + x1 + 0.5 * x1^2 - 1.1 * x1^3 - 0.5 * x2 + 0.8 * x2^2 - 1.1 * x2^3 - rnorm(100,0,sqrt(2))
#' y <- expit(y)
#' df <- data.frame(x1 = x1, x2 = x2, y = y)
#' grid <- data.frame(x1 = seq(-1, 1, length.out = 100), x2 = seq(-1, 1, length.out = 100))
#' model <- "y ~ x1 + I(x1^2) + I(x1^3) + x2 + I(x2^2) + I(x2^3)"
#' SCB_logistic_outcome(df_fit = df, model = model, grid_df = grid)
#'
SCB_logistic_outcome = function(df_fit, model, grid_df, n_boot = 1000, alpha = 0.05){
  formula_ = as.formula(model)
  fit =  suppressWarnings(glm(model, family = binomial(), data = df_fit)) # Suppress warning forfitted probabilities numerically 0 or 1
  y_hat = predict(fit, grid_df, se.fit = TRUE, level = 1 - alpha) # bootstrap true mean
  res_max_v = rep(0,n_boot)
  for(i in 1:n_boot){
    df_boot = df_fit[sample(1:dim(df_fit)[1], replace = T),]
    fit_boot = suppressWarnings(glm(model, family = binomial(), data = df_boot))
    y_hat_boot = predict(fit_boot, grid_df, level = 1 - alpha, se.fit = TRUE	)
    residual = abs(y_hat_boot$fit - y_hat$fit)/y_hat_boot$se.fit
    res_max_v[i] = max(residual)
  }
  thres = quantile(res_max_v, probs = 1 - alpha)
  sim_CB = data.frame(LowerBound = expit(y_hat$fit - thres*y_hat$se.fit), Mean = expit(y_hat$fit),
                      UpperBound = expit(y_hat$fit + thres*y_hat$se.fit), grid_df)
}

#' Construct Simultaneous Confidence Bands for Regression Coefficients
#'
#' This function fits either a linear or logistic regression model and computes simultaneous confidence bands (SCBs)
#' for the model coefficients using a non-parametric bootstrap procedure.
#'
#' @param df_fit A data frame containing the design matrix and response variable used to fit the model.
#' @param model A character string specifying the regression formula (e.g., \code{"y ~ x1 + x2"}).
#' @param n_boot Integer. Number of bootstrap samples to use for constructing the SCBs. Default is 5000.
#' @param alpha Numeric. Significance level for the confidence bands (e.g., 0.05 for 95% SCBs). Default is 0.05.
#' @param type A character string specifying the model type. Either \code{"linear"} (default) or \code{"logistic"}.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{LowerBound}{Lower bound of the simultaneous confidence band for each coefficient.}
#'   \item{Mean}{Estimated value of each coefficient from the fitted model.}
#'   \item{UpperBound}{Upper bound of the simultaneous confidence band for each coefficient.}
#'   \item{Row names}{Correspond to the model coefficients.}
#' }
#'
#' @import stats
#' @export
#'
#' @examples
#' set.seed(262)
#' M <- 50
#' rho <- 0.4
#' n <- 500
#' beta <- rnorm(M, mean = 0, sd = 1)
#' Sigma <- outer(1:M, 1:M, function(i, j) rho^abs(i - j))
#' X <- MASS::mvrnorm(n = n, mu = rep(0, M), Sigma = Sigma)
#' epsilon <- rnorm(n, mean = 0, sd = 1)
#' y <- X %*% beta + epsilon
#' df <- as.data.frame(X)
#' names(df) <- paste0("x", 1:M)
#' df$y <- as.vector(y)
#' model <- "y ~ ."
#' SCB_regression_coef(df, model)
#'
SCB_regression_coef = function(df_fit, model, n_boot = 5000, alpha = 0.05, type = "linear"){
  # type: "linear" or "logistic"
  formula_ = as.formula(model)
  if(type == "linear"){
    fit = lm(model, df_fit)
  }else if(type == "logistic"){
    fit = suppressWarnings(glm(model, family = binomial(), data = df_fit)) # Suppress warning for 0 or 1 probability for simulation
  }
  sum_out = summary(fit)
  coef_hat = fit$coefficients # bootstrap true mean
  coef_sd = sum_out$coefficients[,2]
  res_max_v = rep(0,n_boot)
  for(i in 1:n_boot){
    df_boot = df_fit[sample(1:dim(df_fit)[1], replace = T),]
    if(type == "linear"){
      fit_boot = lm(model, df_boot)
    }else if(type == "logistic"){
      fit_boot = suppressWarnings(glm(model, family = binomial(), data = df_boot))
    }
    sum_out_boot = summary(fit_boot)
    coef_hat_boot =fit_boot$coefficients
    coef_sd_boot = sum_out_boot$coefficients[,2]
    residual = abs(coef_hat_boot - coef_hat)/coef_sd_boot
    res_max_v[i] = max(residual)
  }
  thres = quantile(res_max_v, probs = 1 - alpha)
  sim_CB = data.frame(LowerBound = coef_hat - thres*coef_sd, Mean = coef_hat, UpperBound = coef_hat + thres*coef_sd)
}
