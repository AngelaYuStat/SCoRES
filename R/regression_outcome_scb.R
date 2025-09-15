#' Simultaneous Confidence Bands for Regression Outcomes or Coefficients
#'
#' This function fits a user-specified regression model (linear or logistic),
#' and constructs simultaneous confidence bands (SCBs) for either the mean outcome
#' or the regression coefficients. SCBs are obtained using a nonparametric bootstrap
#' procedure evaluated on a fixed test design matrix, providing simultaneous
#' inference across the entire range of covariates or parameters of interest.
#'
#' @param df_fit A data frame containing the training design matrix used
#' to fit the linear model. Acceptable input format includes numeric and factor.
#' @param model A character string representing the formula for the linear model
#' (e.g., \code{"y ~ x1 + x2"}).
#' @param grid_df A data frame specifying the covariate settings that define the
#' mean outcome or linear combination for which simultaneous confidence bands (SCB)
#' are constructed.
#' Each row represents one covariate combination at which predictions and
#' SCBs are evaluated. Column names should match variables in the fitted model,
#' but `grid_df` may include only the subset of covariates of interest for the SCB
#' (it is not required to cover all model variables).
#' Default is `NULL`, in which case the SCB is constructed over the fitted values
#' based on 'df_fit`.
#' @param n_boot Number of bootstrap samples used in the non-parametric bootstrap
#' procedure to generate the empirical distribution. Default is 1000.
#' @param alpha Significance level for the confidence band
#' (e.g., 0.05 for 95% confidence). Default is 0.05.
#' @param grid_df_boot An optional data frame specifying the input grid at which
#' predictions are evaluated during bootstrap resampling.
#' This allows SCBs to be constructed on a denser set of covariate values
#' if desired. If NULL, uses \code{grid_df}.
#' If `grid_df` is set to NULL, `grid_df_boot` will also be set to `NULL`.
#' This argument is only for `type` = `linear`.
#' @param type A character string specifying the model type. Either \code{"linear"}
#' (default) or \code{"logistic"}.
#' @param fitted Logical. Whether to estimate the simultaneous confidence bands
#' for regression outcomes or coefficients.
#'   \itemize{
#'     \item \code{TRUE} - Estimate the simultaneous confidence bands
#'     for regression outcomes.
#'     \item \code{FALSE} - estimate the simultaneous confidence bands
#'     for regression coefficients.
#'     }
#'   Default is \code{TRUE}.
#' @param w A numeric matrix that specifies the linear combinations of regression coefficients
#' for which simultaneous confidence bands (SCBs) are to be constructed.
#' The number of columns should be equal to the number of coefficients in the
#' regression model fitted.
#' Default is NULL, will return SCBs for all coefficients and the intercept.
#' This argument is only for `fitted` = `FALSE`.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{scb_low}{Lower bound of the simultaneous confidence band.}
#'   \item{Mean}{Predicted mean response from the fitted model.}
#'   \item{scb_up}{Upper bound of the simultaneous confidence band.}
#'   \item{...}{All columns from \code{grid_df}, representing the prediction grid.
#'              Optional, if `fitted` = \code{TRUE}}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(262)
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' epsilon <- rnorm(100,0,sqrt(2))
#' y <- -1 + x1 + 0.5 * x1^2 - 1.1 * x1^3 - 0.5 * x2 + 0.8 * x2^2 - 1.1 * x2^3 + epsilon
#' df <- data.frame(x1 = x1, x2 = x2, y = y)
#' grid <- data.frame(x1 = seq(-1, 1, length.out = 100),
#'                    x2 = seq(-1, 1, length.out = 100))
#' model <- "y ~ x1 + I(x1^2) + I(x1^3) + x2 + I(x2^2) + I(x2^3)"
#' w <- matrix(c(0, 0, 1, 0, 0, 1, 0), ncol = 7)
#' results <- SCB_regression_outcome(df_fit = df, model = model, grid_df = grid, w   = w)
#'
SCB_regression_outcome = function(df_fit, model, grid_df = NULL, n_boot = 1000,
                                  alpha = 0.05, grid_df_boot = NULL, type = "linear",
                                  fitted = TRUE, w = NULL){
  if(fitted == TRUE){
    if(type == "linear"){
      return(SCB_linear_outcome(df_fit, model, grid_df, n_boot,
                                alpha, grid_df_boot))
    }else if(type == "logistic"){
      return(SCB_logistic_outcome(df_fit, model, grid_df, n_boot,
                                alpha))
    }else{
      stop("`type` must be either 'linear' or 'logistic'.")
    }
  }else{
    SCB <- SCB_regression_coef(df_fit, model, n_boot, alpha, type)
    if (is.null(w)) {
      return(SCB)
    }else{
      if(!(is.matrix(w)||is.array(w)) || !is.numeric(w)){
        stop("`w` should be a numeric matrix or array.")
      }
      if(nrow(SCB) != ncol(w)) {
        stop("The number of columns in `w` must be equal to the number of
             the coefficients in the model.")
      }else{
        SCB_coefs <- as.matrix(SCB)
        results <- w %*% SCB_coefs
        return(data.frame(scb_low = results[, 1], Mean = results[, 2], scb_up = results[, 3]))
      }
    }

  }

}
#' Construct Simultaneous Confidence Bands for Linear Regression Outcome
#'
#' This function fits a linear model and constructs simultaneous confidence bands (SCB)
#' using a non-parametric bootstrap method for the mean outcome of regression
#' on a fixed test set design matrix
#'
#' @param df_fit A data frame containing the training design matrix used
#' to fit the linear model. Acceptable input format includes numeric and factor.
#' @param model A character string representing the formula for the linear model
#' (e.g., \code{"y ~ x1 + x2"}).
#' @param grid_df A data frame specifying the covariate settings that define the
#' mean outcome for which simultaneous confidence bands (SCB) are constructed.
#' Each row represents one covariate combination at which predictions and
#' SCBs are evaluated. Column names should match variables in the fitted model,
#' but `grid_df` may include only the subset of covariates of interest for the SCB
#' (it is not required to cover all model variables).
#' Default is `NULL`, in which case the SCB is constructed over the fitted values
#' based on 'df_fit`.
#' @param n_boot Number of bootstrap samples used in the non-parametric bootstrap
#' procedure to generate the empirical distribution. Default is 1000.
#' @param alpha Significance level for the confidence band
#' (e.g., 0.05 for 95% confidence). Default is 0.05.
#' @param grid_df_boot An optional data frame specifying the input grid at which
#' predictions are evaluated during bootstrap resampling.
#' This allows SCBs to be constructed on a denser set of covariate values
#' if desired. If NULL, uses \code{grid_df}.
#' If `grid_df` is set to NULL, `grid_df_boot` will also be set to `NULL`.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{scb_low}{Lower bound of the simultaneous confidence band.}
#'   \item{Mean}{Predicted mean response from the fitted model.}
#'   \item{scb_up}{Upper bound of the simultaneous confidence band.}
#'   \item{...}{All columns from \code{grid_df}, representing the prediction grid.}
#' }
#'
#' @importFrom stats quantile as.formula lm predict
#' @export
#'
#' @examples
#' set.seed(262)
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' epsilon <- rnorm(100,0,sqrt(2))
#' y <- -1 + x1 + 0.5 * x1^2 - 1.1 * x1^3 - 0.5 * x2 + 0.8 * x2^2 - 1.1 * x2^3 + epsilon
#' df <- data.frame(x1 = x1, x2 = x2, y = y)
#' grid <- data.frame(x1 = seq(-1, 1, length.out = 100),
#'                    x2 = seq(-1, 1, length.out = 100))
#' model <- "y ~ x1 + I(x1^2) + I(x1^3) + x2 + I(x2^2) + I(x2^3)"
#' results <- SCB_linear_outcome(df_fit = df, model = model, grid_df = grid)
#'
SCB_linear_outcome = function(df_fit, model, grid_df = NULL, n_boot = 1000,
                              alpha = 0.05, grid_df_boot = NULL){

  if(is.null(df_fit)) stop("`df_fit` must be provided.")
  if(is.null(model)) stop("`model` must be provided.")
  #if(is.null(grid_df)) stop("`grid_df` must be provided.")

  if(!is.data.frame(df_fit)){
    df_fit <- tryCatch(
      as.data.frame(df_fit),
      error = function(e) stop("`df_fit` must be a data.frame or coercible to a data.frame.")
    )
  }
  if((nrow(df_fit) == 0 || ncol(df_fit) == 0)){
    stop("`df_fit` is empty.")
  }

  # transform character variable to factor
  char_vars <- names(df_fit)[sapply(df_fit, is.character)]
  for (v in char_vars) {
    df_fit[[v]] <- factor(df_fit[[v]])
  }

  if(!is.null(grid_df)){
    if(!is.data.frame(grid_df)){
      grid_df <- tryCatch(
        as.data.frame(grid_df),
        error = function(e) stop("`grid_df` must be a data.frame or coercible to a data.frame.")
      )
    }
    if((nrow(grid_df) == 0 || ncol(grid_df) == 0)){
      grid_df <- NULL
    }
  }

  if(!inherits(model, "formula") && !is.character(model)) stop("`model` must be a formula or string.")
  if(!is.numeric(n_boot) || n_boot <= 0 || n_boot %% 1 != 0) stop("`n_boot` must be a positive integer.")
  if(!is.numeric(alpha) || alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).")

  formula_ <- as.formula(model)
  model_vars <- all.vars(formula_)
  if(!("." %in% model_vars)){
    # make sure that no vars in model_vars are missing in df_fit
    if (length(setdiff(model_vars, names(df_fit))) > 0) {
      stop(paste("`df_fit` is missing variables:", paste(setdiff(model_vars, names(df_fit)), collapse = ", ")))
    }
    if(length(model_vars) > 1){
    # for all model_vars, check if they're included in grid_df
    # make sure that the included vars match the format in df_fit
      if(!is.null(grid_df)){
        grid_df <- check_and_align_vars(df_fit, grid_df, model_vars[-1])
    # fill missing variables in grid_df that are included in model_vars
    # if input of grid_df is NULL, fill in all missing values and set them as 0/reference
        grid_df <- fill_missing_with_reference(df_fit, grid_df, model_vars[-1])
      }# else: grid_df is NULL, construct SCB over the fitted values based on 'df_fit'.
    } # else: y~1, no need to check grid_df, but intercept must be included
  }else{
    # y~.: check if y exists in df_fit
    if (length(setdiff(model_vars[1], names(df_fit))) > 0) {
      stop(paste("`df_fit` is missing variables:", model_vars[1]))
    }
    # for all vars in df_fit, check if they're included in grid_df
    # make sure that the included vars match the format in df_fit
    df_fit_no_outcome <- subset(df_fit, select = setdiff(names(df_fit), model_vars[1])) # remove outcome
    if (ncol(df_fit_no_outcome) != 0) {
      if(!is.null(grid_df)){
        grid_df <- check_and_align_vars(df_fit_no_outcome, grid_df)
    # fill missing variables in grid_df that are included in df_fit
    # if input of grid_df is NULL, fill in all missing values and set them as 0/reference
        grid_df <- fill_missing_with_reference(df_fit_no_outcome, grid_df)
      }# else: grid_df is NULL, construct SCB over the fitted values based on 'df_fit'.
    }#else: df_fit only has y, no need to check grid_df
  }

  if(!is.null(grid_df_boot)){
    if(!is.data.frame(grid_df_boot)){
      grid_df_boot <- tryCatch(
        as.data.frame(grid_df_boot),
        error = function(e) stop("`grid_df_boot` must be a data.frame or coercible to a data.frame.")
      )
    }
    if((nrow(grid_df_boot) == 0 || ncol(grid_df_boot) == 0)){
      grid_df_boot <- NULL
    }

    if(!("." %in% model_vars)){
      if(length(model_vars) > 1){
      # for all model_vars, check if they're included in grid_df_boot
      # make sure that the included vars match the format in grid_df
      # grid_df_boot is only for calculate the thres
      # for muneric variable, can be denser
      # for factor/character variable, must be included in grid_df
        if(!is.null(grid_df)){
          grid_df_boot <- check_and_align_vars(grid_df, grid_df_boot, model_vars[-1], grid_boot = TRUE)
      # fill missing variables in grid_df_boot that are included in model_vars
          grid_df_boot <- fill_missing_with_reference(grid_df, grid_df_boot, model_vars[-1])
        }else{
          grid_df_boot <- NULL
        }
      }#else: y~1, no need to check grid_df_boot, but intercept must be included
    }else{
      # y~.
      # for all vars in grid_df, check if they're included in grid_df_boot
      # make sure that the included vars match the format in grid_df
      if (ncol(df_fit_no_outcome) != 0){
        if(!is.null(grid_df)){
          grid_df_boot <- check_and_align_vars(grid_df, grid_df_boot, grid_boot = TRUE)
      # fill missing variables in grid_df_boot that are included in grid_df
          grid_df_boot <- fill_missing_with_reference(grid_df, grid_df_boot)
        }else{
          grid_df_boot <- NULL
        }
      }# else: only y included in the model, no need to check grid_df_boot
    }
  }

  fit <- lm(model, df_fit)
  y_hat <- predict(fit, grid_df, se.fit = TRUE, level = 1 - alpha) # for constructing the whole
  res_max_v <- rep(0,n_boot)
  if(is.null(grid_df_boot)){
    grid_df_boot <- grid_df
    y_hat_level_grid <- y_hat
  }else{
    y_hat_level_grid <- predict(fit, grid_df_boot, se.fit = TRUE, level = 1 - alpha)# bootstrap true mean
  }
  for(i in 1:n_boot){
    df_boot <- df_fit[sample(1:dim(df_fit)[1], replace = T),]
    if (!is.data.frame(df_boot)) {
      df_boot <- as.data.frame(df_boot)
      names(df_boot) <- names(df_fit)
    }
    fit_boot <- lm(model, df_boot)
    y_hat_boot <- predict(fit_boot, grid_df_boot, level = 1 - alpha, se.fit = TRUE	)
    residual <- abs(y_hat_boot$fit - y_hat_level_grid$fit)/y_hat_boot$se.fit
    res_max_v[i] <- max(residual)
  }
  thres <- quantile(res_max_v, probs = 1 - alpha)
  if(!is.null(grid_df)){
    sim_CB <- data.frame(scb_low = y_hat$fit - thres*y_hat$se.fit,
                         Mean = y_hat$fit,
                         scb_up = y_hat$fit + thres*y_hat$se.fit, grid_df)
  }else{
    sim_CB <- data.frame(scb_low = y_hat$fit - thres*y_hat$se.fit,
                         Mean = y_hat$fit,
                         scb_up = y_hat$fit + thres*y_hat$se.fit, df_fit)
  }
  return(sim_CB)
}

#' Expit (Inverse Logit) Function
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

#' Construct Simultaneous Confidence Bands for a Logistic Regression Outcome
#'
#' This function fits a logistic regression model and constructs simultaneous confidence bands (SCB)
#' using a non-parametric bootstrap method for the mean outcome of regression on a fixed test set design matrix
#'
#' @param df_fit A data frame containing the training design matrix used to fit the logistic model.
#' @param model A character string representing the formula for the logistic model (e.g., \code{"y ~ x1 + x2"}).
#' @param grid_df A data frame specifying the covariate settings that define the
#' mean outcome for which simultaneous confidence bands (SCB) are constructed.
#' Each row represents one covariate combination at which predictions and
#' SCBs are evaluated. Column names should match variables in the fitted model,
#' but `grid_df` may include only the subset of covariates of interest for the SCB
#' (it is not required to cover all model variables).
#' Default is `NULL`, in which case the SCB is constructed over the fitted values
#' based on 'df_fit`.
#' @param n_boot Number of bootstrap samples used in the non-parametric bootstrap
#' procedure to generate the empirical distribution. Default is 1000.
#' @param alpha Significance level for the confidence band (e.g., 0.05 for 95% confidence).
#' Default is 0.05.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{scb_low}{Lower bound of the simultaneous confidence band.}
#'   \item{Mean}{Predicted mean response from the fitted model.}
#'   \item{scb_up}{Upper bound of the simultaneous confidence band.}
#'   \item{...}{All columns from \code{grid_df}, representing the prediction grid.}
#' }
#'
#' @importFrom stats quantile as.formula glm predict
#' @export
#'
#' @examples
#' set.seed(262)
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' mu <- -1 + x1 + 0.5 * x1^2 - 1.1 * x1^3 - 0.5 * x2 + 0.8 * x2^2 - 1.1 * x2^3
#' p <- expit(mu)
#' y <- rbinom(100, size = 1, prob = p)
#' df <- data.frame(x1 = x1, x2 = x2, y = y)
#' grid <- data.frame(x1 = seq(-1, 1, length.out = 100), x2 = seq(-1, 1, length.out = 100))
#' model <- "y ~ x1 + I(x1^2) + I(x1^3) + x2 + I(x2^2) + I(x2^3)"
#' results <- SCB_logistic_outcome(df_fit = df, model = model, grid_df = grid)
#'
SCB_logistic_outcome = function(df_fit, model, grid_df = NULL, n_boot = 1000, alpha = 0.05){

  if(is.null(df_fit)) stop("`df_fit` must be provided.")
  if(is.null(model)) stop("`model` must be provided.")
  #if(is.null(grid_df)) stop("`grid_df` must be provided.")

  if(!is.data.frame(df_fit)){
    df_fit <- tryCatch(
      as.data.frame(df_fit),
      error = function(e) stop("`df_fit` must be a data.frame or coercible to a data.frame.")
    )
  }
  if((nrow(df_fit) == 0 || ncol(df_fit) == 0)){
    stop("`df_fit` is empty.")
  }

  # transform character variable to factor
  char_vars <- names(df_fit)[sapply(df_fit, is.character)]
  for (v in char_vars) {
    df_fit[[v]] <- factor(df_fit[[v]])
  }

  if(!is.data.frame(grid_df)){
    grid_df <- tryCatch(
      as.data.frame(grid_df),
      error = function(e) stop("`grid_df` must be a data.frame or coercible to a data.frame.")
    )
    if((nrow(grid_df) == 0 || ncol(grid_df) == 0)){
      grid_df <- NULL
    }
  }

  if(!inherits(model, "formula") && !is.character(model)) stop("`model` must be a formula or string.")
  if(!is.numeric(n_boot) || n_boot <= 0 || n_boot %% 1 != 0) stop("`n_boot` must be a positive integer.")
  if(!is.numeric(alpha) || alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).")

  formula_ <- as.formula(model)
  model_vars <- all.vars(formula_)
  if(!("." %in% model_vars)){
    # make sure that no vars in model_vars are missing in df_fit
    if (length(setdiff(model_vars, names(df_fit))) > 0) {
      stop(paste("`df_fit` is missing variables:", paste(setdiff(model_vars, names(df_fit)), collapse = ", ")))
    }
    if(length(model_vars) > 1){
      # for all model_vars, check if they're included in grid_df
      # make sure that the included vars match the format in df_fit
      if(!is.null(grid_df)){
        grid_df <- check_and_align_vars(df_fit, grid_df, model_vars[-1])
      # fill missing variables in grid_df that are included in model_vars
      # if input of grid_df is NULL, fill in all missing values and set them as 0/reference
        grid_df <- fill_missing_with_reference(df_fit, grid_df, model_vars[-1])
      }# else: grid_df is NULL, construct SCB over the fitted values based on 'df_fit'.
    } #else: y~1, no need to check grid_df, but intercept must be included
  }else{
    # y~.: check if y exists in df_fit
    if (length(setdiff(model_vars[1], names(df_fit))) > 0) {
      stop(paste("`df_fit` is missing variables:", model_vars[1]))
    }
    # for all vars in df_fit, check if they're included in grid_df
    # make sure that the included vars match the format in df_fit
    df_fit_no_outcome <- subset(df_fit, select = setdiff(names(df_fit), model_vars[1])) # remove outcome
    if (ncol(df_fit_no_outcome) != 0) {
      if(!is.null(grid_df)){
        grid_df <- check_and_align_vars(df_fit_no_outcome, grid_df)
      # fill missing variables in grid_df that are included in df_fit
      # if input of grid_df is NULL, fill in all missing values and set them as 0/reference
        grid_df <- fill_missing_with_reference(df_fit_no_outcome, grid_df)
      }# else: grid_df is NULL, construct SCB over the fitted values based on 'df_fit'.
    }#else: df_fit only has y, no need to check grid_df
  }

  fit <- suppressWarnings(glm(model, family = binomial(), data = df_fit)) # Suppress warning forfitted probabilities numerically 0 or 1
  y_hat <- predict(fit, grid_df, se.fit = TRUE, level = 1 - alpha) # bootstrap true mean
  res_max_v <- rep(0,n_boot)
  for(i in 1:n_boot){
    df_boot <- df_fit[sample(1:dim(df_fit)[1], replace = T),]
    if (!is.data.frame(df_boot)) {
      df_boot <- as.data.frame(df_boot)
      names(df_boot) <- names(df_fit)
    }
    fit_boot <- suppressWarnings(glm(model, family = binomial(), data = df_boot))
    y_hat_boot <- predict(fit_boot, grid_df, level = 1 - alpha, se.fit = TRUE	)
    residual <- abs(y_hat_boot$fit - y_hat$fit)/y_hat_boot$se.fit
    res_max_v[i] <- max(residual)
  }
  thres <- quantile(res_max_v, probs = 1 - alpha)
  if(!is.null(grid_df)){
    sim_CB <- data.frame(scb_low = expit(y_hat$fit - thres*y_hat$se.fit),
                         Mean = expit(y_hat$fit),
                        scb_up = expit(y_hat$fit + thres*y_hat$se.fit), grid_df)
  }else{
    sim_CB <- data.frame(scb_low = expit(y_hat$fit - thres*y_hat$se.fit),
                         Mean = expit(y_hat$fit),
                         scb_up = expit(y_hat$fit + thres*y_hat$se.fit), df_fit)
  }

  return(sim_CB)
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
#'   \item{scb_low}{Lower bound of the simultaneous confidence band.
#'   The first row corresponds to the intercept, and subsequent rows correspond to regression coefficients.}
#'   \item{Mean}{Estimated values. The first element is the intercept estimate, and the remaining are coefficient estimates.}
#'   \item{scb_up}{Upper bound of the simultaneous confidence band.
#'   The first row corresponds to the intercept, and subsequent rows correspond to regression coefficients.}
#' }
#'
#' @importFrom stats quantile as.formula lm glm predict binomial
#' @export
#'
#' @examples
#' library(MASS)
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
#' results <- SCB_regression_coef(df, model)
#'
SCB_regression_coef = function(df_fit, model, n_boot = 5000, alpha = 0.05, type = "linear"){

  if(is.null(df_fit)) stop("`df_fit` must be provided.")
  if(is.null(model)) stop("`model` must be provided.")

  if(!is.data.frame(df_fit)){
    df_fit <- tryCatch(
      as.data.frame(df_fit),
      error = function(e) stop("`df_fit` must be a data.frame or coercible to a data.frame.")
    )
  }
  if((nrow(df_fit) == 0 || ncol(df_fit) == 0)){
    stop("`df_fit` is empty.")
  }

  if(!inherits(model, "formula") && !is.character(model)) stop("`model` must be a formula or string.")
  if(!is.numeric(n_boot) || n_boot <= 0 || n_boot %% 1 != 0) stop("`n_boot` must be a positive integer.")
  if(!is.numeric(alpha) || alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).")

  formula_ <- as.formula(model)
  model_vars <- all.vars(formula_)
  if(!("." %in% model_vars)){
    if (length(setdiff(model_vars, names(df_fit))) > 0) {
      stop(paste("`df_fit` is missing variables:",
                 paste(setdiff(model_vars, names(df_fit)), collapse = ", ")))
    }
  }else{
    if (length(setdiff(model_vars[1], names(df_fit))) > 0) {
      stop(paste("`df_fit` is missing variables:", model_vars[1]))
    }
    #if(length(model_vars) > 0 && length(setdiff(names(df_fit), model_vars[1])) == 0){
      #stop("No covariates found in `df_fit`. Please add at least one covariate or adjust the input of `model`.")
    #}
  }

  # type: "linear" or "logistic"
  if(type == "linear"){
    fit <- lm(model, df_fit)
  }else if(type == "logistic"){
    fit <- suppressWarnings(glm(model, family = binomial(), data = df_fit)) # Suppress warning for 0 or 1 probability for simulation
  }else{
    stop("`type` must be either 'linear' or 'logistic'.")
  }

  sum_out <- summary(fit)
  coef_hat <- fit$coefficients # bootstrap true mean
  coef_sd <- sum_out$coefficients[,2]
  res_max_v <- rep(0,n_boot)
  for(i in 1:n_boot){
    df_boot <- df_fit[sample(1:dim(df_fit)[1], replace = T),, drop=F]
    if(type == "linear"){
      fit_boot <- lm(model, df_boot)
    }else if(type == "logistic"){
      fit_boot <- suppressWarnings(glm(model, family = binomial(), data = df_boot))
    }
    sum_out_boot <- summary(fit_boot)
    coef_hat_boot <- fit_boot$coefficients
    coef_sd_boot <- sum_out_boot$coefficients[,2]
    residual <- abs(coef_hat_boot - coef_hat)/coef_sd_boot
    res_max_v[i] <- max(residual)
  }
  thres <- quantile(res_max_v, probs = 1 - alpha)
  sim_CB <- data.frame(scb_low = coef_hat - thres*coef_sd, Mean = coef_hat, scb_up = coef_hat + thres*coef_sd)
  return(sim_CB)
}

#' Validate and align types/levels between training data and prediction grid
#'
#' Ensures that variables in `grid_df` are type-compatible with `df_fit` and that
#' factor (including ordered factor) levels are aligned to those used during model
#' fitting. Character columns in `df_fit` are first coerced to factors. For any
#' factor/ordered variable, `grid_df` is coerced to the same type and levels; any
#' unseen levels in `grid_df` will trigger an error.
#'
#' @param df_fit A data frame used for model fitting. Character columns will be
#'   coerced to factors before alignment.
#' @param grid_df A data frame of covariate settings (newdata) at which
#'   predictions/SCBs are to be evaluated.
#' @param model_vars A vector contained all the interested columns that appear
#' in `df_fit`. Only those variables are aligned in `grid_df`;
#' other columns are left unchanged. Default is NULL.
#' @param grid_boot A logic value for whether the function is for constructing
#' `grid_df_boot` or not.
#' @return \item{grid_df}{The prediction grid with variables aligned to `df_fit`.}
#'
#' @keywords internal
#'
#' @examples
#' # Used internally by SCB_linear_outcome, SCB_logistic_outcome
#'
check_and_align_vars <- function(df_fit, grid_df, model_vars = NULL, grid_boot = FALSE) {

  # find common vars
  if(is.null(model_vars)){
    common_vars <- intersect(names(grid_df), names(df_fit))
  }else{
    common_vars <- intersect(names(grid_df), model_vars)
  }

  for (nm in common_vars) {
    x_fit  <- df_fit[[nm]]
    x_new  <- grid_df[[nm]]

    # ——(A) factor/ordered factor, check levels
    if (is.factor(x_fit)) {
      lv <- levels(x_fit)
      is_ord <- is.ordered(x_fit)

      # if character var, transform it to factor, same levels  as df_fit
      if (is.character(x_new)) {
        unseen <- setdiff(unique(x_new), lv)
        if (length(unseen) > 0) {
          stop(sprintf("Column `%s` in grid_df has unseen level(s): %s",
                       nm, paste(unseen, collapse = ", ")))
        }
        present1 <- unique(as.character(x_fit[!is.na(x_fit)]))
        present2 <- unique(x_new[!is.na(x_new)])
        missing_in_new <- setdiff(present2, present1)
        if (length(missing_in_new) > 0) {
          stop(sprintf(
            "Check input for `%s`: contains extra value(s) not present in training/testing data -> %s",
            var, paste(missing_in_new, collapse = ", ")
          ))
        }
        grid_df[[nm]] <- factor(x_new, levels = lv, ordered = is_ord)
      } else if (is.factor(x_new)) {
        # if factor var, all levels should be included in df_fit
        unseen <- setdiff(unique(as.character(x_new)), lv)
        if (length(unseen) > 0) {
          stop(sprintf("Column `%s` in grid_df has unseen level(s): %s",
                       nm, paste(unseen, collapse = ", ")))
        }
        present1 <- unique(as.character(x_fit[!is.na(x_fit)]))
        present2 <- unique(as.character(x_new[!is.na(x_new)]))
        missing_in_new <- setdiff(present2, present1)
        if(grid_boot == FALSE){
          if (length(missing_in_new) > 0) {
            stop(sprintf(
              "Check input for `%s`: contains extra value(s) not present in training/testing data -> %s",
              var, paste(missing_in_new, collapse = ", ")
            ))
          }
        }else{
          missing_in_fit <- setdiff(present1, present2)
          if (length(missing_in_new) > 0 || length(missing_in_fit) > 0) {
            stop(sprintf(
              "Check input for `%s`: contains extra value(s) not present in training/testing data -> %s",
              var, paste(missing_in_new, collapse = ", ")
            ))
          }
        }
        grid_df[[nm]] <- factor(as.character(x_new), levels = lv, ordered = is_ord)
      } else {
        stop(sprintf("Column `%s` must be factor%sin grid_df (got %s).",
                     nm, if (is_ord) " (ordered) " else " ",
                     paste(class(x_new), collapse = "/")))
      }
      next
    }

    # ——(B) non-factor var
    cls_fit <- class(x_fit)
    cls_new <- class(x_new)

    same_main_class <- identical(cls_fit, cls_new) ||
      (is.numeric(x_fit) && (is.numeric(x_new) || is.integer(x_new))) ||
      (is.integer(x_fit) && is.numeric(x_new))

    if (!same_main_class) {
      stop(sprintf("Column `%s` in grid_df must have the same (or compatible) type as in df_fit. df_fit: %s, grid_df: %s",
                   nm, paste(cls_fit, collapse = "/"), paste(cls_new, collapse = "/")))
    }

    #if (is.numeric(x_fit) && is.integer(x_new)) {
      #grid_df[[nm]] <- as.numeric(x_new)
    #}
  }

  return(grid_df)
}

#' Fill missing variables in grid_df with reference values from df_fit
#'
#' @param df_fit A data frame used for model fitting. Character columns will be
#'   coerced to factors before alignment.
#' @param grid_df A data frame of covariate settings (newdata) at which
#'   predictions/SCBs are to be evaluated.
#' @param model_vars A vector contained all the interested columns that appear
#' in `df_fit`. Only those variables are aligned in `grid_df`;
#' other columns are left unchanged. Default is NULL.
#'
#' @return grid_df with missing columns filled:
#'   - factor/ordered: reference = first level in df_fit
#'   - character: coerced to factor, reference = first unique value
#'   - numeric/integer: reference = numeric_ref (default 0)
#'
#' @keywords internal
#'
#' @examples
#' # Used internally by SCB_linear_outcome, SCB_logistic_outcome
#'
fill_missing_with_reference <- function(df_fit, grid_df, model_vars = NULL) {
  if(is.null(model_vars)){
    missing_vars <- setdiff(names(df_fit), names(grid_df)) # in df_fit but not in grid_df
  }else{
    missing_vars <- setdiff(model_vars, names(grid_df))
  }

  if (length(missing_vars)) {
    nr <- nrow(grid_df)
    for (var in missing_vars) {
      x_fit <- df_fit[[var]]
      if (is.factor(x_fit)) {
        lv <- levels(x_fit)
        ref <- lv[1]
        grid_df[[var]] <- factor(rep(ref, nr),
                                 levels = lv,
                                 ordered = is.ordered(x_fit))
      } else if (is.integer(x_fit)) {
        grid_df[[var]] <- rep(as.integer(0), nr)
      } else if (is.numeric(x_fit)) {
        grid_df[[var]] <- rep(as.numeric(0), nr)
      } else {
        stop(sprintf("Variable `%s`: unsupported type %s",
                     var, paste(class(x_fit), collapse = "/")))
      }
    }
  }
  grid_df
}
