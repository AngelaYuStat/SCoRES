#' Functional Outcome Regression Prediction with Group-Specific Inference
#'
#' This function is an internal function for constructing SCBs for functional data.
#'
#' @param data_df A functional data frame.
#' The input data must have column names, and should contain the functional outcome, time and subject
#' @param object A fitted Function-on-Scalar Regression (FoSR) model object (e.g., from mgcv::bam()/mgcv::gam()).
#' @param fitted Logical. Whether to estimate the simultaneous confidence bands for fitted mean function or fitted parameter function
#'   \itemize{
#'     \item \code{TRUE} - Estimate the simultaneous confidence bands for fitted mean outcome function.
#'     \item \code{FALSE} - estimate the simultaneous confidence bands for fitted parameter function.
#'     }
#'   Default is \code{TRUE}.
#' @param time A character string specifying the name of the time variable used in the model.
#' @param range A numeric vector specifying the range or grid of time points interested. Defaut is the whole time points from the fitted model.
#' @param group_name A character vector that specifies the names of grouping scalar variables to analyze.
#'   Default is \code{NULL}, representing the reference group.
#' @param group_value A numeric vector that specifies the corresponding values of the variables listed in \code{group_name}.
#'   Each entry in \code{group_value} should match the respective variable in \code{group_name}.
#'   Character/factor is not allowed.
#'   Default is \code{NULL}, representing the reference group.
#' @param subject A optional character string specifying the name of the subject-level random effect variable, if included in model fitting.
#'
#' @returns A list containing the following elements:
#' \describe{
#'   \item{s_pred}{Numeric vector of sorted unique time points used for prediction}
#'   \item{pred_df}{Data frame with prediction results, containing:
#'     \itemize{
#'       \item \code{mean}: Predicted mean values
#'       \item \code{se}: Standard errors
#'     }
#'   }
#'   \item{lpmat}{Linear predictor matrix (design matrix) used for confidence interval calculations}
#'   \item{mod_coef}{Vector of model coefficients for selected group}
#'   \item{mod_cov}{Variance-covariance matrix corresponding to the selected group coefficients}
#' }
#'
#' @importFrom stats formula model.frame vcov
#' @importFrom dplyr %>% mutate
#' @importFrom tibble as_tibble
#' @importFrom stats predict
#'
#' @examples
#' library(mgcv)
#' data(pupil)
#' pupil_fpca <- prepare_pupil_fpca(pupil)
#'
#' fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
#'   s(seconds, by = use, k=30, bs = "cr") +
#'   s(id, by = Phi1, bs="re") +
#'   s(id, by = Phi2, bs="re")+
#'   s(id, by = Phi3, bs="re") +
#'   s(id, by = Phi4, bs="re"),
#'   method = "fREML", data = pupil_fpca, discrete = TRUE)
#'
#' results <- mean_response_predict(pupil_fpca, fosr_mod, fitted = TRUE, time = "seconds",
#' group_name = "use", group_value = 1, subject = "id")
#'
#' results <- mean_response_predict(pupil_fpca, fosr_mod, fitted = FALSE, time = "seconds",
#' group_name = "use", group_value = 1, subject = "id")
#'
#' @export
mean_response_predict = function(data_df, object, fitted = TRUE, time, range = NULL, group_name, group_value, subject = NULL){

  mod_coef <- object$coefficients
  mod_cov <- vcov(object) # containing all variance-covariance info for object

  coef_names <- names(mod_coef)
  # get index of the initial terms
  pattern <- paste0("^\\(Intercept\\)$|^s\\(", time, "\\)\\.[0-9]+$")
  idx <- grep(pattern, coef_names)
  intercept_idx <- idx

  # initialize dataframe
  predictors <- all.vars(formula(object))[-1]

  df_pred <- as_tibble(
    setNames(
      lapply(predictors, function(x) 0),
      predictors
    )
  )

  groups_idx <- NULL
  if (!is.null(group_name) && !is.null(group_value)){
    for (i in seq_along(group_name)) {
      var <- group_name[i]
      val <- group_value[i]

      if (!(var %in% names(data_df))) {
        stop(paste0("Variable '", var, "' not found in `data_df`."))
      }

      if (is.character(data_df[[var]])) {
        stop(paste0("The variable '", var, "' in `data_df` is of type character. ",
                    "Please convert it to a numeric and consider refit your functional object."))
      }

      if (is.factor(data_df[[var]])) {
        stop(paste0("The variable '", var, "' in `data_df` is of type factor. ",
                    "Please convert it to a numeric and consider refit your functional object."))
      } else {  # numeric
        unique_vals <- unique(data_df[[var]])
        if (!(val %in% unique_vals)) {
          stop(paste0("Value '", val, "' not found in numeric variable '", var, "'."))
        }
        df_pred[[var]] <- val
      }
    }
    # for group interested
    groups_idx <- unlist(lapply(group_name, function(g) grep(g, coef_names)))
    if (fitted){
      idx <- sort(unique(c(idx, groups_idx)))
    }else{
      idx <- sort(unique(c(groups_idx)))
    }
  }

  if(!is.null(subject)){
    df_pred[subject] <- as.character(model.frame(object)[[subject]][1])
  }

  if(is.null(range)){
    s_pred <- sort(unique(model.frame(object)[[time]]))
  }else{
    s_pred <- sort(unique(range))
  }

  df_pred <- df_pred[rep(1, length(s_pred)), ]
  df_pred[[time]] <- s_pred

  # prepare design matrix
  lpmat <- predict(object, newdata=df_pred, se.fit=TRUE, type = "lpmatrix")

  # get mean response by groups with standard errors
  if (!fitted && !is.null(groups_idx)){
    # fitted = FALSE, fit mean outcome for linear combination
    # of parameter corresponding to the group specified without intercept
    lpmat_no_intercept <- lpmat[, -intercept_idx, drop = FALSE]
    mod_coef_no_intercept <- mod_coef[-intercept_idx]
    mod_cov_no_intercept <- mod_cov[-intercept_idx, -intercept_idx]
    pred_df <- df_pred %>%
      mutate(mean = c(lpmat_no_intercept %*% mod_coef_no_intercept),
             se = c(sqrt(diag(lpmat_no_intercept %*% mod_cov_no_intercept %*% t(lpmat_no_intercept)))))
  }else{
    # fitted = FALSE, fit mean outcome for the group specified
    # fitted = TRUE with no group specified, fit intercept
    pred_df <- df_pred %>% mutate(mean = c(lpmat %*% mod_coef), se = c(sqrt(diag(lpmat %*% mod_cov %*% t(lpmat)))))
  }
  lpmat <- lpmat[, idx]
  # extract means and variance for group interested
  mod_coef <- mod_coef[idx]
  mod_cov <- mod_cov[idx, idx]

  return(list(s_pred = s_pred, pred_df = pred_df, lpmat = lpmat, mod_coef = mod_coef, mod_cov = mod_cov))
}
#' Construct CMA Confidence Intervals via Parametric Method
#'
#' This function computes Correlation and Multiplicity Adjusted (CMA) confidence bands for a specified group in a functional outcome regression model
#' using parameter simulations approach with Gaussian multiplier bootstrap.
#'
#' @param data_df A functional data frame.
#' The input data must have column names, and should contain the functional outcome, time and subject
#' @param object A fitted Function-on-Scalar Regression (FoSR) object (e.g., from mgcv::bam()/mgcv::gam()).
#' @param fitted Logical. Whether to estimate the simultaneous confidence bands for fitted mean function or fitted parameter function
#'   \itemize{
#'     \item \code{TRUE} - Estimate the simultaneous confidence bands for fitted mean outcome function.
#'     \item \code{FALSE} - estimate the simultaneous confidence bands for fitted parameter function.
#'     }
#'   Default is \code{TRUE}.
#' @param alpha Significance level for SCB. Default is 0.05.
#' @param time A character string specifying the name of the time variable used in the model.
#' @param range A numeric vector specifying the range or grid of time points interested. Defaut is the whole time points from the fitted model.
#' @param group_name A character vector that specifies the names of grouping scalar variables to analyze.
#'   Default is \code{NULL}, representing the reference group.
#' @param group_value A numeric vector that specifies the corresponding values of the variables listed in \code{group_name}.
#'   Each entry in \code{group_value} should match the respective variable in \code{group_name}.
#'   Character/factor is not allowed.
#'   Default is \code{NULL}, representing the reference group.
#' @param subject A optional character string specifying the name of the subject-level random effect variable, if included in model fitting.
#' @param nboot An integer specifying the number of bootstrap samples used to construct the confidence bands. Default is 10,000.
#'
#' @returns A list containing:
#'   \item{yhat}{Estimated mean function for the group of interest.}
#'   \item{time}{The time points used.}
#'   \item{se_hat}{Standard errors of the estimated means.}
#'   \item{scb_low}{Lower bound of the simultaneous confidence band.}
#'   \item{scb_up}{Upper bound of the simultaneous confidence band.}
#'   \item{type}{A character description of the output type.}
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats quantile
#'
#' @export
#'
#' @examples
#' # example using pupil data
#' library(mgcv)
#' data(pupil)
#' pupil_fpca <- prepare_pupil_fpca(pupil)
#'
#' fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
#'   s(seconds, by = use, k=30, bs = "cr") +
#'   s(id, by = Phi1, bs="re") +
#'   s(id, by = Phi2, bs="re")+
#'   s(id, by = Phi3, bs="re") +
#'   s(id, by = Phi4, bs="re"),
#'   method = "fREML", data = pupil_fpca, discrete = TRUE)
#'
#' results <- cma(pupil_fpca, fosr_mod, fitted = TRUE, time = "seconds",
#' group_name = "use", group_value = 1, subject = "id")
#'
#' results <- cma(pupil_fpca, fosr_mod, fitted = FALSE, time = "seconds",
#' group_name = "use", group_value = 1, subject = "id")
#'
cma = function(data_df, object, fitted = TRUE, alpha = 0.05, time, range = NULL, group_name = NULL, group_value = NULL, subject = NULL, nboot = NULL){

  if (is.null(data_df)) {
    stop("`data_df` must be provided.")
  }

  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1){
    stop("`alpha` must be in (0, 1).")
  }

  # ---- Ensure 'data' is a named data frame ----
  if(!is.data.frame(data_df)){
    data_df <- tryCatch(
      as.data.frame(data_df),
      error = function(e) stop("`data_df` must be a data.frame or coercible to a data.frame.")
    )
  }

  # Check the existance of column names
  if (is.null(colnames(data_df))) {
    stop("`data_df` must have column names.")
  }

  if (is.null(object)) {
    stop("`object` must be provided.")
  }

  if (is.null(time)) {
    stop("`time` must be provided.")
  }
  if (!is.null(group_name) && !is.null(group_value)){
    if (length(group_name) != length(group_value)) {
      stop("The length of `group_name` and `group_value` must be the same.")
    }
  }else{
    if(!is.null(group_name) && is.null(group_value)|is.null(group_name) && !is.null(group_value)){
      stop("`group_name` and `group_value` must be provided.")
    }
  }
  results <- mean_response_predict(data_df, object, fitted, time, range, group_name, group_value, subject)

  # Number of bootstrap samples (B)
  if(is.null(nboot)){
    nboot <- 1e4
  }else{
    if (!is.numeric(nboot) || nboot <= 0 || nboot %% 1 != 0){
      stop("`nboot` must be a positive integer.")
    }
  }

  # Set up container for bootstrap
  yhat_boot <- matrix(NA, nboot, length(results$s_pred))

  for (i in 1:nboot) {
    beta_boot_i <- mvrnorm(n = 1, mu = results$mod_coef, Sigma = results$mod_cov)
    yhat_boot[i, ] <- results$lpmat %*% beta_boot_i
  }

  pred_df <- results$pred_df
  # Find the max statistic
  dvec <- apply(yhat_boot, 1, function(x) max(abs(x - pred_df$mean) / pred_df$se)) # Substract mean estimate and divided by Df (element by element)

  # Get 95% global confidence band
  Z_global <- quantile(dvec, 1 - alpha)
  y_hat_LB_global <- pred_df$mean - Z_global * pred_df$se
  y_hat_UB_global <- pred_df$mean + Z_global * pred_df$se

  return(list(
    yhat = pred_df$mean,
    time = results$s_pred,
    se_hat = pred_df$se,
    scb_low = y_hat_LB_global,
    scb_up = y_hat_UB_global,
    type = "Global Confidence Interval (CMA)"
  ))

}

#' Prepare Pupil FPCA Dataset
#'
#' Processes data by fitting a mean GAM model, extracting residuals, performing FPCA,
#' and merging the results to create an enhanced dataset for functional regression analysis.
#'
#' @param input_data Raw pupil data
#' @param k_mean Number of basis functions for mean model smooth terms (default: 30)
#' @param k_fpca Number of knots for FPCA estimation (default: 15)
#'
#' @return A tibble containing:
#'   \itemize{
#'     \item Original pupil variables
#'     \item FPCA eigenfunctions (Phi1, Phi2,...)
#'     \item Sorted by subject and time
#'   }
#'
#' @examples
#' library(mgcv)
#' data(pupil)
#' processed_data <- prepare_pupil_fpca(pupil)
#'
#' @importFrom dplyr mutate filter select arrange left_join %>%
#' @importFrom tidyr pivot_wider
#' @importFrom refund fpca.face
#' @importFrom tibble as_tibble
#'
#' @export
prepare_pupil_fpca <- function(input_data, k_mean = 30, k_fpca = 15) {

  # Fit mean model
  mean_mod <- mgcv::gam(
    percent_change ~ s(seconds, k = k_mean, bs = "cr") +
      s(seconds, by = use, k = k_mean, bs = "cr") +
      s(seconds, by = age, k = k_mean, bs = "cr") +
      s(seconds, by = gender, k = k_mean, bs = "cr"),
    data = input_data, method = "REML"
  )

  # Prepare residuals
  resid_df <- input_data %>%
    filter(!is.na(percent_change)) %>%
    select(id, seconds) %>%
    mutate(resid = mean_mod$residuals) %>%
    pivot_wider(
      names_from = seconds,
      values_from = resid,
      names_prefix = "resid."
    )

  resid_mat <- as.matrix(resid_df[, -1])
  rownames(resid_mat) <- resid_df$id

  # FPCA estimation
  fpca_results <- fpca.face(
    Y = resid_mat,
    argvals = unique(input_data$seconds),
    knots = k_fpca
  )

  # Create output dataset
  eigenfunctions <- as.data.frame(fpca_results$efunctions)
  colnames(eigenfunctions) <- paste0("Phi", seq(1, fpca_results$npc))
  eigenfunctions$seconds <- unique(input_data$seconds)

  output_data <- input_data %>%
    left_join(eigenfunctions, by = "seconds") %>%
    as_tibble() %>%
    arrange(id, seconds) %>%
    mutate(id = factor(id))

  return(output_data)
}

