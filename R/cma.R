#' Functional Outcome Regression Prediction with Group-Specific Inference
#'
#' This function is an internal function for constructing SCBs for functional data.
#'
#' @param data A functional data frame, matrix or tibble.
#' The input data must have column names, and should contain the functional outcome, time and subject
#' @param object A fitted Function-on-Scalar Regression (FoSR) model object (e.g., from mgcv::bam()/mgcv::gam()).
#' @param time A character string specifying the name of the time variable used in the model.
#' @param range A numeric vector specifying the range or grid of time points interested. Defaut is the whole time points from the fitted model.
#' @param group_name A character vector that specifies the names of grouping scalar variables to analyze.
#' @param group_value A vector (numeric or character) that specifies the corresponding values of the variables listed in \code{group_name}.
#'   Each entry in \code{group_value} should match the respective variable in \code{group_name}.
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
#' @examples
#' data(ccds)
#' ccds_fpca <- prepare_ccds_fpca(ccds)
#'
#' fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
#'   s(seconds, by = use, k=30, bs = "cr") +
#'   s(subject, by = Phi1, bs="re") +
#'   s(subject, by = Phi2, bs="re")+
#'   s(subject, by = Phi3, bs="re") +
#'   s(subject, by = Phi4, bs="re"),
#'   method = "fREML", data = ccds_fpca, discrete = TRUE)
#'
#' results <- mean_response_predict(ccds_fpca, fosr_mod, time = "seconds",
#' group_name = "use", group_value = 1, subject = "subject")
#'
#' @export
mean_response_predict = function(data, object, time, range = NULL, group_name, group_value, subject = NULL){

  mod_coef <- object$coefficients
  mod_cov <- vcov(object) # containing all variance-covariance info for object

  coef_names <- names(mod_coef)
  # get index of the initial terms
  pattern <- paste0("^\\(Intercept\\)$|^s\\(", time, "\\)\\.[0-9]+$")
  idx <- grep(pattern, coef_names)

  # initialize dataframe
  predictors <- all.vars(formula(object))[-1]

  df_pred <- as_tibble(
    setNames(
      lapply(predictors, function(x) 0),
      predictors
    )
  )

  if(!is.null(group_name)){
    # for group interested
    groups_idx <- unlist(lapply(group_name, function(g) grep(g, coef_names)))
    idx <- sort(unique(c(idx, groups_idx)))
    for (i in seq_along(group_name)) {
      var <- group_name[i]
      val <- group_value[i]

      if (is.character(df_pred[[var]])) {
        stop(paste0("The variable '", var, "' is of type character. ",
                    "Please convert it to a factor to ensure consistent prediction behavior."))
      } else if (is.factor(df_pred[[var]])) {
        df_pred[[var]] <- factor(val, levels = levels(df_pred[[var]]))
      } else {
        df_pred[[var]] <- val
      }
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
  pred_df <- df_pred %>% mutate(mean = c(lpmat %*% mod_coef), se = c(sqrt(diag(lpmat %*% mod_cov %*% t(lpmat)))))

  lpmat <- lpmat[, idx]
  # exract means and variance for group interested
  mod_coef <- mod_coef[idx]
  mod_cov <- mod_cov[idx, idx]

  return(list(s_pred = s_pred, pred_df = pred_df, lpmat = lpmat, mod_coef = mod_coef, mod_cov = mod_cov))
}
#' Construct CMA Confidence Intervals via Parametric Method
#'
#' This function computes Correlation and Multiplicity Adjusted (CMA) confidence bands for a specified group in a functional outcome regression model
#' using parameter simulations approach with Gaussian multiplier bootstrap.
#'
#' @param data A functional data frame, matrix or tibble.
#' The input data must have column names, and should contain the functional outcome, time and subject
#' @param object A fitted Function-on-Scalar Regression (FoSR) object (e.g., from mgcv::bam()/mgcv::gam()).
#' @param alpha Significance level for SCB. Default is 0.05.
#' @param time A character string specifying the name of the time variable used in the model.
#' @param range A numeric vector specifying the range or grid of time points interested. Defaut is the whole time points from the fitted model.
#' @param group_name A character vector that specifies the names of grouping scalar variables to analyze.
#' @param group_value A vector (numeric or character) that specifies the corresponding values of the variables listed in \code{group_name}.
#'   Each entry in \code{group_value} should match the respective variable in \code{group_name}.
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
#' @import stats
#' @export
#'
#' @examples
#' # example using ccds data
#' data(ccds)
#' ccds_fpca <- prepare_ccds_fpca(ccds)
#'
#' fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
#'   s(seconds, by = use, k=30, bs = "cr") +
#'   s(subject, by = Phi1, bs="re") +
#'   s(subject, by = Phi2, bs="re")+
#'   s(subject, by = Phi3, bs="re") +
#'   s(subject, by = Phi4, bs="re"),
#'   method = "fREML", data = ccds_fpca, discrete = TRUE)
#'
#' results <- cma(ccds_fpca, fosr_mod, time = "seconds", group_name = "use", group_value = 1, subject = "subject")
#'
cma = function(data, object, alpha = 0.05, time, range = NULL, group_name, group_value, subject = NULL, nboot = NULL){

  if (is.null(data)) {
    stop("Must provide the origin data!")
  }

  # ---- Ensure 'data' is a named data frame ----
  if (!is.data.frame(data)) {
    if (is.matrix(data)) {
      data <- as.data.frame(data)
    } else {
      stop("Input must be a named data frame or matrix.")
    }
  }

  # Check the existance of column names
  if (is.null(colnames(data))) {
    stop("Input data must have column names.")
  }

  if (is.null(object)) {
    stop("Must provide a fitted functional regression model object.")
  }
  if (is.null(time)) {
    stop("Must provide the time variable name.")
  }
  if (is.null(group_name)) {
    stop("Must provide a specified group variable name.")
  }
  if (is.null(group_value)) {
    stop("Must provide a specified group level.")
  }
  if (length(group_name) != length(group_value)) {
    stop("The length of 'group_name' and 'group_value' must be the same.")
  }
  results <- mean_response_predict(data, object, time, range, group_name, group_value, subject)

  # Number of bootstrap samples (B)
  if(is.null(nboot)){
    nboot <- 1e4
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

  global_df <- list(
    yhat = pred_df$mean,
    time = results$s_pred,
    se_hat = pred_df$se,
    scb_low = y_hat_LB_global,
    scb_up = y_hat_UB_global,
    type = "Global Confidence Interval (CMA)"
  )

  return(global_df)
}

#' Plot Simultaneous Confidence Bands from CMA Output
#'
#' This function visualizes the estimated mean function and simultaneous confidence bands
#' (SCB) from a CMA result object. It supports customization of colors, line types, axis labels,
#' and legend display.
#'
#' @param data A data frame or tibble containing the results of the \code{cma()} function.
#'             Must include columns: \code{time}, \code{yhat}, \code{scb_low}, and \code{scb_up}.
#' @param xlab A character string for the x-axis label. Default is \code{"Time"}.
#' @param ylab A character string for the y-axis label. Default is \code{"Estimated Mean"}.
#' @param line_color Color for the main estimated curve (\code{yhat}). Default is \code{"black"}.
#' @param ci_color Color for the confidence bands (\code{scb_low} and \code{scb_up}). Default is \code{"gray50"}.
#' @param ci_linetype Line type for the confidence bands. Default is \code{"dashed"}.
#' @param group_label A character label for the color legend. Default is \code{"Group"}.
#' @param title A character string for the plot title. Default is \code{"Simultaneous Confidence Band"}.
#' @param show_legend Logical. Whether to display the legend. Default is \code{FALSE}.
#'
#' @returns A \code{ggplot2} object representing the confidence band plot.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#'
#' # example using ccds data
#' data(ccds)
#' ccds_fpca <- prepare_ccds_fpca(ccds)
#'
#' fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
#'   s(seconds, by = use, k=30, bs = "cr") +
#'   s(subject, by = Phi1, bs="re") +
#'   s(subject, by = Phi2, bs="re")+
#'   s(subject, by = Phi3, bs="re") +
#'   s(subject, by = Phi4, bs="re"),
#'   method = "fREML", data = ccds_fpca, discrete = TRUE)
#'
#' ccds_cma <- cma(fosr_mod, time = "seconds", groups = "use", subject = "subject")
#' ccds_cma <- tibble::as_tibble(ccds_cma)
#' plot_cma(ccds_cma, xlab = "Time (s)", ylab = "Mean",
#'          line_color = "blue", ci_color = "gray", ci_linetype = "dotted",
#'          title = "Global SCB for Treatment Group")
#'
plot_cma <- function(data,
                     xlab = "Time",
                     ylab = "Estimated Mean",
                     line_color = "black",
                     ci_color = "gray50",
                     ci_linetype = "dashed",
                     group_label = "Group",
                     title = "Simultaneous Confidence Band",
                     show_legend = FALSE) {

  p <- ggplot(data, aes(x = time, color = type)) +
    geom_line(aes(y = yhat), size = 1) +
    geom_line(aes(y = scb_low), linetype = ci_linetype) +
    geom_line(aes(y = scb_up), linetype = ci_linetype) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
    theme_minimal() +
    labs(
      title = title,
      x = xlab,
      y = ylab,
      color = group_label
    )

  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }

  return(p)
}

#' Prepare ccds FPCA Dataset
#'
#' Processes ccds data by fitting a mean GAM model, extracting residuals, performing FPCA,
#' and merging the results to create an enhanced dataset for functional regression analysis.
#'
#' @param input_data Raw ccds data frame containing:
#'   \itemize{
#'     \item \code{percent_change}: Functional response variable
#'     \item \code{seconds}: Time variable
#'     \item \code{use}: Binary grouping variable
#'     \item \code{subject}: Subject identifier
#'   }
#' @param k_mean Number of basis functions for mean model smooth terms (default: 30)
#' @param k_fpca Number of knots for FPCA estimation (default: 15)
#'
#' @return A tibble containing:
#'   \itemize{
#'     \item Original ccds variables
#'     \item FPCA eigenfunctions (Phi1, Phi2,...)
#'     \item Sorted by subject and time
#'   }
#'
#' @examples
#' data(ccds)
#' processed_data <- prepare_ccds_fpca(ccds)
#'
#' @importFrom dplyr mutate filter select arrange left_join
#' @importFrom tidyr pivot_wider
#' @importFrom refund fpca.face
#' @importFrom tibble as_tibble
#'
#' @export
prepare_ccds_fpca <- function(input_data, k_mean = 30, k_fpca = 15) {

  # Fit mean model
  mean_mod <- gam(
    percent_change ~ s(seconds, k = k_mean, bs = "cr") +
      s(seconds, by = use, k = k_mean, bs = "cr"),
    data = input_data, method = "REML"
  )

  # Prepare residuals
  resid_df <- input_data %>%
    filter(!is.na(percent_change)) %>%
    select(subject, seconds) %>%
    mutate(resid = mean_mod$residuals) %>%
    pivot_wider(
      names_from = seconds,
      values_from = resid,
      names_prefix = "resid."
    )

  resid_mat <- as.matrix(resid_df[, -1])
  rownames(resid_mat) <- resid_df$subject

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
    arrange(subject, seconds) %>%
    mutate(subject = factor(subject))

  return(output_data)
}
