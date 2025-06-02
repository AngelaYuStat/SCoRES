#' Construct CMA Confidence Intervals via Parametric Method
#'
#' This function computes Correlation and Multiplicity Adjusted (CMA) confidence bands for a specified group in a functional outcome regression model
#' using parameter simulations approach with Gaussian multiplier bootstrap.
#'
#' @param object A fitted functional outcome regression model object (e.g., from mgcv::bam()/mgcv::gam()).
#' @param model A character string specifying the model type used (e.g., gam or bam).
#' @param time A character string specifying the name of the time variable used in the model.
#' @param range A numeric vector specifying the range or grid of time points interested. Defauts is the whole time points from the fitted model.
#' @param groups A character vector specifying the names of group-related covariates of interest. If NULL, CIs are for the reference group.
#' @param subject A character string specifying the name of the subject-level random effect variable.
#' @param nboot An integer specifying the number of bootstrap samples used to construct the confidence bands. Defaults to 10,000.
#'
#' @returns A list containing:
#'   \item{yhat}{Estimated mean function for the group of interest.}
#'   \item{time}{The time points used.}
#'   \item{se_hat}{Standard errors of the estimated means.}
#'   \item{scb_low}{Lower bound of the simultaneous confidence band.}
#'   \item{scb_uo}{Upper bound of the simultaneous confidence band.}
#'   \item{type}{A character description of the output type.}
#'
#' @import stats
#' @export
#'
#' @examples
#' # example using ccds data
#' data(ccds_fpca)
#'
#' fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
#'   s(seconds, by = use, k=30, bs = "cr") +
#'   s(subject, by = Phi1, bs="re") +
#'   s(subject, by = Phi2, bs="re")+
#'   s(subject, by = Phi3, bs="re") +
#'   s(subject, by = Phi4, bs="re"),
#'   method = "fREML", data = ccds_fpca, discrete = TRUE)

#' cma(fosr_mod, "gam", "seconds", groups = "use", subject = "subject")
#'
cma = function(object, model, time, range = NULL, groups = NULL, subject = NULL, nboot = NULL){

  mod_coef <- object$coefficients

  if(model %in% c("gam", "bam")){
    mod_cov <- vcov(object) # containing all variance-covariance info for object
  }

  # get all names for terms
  coef_names <- names(mod_coef)

  # get index of the initial terms
  idx <- grep("^\\(Intercept\\)$|^s\\(seconds\\)\\.[0-9]+$", coef_names)

  # initialize dataframe
  predictors <- all.vars(formula(object))[-1]

  df_pred <- as_tibble(
    setNames(
      lapply(predictors, function(x) 0),
      predictors
    )
  )

  if(!is.null(groups)){
    # for group interested
    groups_idx <- unlist(lapply(groups, function(g) grep(g, coef_names)))
    idx <- sort(unique(c(idx, groups_idx)))
    df_pred[groups] <- 1
  }

  if(!is.null(subject)){
    df_pred[subject] <- as.character(model.frame(object)[[subject]][1])
  }

  if(is.null(range)){
    s_pred <- unique(model.frame(object)[[time]])
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

  # Number of bootstrap samples (B)
  if(is.null(nboot)){
    nboot <- 1e4
  }

  # Set up container for bootstrap
  yhat_boot <- matrix(NA, nboot, length(s_pred))

  for (i in 1:nboot) {
    beta_boot_i <- mvrnorm(n = 1, mu = mod_coef, Sigma = mod_cov)
    yhat_boot[i, ] <- lpmat %*% beta_boot_i
  }

  # Find the max statistic
  dvec <- apply(yhat_boot, 1, function(x) max(abs(x - pred_df$mean) / pred_df$se)) # Substract mean estimate and divided by Df (element by element)

  # Get 95% global confidence band
  Z_global <- quantile(dvec, 0.95)
  y_hat_LB_global <- pred_df$mean - Z_global * pred_df$se
  y_hat_UB_global <- pred_df$mean + Z_global * pred_df$se

  global_df <- list(
    yhat = pred_df$mean,
    time = s_pred,
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
#' data(ccds_cma)
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


