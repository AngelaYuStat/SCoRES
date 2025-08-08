#' Construct Simultaneous Confidence Bands (SCB) For One Dimensional Functional Data
#'
#' This function builds simultaneous confidence bands through two distinct approaches..
#'
#' @param data_df A functional data frame.
#' The input data must have column names, and should contain the functional outcome, time, subject
#' @param object A fitted Function-on-Scalar Regression (FoSR) object (e.g., from mgcv::bam()/mgcv::gam()). Default is \code{NULL}
#' @param method A character string specifying the approach:
#'   \itemize{
#'     \item \code{"cma"} - Correlation and Multiplicity Adjusted (CMA) confidence bands
#'           via parametric approach (requires a fitted functional regression model)
#'     \item \code{"multiplier"} - Dense confidence bands via Multiplier-t Bootstrap method
#'     For \code{method = "multiplier"} option, the user should check whether their input data (\code{data}) has zero entries on every subject for a specific time point.
#'     The function only allow all zero entries for time point zero. Otherwise, the function will return an error.
#'     Besides, for Multiplier Bootstrap, if there exist NA in data, the function will impute these NA's via fpca.face.
#'   }
#' @param fitted Logical. Whether to estimate the simultaneous confidence bands for fitted mean function or fitted parameter function
#'   \itemize{
#'     \item \code{TRUE} - Estimate the simultaneous confidence bands for fitted mean outcome function.
#'     \item \code{FALSE} - estimate the simultaneous confidence bands for fitted parameter function.
#'     }
#'   Default is \code{TRUE}.
#' @param est_mean Logical. Whether to use the fitted mean from a functional regression model in Multiplier Bootstrap (\code{method = "multiplier"}) for estimating the quantile.
#' This argument can be ignored if choose \code{method = "cma"}.
#'   If \code{TRUE}, the function will use:
#'   \itemize{
#'     \item \strong{Fitted mean:} The estimated mean function from a fitted model (e.g., \code{mgcv::bam}).
#'     Need to provide a fitted functional outcome regression model object. Otherwise, the sample mean will be calculated.
#'   }
#'   If \code{FALSE}, the function will instead use:
#'   \itemize{
#'     \item \strong{Sample mean:}
#'           \deqn{\hat{\beta}(s) = \frac{1}{n} \sum_{i=1}^n \tilde{Y}_i(s)}
#'           where \eqn{\tilde{Y}_i(s)} is a smoothed version of the observed functional response.
#'     For this sample mean option, if there exists NA in \code{data}, will perform imputation to fill in those values.
#'     }
#'   Default is \code{FALSE}.
#' @param alpha Significance level for SCB. Default is 0.05.
#' @param outcome A character string specifying the name of the outcome variable used in the model. Must be provided if choose \code{method = "multiplier"}.
#' @param time A character string specifying the name of the time variable used in the model.
#' @param range A numeric vector specifying the range or grid of time points interested. Defaut is the whole time points from the fitted model.
#' @param group_name A character vector that specifies the names of grouping scalar variables to analyze.
#'   Default is \code{NULL}, representing the reference group if \code{method = "cma"}.
#'   For \code{method = "multiplier"}, \code{NuULL} is not allowed.
#' @param group_value A numeric vector that specifies the corresponding values of the variables listed in \code{group_name}.
#'   Each entry in \code{group_value} should match the respective variable in \code{group_name}.
#'   Character/factor is not allowed.
#'   Default is \code{NULL}, representing the reference group if \code{"cma"}.
#'   For \code{method = "multiplier"}, \code{NuULL} is not allowed.
#' @param subject A character string specifying the name of the subject variable.
#' @param nboot An integer specifying the number of bootstrap samples used to construct the confidence bands. Default is 10,000 for cma, 5000 for Multiplier Bootstrap.
#' @param method_SD Method for SD estimation: "t" or "regular". Default is "t".
#' @param weights Multiplier type: "rademacher", "gaussian", or "mammen". Default is "rademacher"."
#'
#' For CMA bands (\code{method = "cma"}), please provide a pre-fitted functional regression model
#' For Multiplier Bootstrap bands (\code{method = "multiplier"}), the original input data is required.
#'
#' @returns A list containing:
#'   \item{yhat}{Estimated mean function for the group of interest.}
#'   \item{time}{The time points used.}
#'   \item{se_hat}{Standard errors of the estimated means.}
#'   \item{scb_low}{Lower bound of the simultaneous confidence band.}
#'   \item{scb_up}{Upper bound of the simultaneous confidence band.}
#'   \item{type}{A character description of the output type.}
#'
#'@importFrom tidyr pivot_wider
#'@importFrom dplyr select mutate all_of %>%
#'@importFrom refund fpca.face
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
#' # CMA approach
#' results <- SCB_functional_outcome(data_df = pupil, object = fosr_mod, method = "cma", fitted = TRUE,
#'                                   est_mean = TRUE, outcome = "percent_change", time = "seconds",
#'                                   group_name = "use", group_value = 1, subject = "id")
#'
#'
#' # multiplier bootstrap
#' results <- SCB_functional_outcome(data_df = pupil, object = fosr_mod,
#'                                   method = "multiplier", fitted = TRUE,
#'                                   est_mean = TRUE, outcome = "percent_change",
#'                                   time = "seconds", group_name = "use",
#'                                   group_value = 1, subject = "id")
#' results <- SCB_functional_outcome(data_df = pupil, object = fosr_mod,
#'                                   method = "multiplier", fitted = TRUE,
#'                                   est_mean = FALSE, outcome = "percent_change",
#'                                   time = "seconds", group_name = "use",
#'                                   group_value = 1, subject = "id")
#'
#'
#' @export
SCB_functional_outcome = function(data_df, object = NULL, method, fitted = TRUE, est_mean = TRUE, alpha = 0.05,
                                  outcome, time, range = NULL, group_name = NULL, group_value = NULL,
                                  subject, nboot = NULL, method_SD = "t", weights = "rademacher") {
  if (is.null(data_df)) {
    stop("`data_df` must be provided.")
  }

  if (!is.null(nboot)){
    if (!is.numeric(nboot) || nboot <= 0 || nboot %% 1 != 0){
      stop("`nboot` must be a positive integer.")
    }
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

  if (method == "cma") {
    if (is.null(object)) {
      stop("`object` must be provided when `method = 'cma'`.")
    }
    return(cma(data_df = data_df, object = object, fitted = fitted, alpha = alpha, time = time, range = range, group_name = group_name, group_value = group_value, subject = subject, nboot = nboot))
  }

  if (method != "multiplier") {
    stop("`method` must be either 'cma' or 'multiplier'.")
  }

  # method == "multiplier"
  if (is.null(outcome)) {
    stop("`outcome` must be provided.")
  }
  if (is.null(time)) {
    stop("`time` must be provided.")
  }
  if (is.null(subject)) {
    stop("`subject` must be provided.")
  }

  # Identify covariates
  vars <- c(outcome, time)
  if (!is.null(subject)) {
    vars <- c(vars, subject)
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

  # keep_rows <- rep(TRUE, nrow(data))
  # if (!is.null(group_name) && !is.null(group_value)){
    # for (i in seq_along(group_name)) {
      # var <- group_name[i]
      # val <- group_value[i]

      # if (!var %in% names(data)) {
        # stop(paste0("Variable '", var, "' not found in data."))
      # }

      # if (is.character(data[[var]])) {
        # stop(paste0("The variable '", var, "' is of type character. ",
                    # "Please convert it to a numerical variable and consider refit your functional object if necessary."))
      # }

      # if (is.factor(data[[var]])) {
        # stop(paste0("The variable '", var, "' is of type factor. ",
                    # "Please convert it to a numerical variable and refit your functional object if necessary."))
      # }else {  # numeric
        # unique_vals <- unique(data[[var]])
        # if (!(val %in% unique_vals)) {
          # stop(paste0("Value '", val, "' not found in numeric variable '", var, "'. "))
        # }
      # }
    # }

    # for (i in seq_along(group_name)) {
      # var <- group_name[i]
      # val <- group_value[i]
      # keep_rows <- keep_rows & (data[[var]] == val)
    # }
  # }else{
    # stop("Must provide both group variable name and value for multiplier bootstrap!")
  # }

  # Compute SCB
  if (!is.null(object)) {

    df <- mean_response_predict(data_df, object, fitted = fitted, time, range, group_name = group_name, group_value = group_value, subject)
    pred_df <- df$pred_df

    # impute NA if necessary
    if (anyNA(data_df)){
      data_wide <- data_df %>%
        select(all_of(c(subject, time, outcome))) %>%
        pivot_wider(
          names_from = !!sym(time),
          values_from = !!sym(outcome),
          names_prefix = "time_"
        ) %>%
        select(-all_of(subject)) # every row is a subject

      fpc_obj <- fpca.face(as.matrix(data_wide))

      data_df <- data_df %>%
        mutate(Yhat = as.vector(t(fpc_obj$Yhat)),
               !!sym(outcome) := ifelse(is.na(.data[[outcome]]), Yhat, .data[[outcome]])
        ) %>%
        select(-Yhat)

    }

    # Subset outcome values
    # if (is.null(range)) {
    # Y_samples <- data[[outcome]][keep_rows]
    Y_samples <- data_df[[outcome]][order(data_df[[subject]], data_df[[time]])]
    s_pred <- unique(data_df[[time]])
    n_time <- length(s_pred)
    # } else {
    # df <- data[data[[time]] %in% sort(unique(range)), ]
    # keep_rows_sub <- keep_rows[data[[time]] %in% sort(unique(range))]
    # Y_samples <- df[[outcome]][keep_rows_sub]
    # Y_samples <- df[[outcome]]
    # s_pred <- unique(range)
    # n_time <- length(s_pred)
    # }

    # Construct outcome matrix (T x N)
    n_subject <- length(Y_samples) / n_time
    Y_mat <- matrix(Y_samples, nrow = n_time, ncol = n_subject)
    # Y_mat[is.na(Y_mat)] <- 0  # replace NA with 0

    if(est_mean){
      thres <- SCB_dense(A = Y_mat, mean_A = pred_df$mean, alpha = alpha,
                         Mboots = nboot, method = method_SD, weights = weights, SCB = FALSE)
      scb_up = pred_df$mean + thres*pred_df$se
      scb_low = pred_df$mean - thres*pred_df$se
      # return index, scb_up, scb_low
      results <- list(
        yhat = pred_df$mean,
        se_hat = pred_df$se,
        scb_low = scb_low,
        scb_up = scb_up,
        type = "Dense Confidence Interval"
      )
    }else{
      # Check if exists zero entries
      if (s_pred[1] == 0 && all(Y_mat[1, ] == 0)) {
        Y_mat <- Y_mat[-1, , drop = FALSE]
        s_pred <- s_pred[-1]
        yhat <- pred_df$mean[-1]
        se_hat <- pred_df$se[-1]
      }else{
        yhat <- pred_df$mean
        se_hat <- pred_df$se
      }
      # Identify rows where all entries are zero, if exist, throw an error.
      # zero_rows <- rowSums(Y_mat == 0) == ncol(Y_mat)
      # Y_mat <- Y_mat[!zero_rows, ]
      zero_rows <- rowSums(Y_mat == 0) == ncol(Y_mat)
      if(any(zero_rows)){
        stop(paste0("Expected no rows where all entries are zero in `data_df` for `est_mean = 'FALSE'`."))
      }
      thres <- SCB_dense(A = Y_mat, alpha = alpha,
                            Mboots = nboot, method = method_SD, weights = weights, SCB = FALSE)
      scb_up = yhat + thres*se_hat
      scb_low = yhat - thres*se_hat
      # return index, scb_up, scb_low
      results <- list(
        yhat = yhat,
        se_hat = se_hat,
        scb_low = scb_low,
        scb_up = scb_up,
        type = "Dense Confidence Interval"
      )
    }
  }else{
    print("No Functional Regression Object provided, will only compute an overall SCB for the outcome regardless of the group specified.")

    # impute NA if necessary
    if (anyNA(data_df)){
      data_wide <- data_df %>%
        select(all_of(c(subject, time, outcome))) %>%
        pivot_wider(
          names_from = !!sym(time),
          values_from = !!sym(outcome),
          names_prefix = "time_"
        ) %>%
        select(-all_of(subject)) # every row is a subject

      fpc_obj <- fpca.face(as.matrix(data_wide))

      data_df <- data_df %>%
        mutate(Yhat = as.vector(t(fpc_obj$Yhat)),
               !!sym(outcome) := ifelse(is.na(.data[[outcome]]), Yhat, .data[[outcome]])
        ) %>%
        select(-Yhat)

    }

    if (is.null(range)) {
      # Y_samples <- data[[outcome]][keep_rows]
      Y_samples <- data_df[[outcome]][order(data_df[[subject]], data_df[[time]])]
      s_pred <- unique(data_df[[time]])
      n_time <- length(s_pred)
    } else {
      df <- data_df[data_df[[time]] %in% sort(unique(range)), ]
      # keep_rows_sub <- keep_rows[data[[time]] %in% sort(unique(range))]
      # Y_samples <- df[[outcome]][keep_rows_sub]
      Y_samples <- df[[outcome]][order(data_df[[subject]], data_df[[time]])]
      s_pred <- unique(range)
      n_time <- length(s_pred)
    }

    # Construct outcome matrix (T x N)
    n_subject <- length(Y_samples) / n_time
    Y_mat <- matrix(Y_samples, nrow = n_time, ncol = n_subject)

    # Check if exists zero entries
    if (s_pred[1] == 0 && all(Y_mat[1, ] == 0)) {
      Y_mat <- Y_mat[-1, , drop = FALSE]
      s_pred <- s_pred[-1]
    }

    # Identify rows where all entries are zero, if exist, throw an error.
    # zero_rows <- rowSums(Y_mat == 0) == ncol(Y_mat)
    # Y_mat <- Y_mat[!zero_rows, ]
    zero_rows <- rowSums(Y_mat == 0) == ncol(Y_mat)
    if(any(zero_rows)){
      stop(paste0("Expected no rows where all entries are zero in `data_df`."))
    }
    results <- SCB_dense(A = Y_mat, mean_A = NULL, alpha = alpha,
                         Mboots = nboot, method = method_SD, weights = weights)
  }

  results$time <- s_pred
  return(results)
}

