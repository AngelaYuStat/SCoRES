#' Construct Simultaneous Confidence Bands (SCB) For One Dimensional Functional Data Using Specified Methods
#'
#' This function builds simultaneous confidence bands through two distinct approaches..
#'
#' @param data A functional data frame, matrix or tibble.
#' The input data must have column names, and should contain the functional outcome, time and subject
#' @param object A fitted Function-on-Scalar Regression (FoSR) object (e.g., from mgcv::bam()/mgcv::gam()). Default is NULL.
#' If provided when \code{method = "wild"}, the mean estimation from the output of functional regression model wil be used to calculate the SCBs, instead of sample mean
#' @param method A character string specifying the approach:
#'   \itemize{
#'     \item \code{"cma"} - Correlation and Multiplicity Adjusted (CMA) confidence bands
#'           via parametric approach (requires a fitted functional regression model)
#'     \item \code{"wild"} - Dense confidence bands via Multiplier-t Bootstrap method
#'           (requires raw data)
#'   }
#' @param est_mean Logical. Whether to use the fitted mean from a functional regression model.
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
#'   }
#'   Default is \code{FALSE}.
#' @param alpha Significance level for SCB. Default is 0.05.
#' @param outcome A character string specifying the name of the outcome variable used in the model. Must be provided if choose \code{method = "wild"}.
#' @param time A character string specifying the name of the time variable used in the model.
#' @param range A numeric vector specifying the range or grid of time points interested. Defaut is the whole time points from the fitted model.
#' @param group_name A character vector that specifies the names of grouping scalar variables to analyze.
#' @param group_value A vector (numeric or character) that specifies the corresponding values of the variables listed in \code{group_name}.
#'   Each entry in \code{group_value} should match the respective variable in \code{group_name}.
#' @param subject A character string specifying the name of the subject-level random effect variable.
#' @param nboot An integer specifying the number of bootstrap samples used to construct the confidence bands. Default is 10,000 for cma, 5000 for wild bootstrap.
#' @param method_SD Method for SD estimation: "t" or "regular". Default is "t".
#' @param weights Multiplier type: "rademacher", "gaussian", or "mammen". Default is "rademacher"."
#'
#' For CMA bands (\code{method = "cma"}), please provide a pre-fitted functional regression model
#' For wild bootstrap bands (\code{method = "wild"}), the original input data is required.
#'
#' @returns A list containing:
#'   \item{yhat}{Estimated mean function for the group of interest.}
#'   \item{time}{The time points used.}
#'   \item{se_hat}{Standard errors of the estimated means.}
#'   \item{scb_low}{Lower bound of the simultaneous confidence band.}
#'   \item{scb_up}{Upper bound of the simultaneous confidence band.}
#'   \item{type}{A character description of the output type.}
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
#' # CMA approach
#' results <- SCB_functional_outcome(data = ccds, object = fosr_mod, method = "cma", est_mean = TRUE,
#'                        outcome = "percent_change", time = "seconds",
#'                        group_name = "use", group_value = 1, subject = "subject")
#'
#'
#' # Wild bootstrap
#' results <- SCB_functional_outcome(data = ccds, object = fosr_mod, method = "wild",
#'                        est_mean = TRUE, outcome = "percent_change",
#'                        time = "seconds", group_name = "use", group_value = 1, subject = "subject")
#'
#'
#' @export
SCB_functional_outcome = function(data, object = NULL, method, est_mean = FALSE, alpha = 0.05,
                                  outcome, time, range = NULL, group_name, group_value,
                                  subject = NULL, nboot = NULL, method_SD = "t", weights = "rademacher") {
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

  if (method == "cma") {
    if (is.null(object)) {
      stop("Must provide a functional regression model object!")
    }
    return(cma(data, object, alpha = alpha, time = time, range = range, group_name = group_name, group_value = group_value, subject = subject, nboot = nboot))
  }

  if (method != "wild") {
    stop("Must choose between 'cma' and 'wild' as the method.")
  }

  # method == "wild"
  if (is.null(outcome)) {
    stop("Must provide the outcome variable name.")
  }
  if (is.null(time)) {
    stop("Must provide the time variable name.")
  }

  # Identify covariates
  vars <- c(outcome, time)
  if (!is.null(subject)) {
    vars <- c(vars, subject)
  }

  if (length(group_name) != length(group_value)) {
    stop("The length of 'group_name' and 'group_value' must be the same.")
  }


  for (i in seq_along(group_name)) {
    var <- group_name[i]
    val <- group_value[i]

    if (!var %in% names(data)) {
      stop(paste0("Variable '", var, "' not found in data."))
    }

    if (is.character(data[[var]])) {
      stop(paste0("The variable '", var, "' is of type character. ",
                  "Please convert it to a factor to ensure consistent comparison behavior."))
    }

    if (is.factor(data[[var]])) {
      if (!val %in% levels(data[[var]])) {
        stop(paste0("The value '", val, "' is not a level of factor variable '", var, "'."))
      }
    }
  }

  keep_rows <- rep(TRUE, nrow(data))
  for (i in seq_along(group_name)) {
    var <- group_name[i]
    val <- group_value[i]
    keep_rows <- keep_rows & (data[[var]] == val)
  }


  # Subset outcome values
  if (is.null(range)) {
    Y_samples <- data[[outcome]][keep_rows]
    s_pred <- unique(data[[time]])
    n_time <- length(s_pred)
  } else {
    df <- data[data[[time]] %in% sort(unique(range)), ]
    keep_rows_sub <- keep_rows[data[[time]] %in% sort(unique(range))]
    Y_samples <- df[[outcome]][keep_rows_sub]
    s_pred <- unique(range)
    n_time <- length(s_pred)
  }

  # Construct outcome matrix (T x N)
  n_subject <- length(Y_samples) / n_time
  Y_mat <- matrix(Y_samples, nrow = n_time, ncol = n_subject)
  Y_mat[is.na(Y_mat)] <- 0  # replace NA with 0

  # Compute SCB
  if (est_mean && !is.null(object)) {
    df <- mean_response_predict(data, object, time, range, group_name = group_name, group_value = group_value, subject)
    pred_df <- df$pred_df
    results <- SCB_dense(A = Y_mat, mean_A = pred_df$mean, alpha = alpha,
                         Mboots = nboot, method = method_SD, weights = weights)
  } else {
    # Identify rows where all entries are zero, if exist, then drop those rows.
    zero_rows <- rowSums(Y_mat == 0) == ncol(Y_mat)
    Y_mat <- Y_mat[!zero_rows, ]
    results <- SCB_dense(A = Y_mat, mean_A = NULL, alpha = alpha,
                         Mboots = nboot, method = method_SD, weights = weights)
  }

  results$time <- s_pred
  return(results)
}

