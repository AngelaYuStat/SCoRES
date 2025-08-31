#' Functional Outcome Regression Prediction with Group-Specific Inference
#'
#' This function is an internal function for constructing SCBs for functional data.
#'
#' @param data_df A functional data frame that contain both name and values for
#' variables including functional outcome, domain (e.g. time) and ID (e.g. subject names)
#' used to fit `object`.
#' @param object A fitted Function-on-Scalar Regression (FoSR) model object
#' (e.g., from mgcv::gam()/mgcv::bam()).
#' @param fitted Logical. Whether to estimate the simultaneous confidence bands
#' for fitted mean function or fitted parameter function
#'   \itemize{
#'     \item \code{TRUE} - Estimate the simultaneous confidence bands
#'     for fitted mean outcome function.
#'     \item \code{FALSE} - estimate the simultaneous confidence bands
#'     for fitted parameter function.
#'     }
#'   Default is \code{TRUE}.
#' @param outcome A character string specifying the name of the outcome variable
#' used in the model.
#' @param domain A character string specifying the name of the domain variable
#' (e.g. time) used in the model.
#' @param subset An atomic character vector (e.g., c("user = 1", "age = 30"))
#' specified the target function for constructing the SCB.
#' Each element must be of the form <name> = <value>, where <name> is the name
#' of a scalar grouping variable and <value> is the desired value.
#' Whitespace is ignored. Binary or categorical character variable should be
#' transformed into numeric. Factor is not allowed here because if the input
#' data contains factor variables, they will be automatically expanded into
#' dummy (indicator) variables when constructing the design matrix, and
#' the resulting variable names may differ from the original factor names.
#' Default is \code{NULL}, representing the reference group.
#' @param id A character string specifying the name of the ID variable.
#'
#' @returns A list containing the following elements:
#' \describe{
#'   \item{s_pred}{Numeric vector of sorted unique domain used for prediction}
#'   \item{pred_df}{Data frame with prediction results, containing:
#'     \itemize{
#'       \item \code{mean}: Predicted mean values
#'       \item \code{se}: Standard errors
#'     }
#'   }
#'   \item{lpmat}{Linear predictor matrix (design matrix)
#'   used for confidence interval calculations}
#'   \item{mod_coef}{Vector of model coefficients for selected group}
#'   \item{mod_cov}{Variance-covariance matrix corresponding to
#'   the selected group coefficients}
#' }
#'
#' @importFrom stats formula model.frame vcov
#' @importFrom dplyr mutate
#' @importFrom tibble as_tibble
#' @importFrom stats predict
#' @importFrom magrittr %>%
#'
#' @examples
#' library(mgcv)
#' data(pupil)
#' pupil_fpca <- prepare_pupil_fpca(pupil)
#'
#' fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
#'   s(seconds, by = use, k=30, bs = "cr") +
#'   s(seconds, by = age, k = 30, bs = "cr") +
#'   s(seconds, by = gender, k = 30, bs = "cr") +
#'   s(id, by = Phi1, bs="re") +
#'   s(id, by = Phi2, bs="re")+
#'   s(id, by = Phi3, bs="re") +
#'   s(id, by = Phi4, bs="re"),
#'   method = "fREML", data = pupil_fpca, discrete = TRUE)
#'
#' results <- mean_response_predict(pupil_fpca, fosr_mod, fitted = TRUE,
#' outcome = "percent_change", domain = "seconds", subset = c("use = 1"), id = "id")
#'
#' @export
#'
mean_response_predict = function(data_df, object, fitted = TRUE, outcome, domain, subset = NULL, id){

  mod_coef <- object$coefficients
  mod_cov <- vcov(object) # containing all variance-covariance info for object

  coef_names <- names(mod_coef)
  # get index of the initial terms
  pattern <- paste0("^\\(Intercept\\)$|^s\\(", domain, "\\)\\.[0-9]+$")
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

  if (!is.null(subset)){
    # Identify covariates
    #vars <- c(outcome, domain, id)

    m <- regexec("^\\s*([^=]+?)\\s*=\\s*(.+)$", subset)
    res <- regmatches(subset, m)

    ok <- lengths(res) > 0
    if (any(!ok)) {
      warning("Some elements did not match the required <name> = <value> pattern.
              Please check your input for `subset`.")
    }

    group_name <- sapply(res, `[`, 2)
    group_value <- sapply(res, `[`, 3)
    #safely converge to numeric if possible
    group_value <- type.convert(trimws(group_value), as.is = TRUE)

    if (!is.null(group_name) && !is.null(group_value)){
      for (i in seq_along(group_name)) {
        var <- group_name[i]
        val <- group_value[i]

        if (!var %in% names(data_df)) {
          stop(paste0("Variable '", var, "' not found in `data_df`."))
        }

        if (is.character(data_df[[var]])) {
          stop(paste0("The variable '", var, "' is of type character. ",
                      "Please convert it to a numerical variable and consider refit
          your functional object if necessary."))
        }

        if (is.factor(data_df[[var]])) {
          stop(paste0("The variable '", var, "' is of type factor. ",
                      "Please convert it to a numerical variable and consider refit
          your functional object if necessary."))
        }

        if (is.numeric(data_df[[var]])) { # numeric
          #unique_vals <- unique(data_df[[var]])
          #if (!(val %in% unique_vals)) {
            #stop(paste0("Value '", val, "' not found in numeric variable '", var, "'. "))
          #}
          if(!is.numeric(val)) stop("Value '", val, "' is not numeric.")
          df_pred[[var]] <- val
        }else{
          stop(paste0("The variable '", var, "' is not of type numeric. ",
                      "Please convert it to a numerical variable and consider refit
          your functional object if necessary."))
        }
      }
      # for group interested
      groups_idx <- unlist(lapply(group_name, function(g) grep(g, coef_names)))
      if (fitted){
        idx <- sort(unique(c(idx, groups_idx)))
      }else{
        idx <- sort(unique(c(groups_idx)))
      }
    }else{
      stop("Must provide valid input for `subset`.
           Each element must be of the form <name> = <value>.")
    }
  }

  #if(!is.null(id)){
    #df_pred[id] <- as.character(model.frame(object)[[id]][1])
  #}

  if (id %in% names(model.frame(object))) {
    df_pred[id] <- as.character(model.frame(object)[[id]][1])
  }

  s_pred <- sort(unique(model.frame(object)[[domain]]))

  df_pred <- df_pred[rep(1, length(s_pred)), ]
  df_pred[[domain]] <- s_pred

  # prepare design matrix
  lpmat <- predict(object, newdata=df_pred, se.fit=TRUE, type = "lpmatrix")

  # get mean response by groups with standard errors
  if (!fitted && !is.null(groups_idx)){
    # fitted = FALSE with group specified, fit mean outcome for linear combination
    # of parameter corresponding to the group specified without intercept
    lpmat_no_intercept <- lpmat[, -intercept_idx, drop = FALSE]
    mod_coef_no_intercept <- mod_coef[-intercept_idx]
    mod_cov_no_intercept <- mod_cov[-intercept_idx, -intercept_idx]
    pred_df <- df_pred %>%
      mutate(mean = c(lpmat_no_intercept %*% mod_coef_no_intercept),
             se = c(sqrt(diag(lpmat_no_intercept %*% mod_cov_no_intercept %*% t(lpmat_no_intercept)))))
  }else{
    # fitted = FALSE with no group specified, fit intercept
    # fitted = TRUE with group specified, fit mean outcome for the group specified
    # fitted = TRUE with no group specified, fit intercept
    pred_df <- df_pred %>% mutate(mean = c(lpmat %*% mod_coef),
                                  se = c(sqrt(diag(lpmat %*% mod_cov %*% t(lpmat)))))
  }
  lpmat <- lpmat[, idx]
  # extract means and variance for group interested
  mod_coef <- mod_coef[idx]
  mod_cov <- mod_cov[idx, idx]

  return(list(s_pred = s_pred, pred_df = pred_df, lpmat = lpmat,
              mod_coef = mod_coef, mod_cov = mod_cov))
}
#' Construct CMA Confidence Intervals via Parametric Method
#'
#' This function computes Correlation and Multiplicity Adjusted (CMA) confidence bands
#'  for a specified group in a functional outcome regression model
#' using parameter simulations approach with Gaussian multiplier bootstrap.
#'
#' @param data_df A functional data frame that contain both name and values for
#' variables including functional outcome, domain (e.g. time) and ID (e.g. subject names)
#' used to fit `object`.
#' @param object A fitted Function-on-Scalar Regression (FoSR) object
#' (e.g., from mgcv::gam()/mgcv::bam()).
#' @param fitted Logical. Whether to estimate the simultaneous confidence bands
#' for fitted mean function or fitted parameter function
#'   \itemize{
#'     \item \code{TRUE} - Estimate the simultaneous confidence bands
#'     for fitted mean outcome function.
#'     \item \code{FALSE} - estimate the simultaneous confidence bands
#'     for fitted parameter function.
#'     }
#'   Default is \code{TRUE}.
#' @param alpha Significance level for SCB. Default is 0.05.
#' @param outcome A character string specifying the name of the outcome variable
#' used in the model.
#' @param domain A character string specifying the name of the domain variable
#' (e.g. time) used in the model.
#' @param subset An atomic character vector (e.g., c("user = 1", "age = 30"))
#' specified the target function for constructing the SCB.
#' Each element must be of the form <name> = <value>, where <name> is the name
#' of a scalar grouping variable and <value> is the desired value.
#' Whitespace is ignored. Binary or categorical character variable should be
#' transformed into numeric. Factor is not allowed here because if the input
#' data contains factor variables, they will be automatically expanded into
#' dummy (indicator) variables when constructing the design matrix, and
#' the resulting variable names may differ from the original factor names.
#' Default is \code{NULL}, representing the reference group.
#' @param id A character string specifying the name of the ID variable.
#' @param nboot An integer specifying the number of bootstrap samples used to
#' construct the confidence bands. Default is 10,000.
#'
#' @returns A list containing:
#'   \item{mu_hat}{Estimated mean function for the group of interest.}
#'   \item{domain}{The domain used.}
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
#'   s(seconds, by = age, k = 30, bs = "cr") +
#'   s(seconds, by = gender, k = 30, bs = "cr") +
#'   s(id, by = Phi1, bs="re") +
#'   s(id, by = Phi2, bs="re")+
#'   s(id, by = Phi3, bs="re") +
#'   s(id, by = Phi4, bs="re"),
#'   method = "fREML", data = pupil_fpca, discrete = TRUE)
#'
#' results <- cma(pupil_fpca, fosr_mod, fitted = TRUE, outcome = "percent_change",
#'                domain = "seconds", subset = c("use = 1"), id = "id")
#'
cma = function(data_df, object, fitted = TRUE, alpha = 0.05, outcome, domain,
               subset = NULL, id, nboot = NULL){

  results <- mean_response_predict(data_df, object, fitted, outcome = outcome,
                                   domain = domain, subset = subset, id = id)

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

  return(list(
    mu_hat = pred_df$mean,
    domain = results$s_pred,
    se_hat = pred_df$se,
    scb_low = y_hat_LB_global,
    scb_up = y_hat_UB_global,
    type = "CMA Confidence Band"
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
#'     \item Sorted by ID and domain
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

