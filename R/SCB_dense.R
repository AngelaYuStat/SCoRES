#' Construct Simultaneous Confidence Bands (SCB) for Dense Functional Data
#'
#' @param A A data array of dimension (D₁, D₂, ..., N), where N is the number of repetition/subjects. There should be no NA in A.
#' @param mean_A Optional array of same shape as \code{A[,,1]}, representing the estimated mean of the data.
#' @param alpha Significance level for SCB. Default is 0.05.
#' @param Mboots Number of bootstrap replications. Default is 5000.
#' @param method Method for SD estimation: "t" or "regular". Default is "t".
#' @param weights Multiplier type: "rademacher", "gaussian", or "mammen". Default is "rademacher"."
#' @param SCB Logical value for whether to calculate the SCB or not. Default is TRUE.
#'
#' @returns If `SCB = TRUE`, returns a list containing:
#' \item{yhat}{Estimated mean function for the group of interest.}
#' \item{se_hat}{Standard errors of the estimated means.}
#' \item{scb_low}{Lower bound of the simultaneous confidence band.}
#' \item{scb_up}{Upper bound of the simultaneous confidence band.}
#' \item{type}{A character description of the output type.}
#'
#' If `SCB = FALSE`, returns:
#' \item{thres}{The alpha quantile estimated by Multiplier Bootstrap}
#'
#' @importFrom stats var
#'
#' @keywords internal
#' # Used internally by SCB_functional_outcome
#' @examples
#'
#' data(ccds)
#' ccds_fpca <- prepare_ccds_fpca(ccds)
#'
#' Y_samples_use <- ccds_fpca$percent_change[ccds_fpca$use == 1]
#' n_time <- 401
#' n_subject_use <- length(Y_samples_use) / n_time
#' # transform to matrix (T, N)
#' Y_mat_use <- matrix(Y_samples_use, nrow = n_time, ncol = n_subject_use)
#' # using sample mean for bootstrap
#' Y_mat_use[is.na(Y_mat_use)] <- 0
#' SCB_dense(A = Y_mat_use[-1,])
#'
SCB_dense = function(A, mean_A = NULL, alpha = 0.05, Mboots  = NULL,
                     method  = "t", weights = "rademacher", SCB = TRUE){
  # require(SIRF)
  # Get the number of repeats and dimension
  # if(any(is.na(A))){A[is.na(A)] <- 0}
  dimA = dim(A)
  N = dimA[length(dimA)]
  D = length(dimA) - 1
  # getting the estimated mean
  if(is.null(mean_A)){mean_A <- apply(A, 1:D, mean)}
  mean_A_array = array(rep(mean_A, N), dim = c(dimA[1:D], N))
  if(is.null(Mboots)){Mboots <- 5e3}
  # Get the threshold
  thres = MultiplierBootstrap(sqrt(N/(N - 1))*(A-mean_A_array), alpha = alpha,
                              Mboots  = Mboots, method = method, weights = weights)$q

  if(SCB == TRUE){
    # getting the estimated sd
    sd_A = sqrt(apply(A, 1:D, var))
    scb_up = mean_A + thres*sd_A/sqrt(N)
    scb_low = mean_A - thres*sd_A/sqrt(N)
    # return index, scb_up, scb_low
    dense_df <- list(
      yhat = mean_A,
      se_hat = sd_A,
      scb_low = scb_low,
      scb_up = scb_up,
      type = "Dense Confidence Interval"
    )
    return(dense_df)
  }else{
    return(thres)
  }
}
