#' Construct Simultaneous Confidence Bands (SCB) for Dense Functional Data
#'
#' @param A An array of dimension (D₁, D₂, ..., N), where N is the number of repetitions.
#' @param mean_A Optional array of same shape as \code{A[,,1]}, representing the mean function.
#' @param alpha_level Significance level for SCB. Default is 0.05.
#'
#' @return A list with upper and lower SCB bounds and the threshold used.
#' @export
#'
#' @examples
#' A <- array(rnorm(100 * 50), dim = c(100, 50))
#' result <- SCB_dense(A)
#'
#' data(ccds_fpca)
#' data(ccds_cma)
#'
#' Y_samples <- ccds_fpca$percent_change
#' n_time <- length(unique(ccds_fpca$seconds))
#' n_subject <- length(Y_samples) / n_time
#' # transform to matrix (T, N)
#' Y_mat <- matrix(Y_samples, nrow = n_time, ncol = n_subject)
#' Y_mat[is.na(Y_mat)] <- 0 # replace NA
#' SCB_dense(Y_mat, ccds_cma$yhat)
#'
SCB_dense = function(A, mean_A = NULL, alpha_level = 0.05){
  # require(SIRF)
  # Get the number of repeats and dimension
  dimA = dim(A)
  N = dimA[length(dimA)]
  D = length(dimA) - 1
  # getting the estimated sd
  sd_A = sqrt(apply(A, 1:D, var))
  # getting the estimated mean
  if(is.null(mean_A)){mean_A <- apply(A, 1:D, mean)}
  mean_A_array = array(rep(mean_A, N), dim = c(dimA[1:D], N))
  # Get the threshold
  thres = MultiplierBootstrap(sqrt(N/(N - 1))*(A-mean_A_array), alpha = alpha_level)$q # the default is the multiplier-t with Rademacher multipliers, here alpha is a
  # construct the SCB
  scb_up = mean_A + thres*sd_A/sqrt(N)
  scb_low = mean_A - thres*sd_A/sqrt(N)
  # return index, scb_up, scb_low
  return(list(scb_up = scb_up, scb_low = scb_low, thres = thres))
}
