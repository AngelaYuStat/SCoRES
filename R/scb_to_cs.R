#' Construct Inverse Confidence Sets from Simultaneous Confidence Bands
#'
#' This function constructs inverse confidence sets (CS) from simultaneous confidence bands (SCB),
#' allowing visualization and containment check of level sets of true or estimated functions.
#'
#' @param scb_up A matrix or array containing the upper simultaneous confidence interval.
#'               Must have the same dimensions as defined by the grid from `x1` and `x2`.
#' @param scb_low A matrix or array containing the lower bounds of the simultaneous confidence bands.
#'               Must have the same dimensions as defined by the grid from `x1` and `x2`.
#' @param levels A list of scalers for different levels or matrix containing interval sets to construct the confidence sets.
#' @param true_mean Optional matrix of the true mean function (for simulation study or evaluation purposes).
#'                  Should have the same dimension as `scb_up`.
#' @param est_mean Optional matrix of the estimated mean function, used for plotting if `true_mean` is not available.
#'                  Should have the same dimension as `scb_up`.
#' @param x1 A numeric vector of coordinates for the first dimension (e.g., time or x-axis).
#' @param x2 A numeric vector of coordinates for the second dimension (e.g., space or y-axis).
#' @param type A character string specifying the type of inverse set to construct if levels are not a matrix.
#'        Choices are `"upper"`, `"lower"`, `"two-sided"` or `"interval"`.
#'        Notice that `"two-sided"` type is not available for plotting (\code{return_plot = TRUE}).
#' @param return_contain_only Logical. If `TRUE`, only return a matrix/logical map indicating whether the level set is contained.
#' @param return_plot Logical. If `TRUE`, return a ggplot object for visualization.
#' @param xlab A character for the name of the x axis for the returned ggplot object.
#' @param ylab A character for the name of the y axis for the returned ggplot object.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{levels}{A vector (or matrix) of threshold levels used to define the confidence sets.}
#'   \item{U_in}{(Optional) A list of logical matrices indicating whether each point is within the conservative inner confidence set for each level. Returned only when \code{return_contain_only = FALSE} and \code{type != "two-sided"}.}
#'   \item{U_out}{(Optional) A list of logical matrices indicating whether each point is within the liberal outer confidence set for each level. Returned only when \code{return_contain_only = FALSE} and \code{type != "two-sided"}.}
#'   \item{L_out}{(Two-sided only) A list of logical matrices indicating lower bound containment (for \code{type = "two-sided"}).}
#'   \item{U_out}{(Two-sided only) A list of logical matrices indicating upper bound containment (for \code{type = "two-sided"}).}
#'   \item{contain_individual}{A logical vector indicating whether the true mean is fully contained within each level's confidence set. Returned only if \code{true_mean} is provided.}
#'   \item{contain_all}{A single logical value indicating whether the true mean is contained in all levels' confidence sets.}
#'   \item{plot_cs}{(Optional) A list of ggplot2 objects for visualizing the confidence sets, returned when \code{return_plot = TRUE}. Includes both a combined plot and individual plots per level.}
#' }
#'
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
#' result <- SCB_linear_outcome(df_fit = df, model = model, grid_df = grid)
#' scb_to_cs(result$UpperBound, result$LowerBound, c(-1, -0.5, 0.5, 1),
#' x1 = grid$x1, x2 = grid$x2, est_mean = results$Mean)
#'
scb_to_cs = function(scb_up, scb_low, levels, true_mean = NULL,est_mean = NULL, x1, x2, type = "upper", return_contain_only = F, return_plot = F, xlab = NULL, ylab = NULL)
{
  if(return_plot){
    p_para <- list(xlab = xlab, ylab = ylab)
    if(is.null(true_mean)){
      pl_together = plot_cs(SCB = list(scb_up = scb_up, scb_low = scb_low),
                            levels = levels, type = type, x = x1, y = x2, mu_hat = est_mean,
                            together = T, xlab = p_para$xlab,ylab = p_para$ylab)
      pl = plot_cs(SCB = list(scb_up = scb_up, scb_low = scb_low),
                   levels = levels, type = type, x = x1, y = x2, mu_hat = est_mean,
                   together = F, xlab = p_para$xlab,ylab = p_para$ylab)
      pl = list(pl_together, pl)
    }else{
      pl_together = plot_cs(SCB = list(scb_up = scb_up, scb_low = scb_low),
                            levels = levels, type = type, x = x1, y = x2,mu_true = true_mean,
                            mu_hat = est_mean, together = T,xlab = p_para$xlab,ylab = p_para$ylab)
      pl = plot_cs(SCB = list(scb_up = scb_up, scb_low = scb_low), levels = levels, type = type,
                   x = x1, y = x2,mu_true = true_mean, mu_hat = est_mean,
                   together = F,
                   xlab = p_para$xlab,ylab = p_para$ylab)
      pl = list(pl_together, pl)
    }
  }else{
    pl = NULL
  }
  contain_v = c()
  in_list = list()
  out_list = list()
  j = 1
  if(is.vector(levels)){
    if(type == "upper"){
      for(c in levels){
        in_set = scb_low >= c
        out_set = scb_up >= c
        if(!return_contain_only){
          in_list[[j]] = in_set
          out_list[[j]]= out_set
        }
        if(!is.null(true_mean)){
          true_set = true_mean >= c
          contain_v = c(contain_v,incl_f(in_set, true_set) & incl_f(true_set, out_set))
        }
        j = j + 1
      }
    }else if(type == "lower"){
      for(c in levels){
        in_set = scb_up <= c
        out_set = scb_low <= c
        if(!return_contain_only){
          in_list[[j]] = in_set
          out_list[[j]]= out_set
        }
        if(!is.null(true_mean)){
          true_set = true_mean <= c
          contain_v = c(contain_v,incl_f(in_set, true_set) & incl_f(true_set, out_set))
        }
        j = j + 1
      }
    }else if(type == "two-sided"){
      out_set_low_list = list()
      out_set_up_list = list()
      for(c in levels){
        out_set_low = scb_low <= c
        out_set_up = scb_up >= c
        if(!return_contain_only){
          out_set_low_list[[j]] = out_set_low
          out_set_up_list[[j]]= out_set_up
        }
        if(!is.null(true_mean)){
          true_set_low = true_mean <= c
          true_set_up = true_mean >= c
          contain_v = c(contain_v,incl_f(true_set_low, out_set_low) & incl_f(true_set_up, out_set_up))
        }
        j = j + 1
      }
      return(list(levels = levels, L_out = out_set_low_list, U_out = out_set_up_list,
                  contain_individual = contain_v, contain_all = all(contain_v), plot_cs = plot_cs))
    }

  }else if(type == "interval"){# When we have interval level sets
    l_dim = dim(levels)
    for(i in 1:l_dim[1]){
      c = c(levels[i,])
      in_set = scb_low >= c$low & scb_up <= c$up
      out_set = scb_up >= c$low & scb_low <= c$up
      if(!return_contain_only){
        in_list[[i]] = in_set
        out_list[[i]]= out_set
      }
      if(!is.null(true_mean)){
        true_set = true_mean >= c$low & true_mean <= c$up
        contain_v = c(contain_v,incl_f(in_set, true_set) & incl_f(true_set, out_set))
      }
    }
  }else{
    stop("Type must be chosen between 'upper', 'lower', 'two-sided' or 'interval'.")
  }
  return(list(levels = levels, U_in = in_list, U_out = out_list,
              contain_individual = contain_v, contain_all = all(contain_v), plot_cs = pl))
}

#' Test inclusion relationship between two logical vectors
#'
#' This function tests whether all \code{TRUE} positions in vector \code{A} are also \code{TRUE} in vector \code{B},
#' which is equivalent to testing whether \code{A} is element-wise included in \code{B}.
#'
#' @param A A logical vector.
#' @param B A logical vector of the same length as \code{A}.
#'
#' @return A logical value: \code{TRUE} if all \code{TRUE} values in \code{A} are also \code{TRUE} in \code{B}; otherwise \code{FALSE}.
#' @keywords internal
#'
#' @examples
#' # Used internally by scb_to_cs
#'
incl_f <- function(A, B) {
  min(B - A) >= 0 #if B includes A
}
# Make a data object

