#' Construct Simultaneous Confidence Region for Excursion/Interval Sets from Simultaneous Confidence Bands
#'
#' This function constructs simultaneous confidence region for upper and lower
#' excursion sets, and interval sets from simultaneous confidence bands (SCB).
#' It allows estimation of inner and outer confidence region under
#' single or multiple thresholds. Visualization of the confidence region is also
#' included, along with a containment check for the coverage of true or estimated functions.
#'
#' @param scb_up A numeric vector (1D) or matrix (2D) containing
#' the upper simultaneous confidence interval.
#' @param scb_low A numeric vector (1D) or matrix (2D) containing
#' the lower bounds of the simultaneous confidence bands.
#' Dimensions of `scb_up` and `scb_low` must match.
#' @param levels A numeric vector or list of scalers for different levels or matrix
#' containing interval sets to construct the confidence sets.
#' If \code{type} = "upper", "lower", or "two-sided", `levels` should be a vector.
#' "upper" represents upper excursion sets, and "lower" represents lower excursion sets.
#' If "two-sided" option is chosen, will estimate only outer CSs for both upper
#' and lower excursion sets.
#' If \code{type = "interval"}, then \code{levels} should be a \code{list}
#' with two named elements: \code{low} and \code{up},
#' corresponding to the bounds of the interval \code{[low, up]}.
#' @param true_mean Optional matrix of the true mean function.
#' Should have the same dimension as `scb_up` and `scb_low`.
#' @param est_mean Optional matrix of the estimated mean function,
#' used for plotting if `true_mean` is not available.
#' Should have the same dimension as `scb_up` and `scb_low`.
#' @param x1 A numeric vector of coordinates for the first dimension used
#' for plotting the inner and outer confidence region. Default is NULL.
#' Dimension of `x1` must match the first dimension of `scb_up` and `scb_low`.
#' @param x2 A numeric vector of coordinates for the second dimension used
#' for plotting inner and outer confidence region. Default is NULL.
#' Dimension of `x1` must match the second dimension of `scb_up` and `scb_low`.
#' @param type A character string specifying the type of inverse set to construct
#' if levels are not a matrix. Choices are `"upper"`, `"lower"`, `"two-sided"`
#' or `"interval"`. Notice that `"two-sided"` and `"interval"` type is not available for
#' plotting (\code{return_plot = TRUE}).
#' @param return_contain_only Logical. If `TRUE`, only return a matrix/logical
#' map indicating which point is contained within two types of CSs across all levels.
#' @param return_plot Logical. If `TRUE`, return a ggplot object for visualizing
#' the inner and outer confidence region.
#' @param xlab A character for the name of the x axis used for plotting the inner
#' and outer confidence region. Default is NULL.
#' @param ylab A character for the name of the y axis used for plotting the inner
#' and outer confidence region. Default is NULL.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{levels}{A vector (or list) of threshold levels used to define the
#'   confidence sets. Same as the input `levels`.}
#'   \item{U_in}{(Optional) A list of logical matrices indicating whether each
#'   point is within the simultaneous inner confidence set for each level.
#'   Returned only when \code{return_contain_only = FALSE} and \code{type != "two-sided"}.}
#'   \item{U_out}{(Optional) A list of logical matrices indicating whether each
#'   point is within the simultaneous outer confidence set for each level.
#'   Returned only when \code{return_contain_only = FALSE} and \code{type != "two-sided"}.}
#'   \item{L_out}{(Two-sided only) A list of logical matrices
#'   indicating lower bound containment (for \code{type = "two-sided"} and \code{return_contain_only = FALSE}).}
#'   \item{U_out}{(Two-sided only) A list of logical matrices
#'   indicating upper bound containment (for \code{type = "two-sided"} and \code{return_contain_only = FALSE}).}
#'   \item{contain_individual}{A logical vector indicating
#'   whether the true mean is fully contained within each level's simultaneous
#'   inner and outer confidence region. Returned only if \code{true_mean} is provided.}
#'   \item{contain_all}{A single logical value indicating whether the true mean
#'   is contained in all levels' simultaneous inner and outer confidence region.
#'   Returned only if \code{true_mean} is provided.}
#'   \item{plot_cs}{(Optional) A list of ggplot2 objects for visualizing the
#'   SCBs and simultaneous confidence region across all levels,
#'   returned when \code{return_plot = TRUE}. Includes both a combined plot and
#'   individual plots per level.}
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
#' scb_to_cs(result$scb_up, result$scb_low, c(-1, -0.5, 0.5, 1),
#' x1 = grid$x1, x2 = grid$x2, est_mean = results$Mean)
#'
scb_to_cs = function(scb_up, scb_low, levels, true_mean = NULL, est_mean = NULL, x1 = NULL, x2 = NULL, type = "upper", return_contain_only = F, return_plot = F, xlab = NULL, ylab = NULL)
{
  if(is.null(scb_up)||is.null(scb_low)){
    stop("Must provide input for `scb_up` and `scb_low`.")
  }

  if(!is.numeric(scb_up) || !is.numeric(scb_low)) {
    stop("Values of `scb_up` and `scb_low` must be numeric.")
  }

  if (!all(sapply(list(scb_up, scb_low),
                  function(x) (is.atomic(x) && is.null(dim(x)) || is.matrix(x) || is.array(x))))) {
    stop("`scb_up` and `scb_low` must each be a vector, array or matrix.")
  }

  nd <- length(dim(scb_up))
  if (nd == 0L) {
    if(!(length(scb_up) == length(scb_low))) {
      stop("Dimensions of `scb_up` and `scb_low` must match.")
    }
    if(!is.null(true_mean)){
      if(!is.numeric(true_mean)) stop("Values of `true_mean` must be numeric.")
      if(!(length(scb_up) == length(true_mean))) {
        stop("Dimensions of `scb_up`, `scb_low` and `true_mean` must match.")
      }
      if (!all(sapply(list(true_mean),
                      function(x) (is.atomic(x) && is.null(dim(x))) || is.matrix(x) || is.array(x)))) {
        stop("`true_mean` must each be a vector, array or matrix.")
      }
    }
  }else{
    if(!identical(dim(scb_up), dim(scb_low))) {
      stop("Dimensions of `scb_up` and `scb_low` must match.")
    }
    if(!is.null(true_mean)){
      if(!is.numeric(true_mean)) stop("Values of `true_mean` must be numeric.")
      if(!identical(dim(scb_up), dim(true_mean))) {
        stop("Dimensions of `scb_up`, `scb_low` and `true_mean` must match.")
      }
      if (!all(sapply(list(true_mean),
                      function(x) (is.atomic(x) && is.null(dim(x))) || is.matrix(x) || is.array(x)))) {
        stop("`true_mean` must each be a vector, array or matrix.")
      }
    }
  }

  if(is.null(levels)) {
    stop("Must provide input for `levels`.")
  }else{
    if(type %in% c("upper", "lower", "two-sided")){
      if (!is.null(dim(levels)) && length(dim(levels)) == 1) {
        levels <- as.vector(levels)
      }
      if(!(is.atomic(levels) && is.null(dim(levels)))) stop("`levels` should be a vector if `type` = upper or lower.")
      if(!is.numeric(levels)) {
        stop("Values of `levels` must be numeric.")
      }
    }else if(type == "interval"){
      if(!is.list(levels)) stop("`levels` should be a list if `type` = interval")
      if (!all(c("low", "up") %in% names(levels))) {
        stop("`levels` must have elements named 'low' and 'up'.")
      }
      if(!is.numeric(levels$low)||!is.numeric(levels$up)) {
        stop("All elements in `levels` must be numeric.")
      }
    }else{
      stop("`type` must be chosen between 'upper', 'lower', 'two-sided' or 'interval'.")
    }
  }

  if(!is.logical(return_contain_only)) stop("`return_contain_only` must be logical.")
  if(!is.logical(return_plot)) stop("`return_plot` must be logical.")
  if(return_plot && type %in% c("two-sided", "interval")) stop("Current function doesn't support plotting for `two-sided` and `interval`.")

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
      if(return_contain_only && !is.null(true_mean)){
        return(list(levels = levels, contain_individual = contain_v, contain_all = all(contain_v), plot_cs = pl))

      }else if(!return_contain_only && !is.null(true_mean)){
        return(list(levels = levels, L_out = out_set_low_list, U_out = out_set_up_list,
                    contain_individual = contain_v, contain_all = all(contain_v), plot_cs = pl))
      }else if(!return_contain_only && is.null(true_mean)){
        return(list(levels = levels, L_out = out_set_low_list, U_out = out_set_up_list, plot_cs = pl))
      }else{
        return(list(levels = levels, plot_cs = pl))
      }
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
  }
  if(return_contain_only && !is.null(true_mean)){
    return(list(levels = levels, contain_individual = contain_v, contain_all = all(contain_v), plot_cs = pl))

  }else if(!return_contain_only && !is.null(true_mean)){
    return(list(levels = levels, U_in = in_list, U_out = out_list,
                contain_individual = contain_v, contain_all = all(contain_v), plot_cs = pl))
  }else if(!return_contain_only && is.null(true_mean)){
    return(list(levels = levels, U_in = in_list, U_out = out_list, plot_cs = pl))
  }else{
    return(list(levels = levels, plot_cs = pl))
  }
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
# remember to explain the 4 choose of types, especially "two-sided"
