#' Plot Inversion of Simultaneous Confidence Bands (SCBs) into Inner and Outer Simultaneous Confidence Regions (SCRs)
#'
#' Visualizes simultaneous confidence regions of upper and lower excursion sets for
#' discrete, 1D or 2D data, using contour or band plots.
#' Supports plotting confidence regions at multiple levels and labeling contours.
#'
#' @param SCB A numeric list returned by `regression_outcome_scb()`,
#' `functional_outcome_scb()` or a custom list with two arrays of the same dimension:
#' `scb_up` and `scb_low`, representing the upper and lower confidence bounds
#' respectively. \verb{SCB$scb_up} and \verb{SCB$scb_low} should be numeric vectors (1D)
#' or matrices (2D) containing the upper simultaneous confidence interval.
#' Dimensions of `SCB$scb_up` and `SCB$scb_low` must match.
#' @param levels A numeric vector or list of scalers for different levels or matrix
#' containing interval sets to construct the confidence regions.
#' If \code{type} = "upper" or "lower", `levels` should be a vector.
#' "upper" represents upper excursion sets, and "lower" represents lower excursion sets.
#' @param type A character specifying the type of inverse sets to fit.
#' Choices are `"upper"` and `"lower"`. Default is `"upper"`.
#' @param x A numerical vector of x-axis coordinates for 1D and 2D cases.
#' For discrete coordinates, use a character vector.
#' The order of x should correspond to the order of `scb_up` and `scb_low` in `SCB`.
#' @param y Optional vector of y-axis coordinates for 2D data.
#' @param mu_hat A numeric array (1D) or matrix (2D) of estimated means.
#' If `mu_true` is provided, this will be overwritten by the true mean. Default is NULL.
#' An input must be provided for either `mu_hat` or `mu_true`.
#' @param mu_true Optional numeric array (1D) or matrix (2D) of true means,
#' which overrides `mu_hat` if provided. Default is NULL.
#' @param together Optional logical value for plotting option.
#' If `TRUE`, plots all confidence levels on the same figure;
#' otherwise, generates one plot per level. Default is `TRUE`.
#' @param xlab Optional character for the label of the x-axis. Default is `"x1"`.
#' @param ylab Optional character for the label of the y-axis. Default is `"x2"`.
#' @param level_label Optional logical input for level displaying option.
#' If `TRUE`, displays numeric level labels on contour lines for 2D confidence sets.
#' Default is `TRUE`.
#' @param min.size Optional logical input for minimum number of points
#' required for a contour to be labeled. Default is `5`.
#' @param palette Optional character value for the name of the HCL color palette
#' to use when plotting multiple levels together. Default is `"gray"`.
#' @param color_level_label Optional character value for the color used
#' for contour level labels. Default is `"black"`.
#'
#' @returns A \code{ggplot2} object that includes both simultaneous confidence intervals
#' and simultaneous confidence region of excursion sets corresponding to levels assigned.
#'
#' @importFrom metR geom_text_contour
#' @import ggplot2
#' @importFrom patchwork plot_layout wrap_plots
#' @importFrom reshape melt
#' @importFrom dplyr mutate desc
#' @importFrom magrittr %>%
#' @importFrom forcats fct_reorder
#' @importFrom grDevices hcl.colors
#'
#' @export
#'
#' @references
#' Ren, J., Telschow, F. J. E., & Schwartzman, A. (2024).
#' Inverse set estimation and inversion of simultaneous confidence intervals.
#' \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, 73(4), 1082–1109.
#' \doi{10.1093/jrsssc/qlae027}
#'
#' @examples
#'
#' # example using pupil data
#' library(mgcv)
#' data(pupil)
#' \dontrun{
#' pupil_fpca <- prepare_pupil_fpca(pupil)
#'
#' fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
#'   s(seconds, by = use, k=30, bs = "cr") +
#'   s(id, by = Phi1, bs="re") +
#'   s(id, by = Phi2, bs="re") +
#'   s(id, by = Phi3, bs="re") +
#'   s(id, by = Phi4, bs="re"),
#'   method = "fREML", data = pupil_fpca, discrete = TRUE)
#'
#' pupil_multiplier <- SCB_functional_outcome(data = pupil_fpca, object = fosr_mod,
#'                                    method = "multiplier",
#'                                    outcome = "percent_change",
#'                                    domain = "seconds", subset= c("use = 1"),
#'                                    id = "id")
#'
#' pupil_multiplier <- tibble::as_tibble(pupil_multiplier)
#'
#' plot_cs(pupil_multiplier,levels = c(-18), x = pupil_multiplier$domain,
#'         mu_hat = pupil_multiplier$mu_hat, xlab = "", ylab = "",
#'         level_label = T, min.size = 40, palette = "Spectral",
#'         color_level_label = "black")
#' }
#'
#' mean_mod <- mgcv::gam(percent_change ~ s(seconds, k = 5, bs = "cr") +
#' s(seconds, by = use, k = 5, bs = "cr"),
#' data = pupil, method = "REML")
#'
#' pupil_multiplier <- SCB_functional_outcome(data = pupil, object = mean_mod,
#'                                    method = "multiplier",
#'                                    outcome = "percent_change",
#'                                    domain = "seconds", subset= c("use = 1"),
#'                                    id = "id", nboot = 50)
#'
#' pupil_multiplier <- tibble::as_tibble(pupil_multiplier)
#'
#' plot_cs(pupil_multiplier,levels = c(-18), x = pupil_multiplier$domain,
#'         mu_hat = pupil_multiplier$mu_hat, xlab = "", ylab = "",
#'         level_label = T, min.size = 40, palette = "Spectral",
#'         color_level_label = "black")

plot_cs = function(SCB, levels, type = "upper", x, y = NULL, mu_hat = NULL, mu_true = NULL, together = TRUE, xlab = "X1", ylab = "X2", level_label = TRUE,
                   min.size = 5, palette = "gray", color_level_label = "black"){

  if(!is.list(SCB)) stop("`SCB` should be a list.")
  if(!all(c("scb_low", "scb_up") %in% names(SCB))) {
    stop("`SCB` must have elements named 'scb_low' and 'scb_up'.")
  }else{
    if(!is.numeric(SCB$scb_up) || !is.numeric(SCB$scb_low)) {
      stop("Values of `SCB$scb_up` and `SCB$scb_low` must be numeric.")
    }
  }

  if(is.null(levels)) {
    stop("Must provide input for `levels`.")
  }else{
    if(type %in% c("upper", "lower")){
      if (!is.null(dim(levels)) && length(dim(levels)) == 1) {
        levels <- as.vector(levels)
      }
      if(!(is.atomic(levels) && is.null(dim(levels)))) stop("`levels` should be a vector if `type` = upper or lower.")
      if(!is.numeric(levels)) {
        stop("Values of `levels` must be numeric.")
      }
    }else if(type == "interval"){
      #if(!is.list(levels)) stop("`levels` should be a list if `type` = interval")
      #if (!all(c("low", "up") %in% names(levels))) {
        #stop("`levels` must have elements named 'low' and 'up'.")
      #}
      #if(!is.numeric(levels$low)||!is.numeric(levels$up)) {
        #stop("All elements in `levels` must be numeric.")
      #}
      stop("'interval' is not avaliable for plotting, please choose between 'upper' and 'lower'.")
    }else if(type == "two-sided"){
      #stop("'two-sided' is not avaliable for plotting, please choose between 'upper', 'lower' or 'interval'.")
      stop("'two-sided' is not avaliable for plotting, please choose between 'upper' and 'lower'.")
    }else{
      #stop("`type` must be chosen between 'upper', 'lower' or 'interval'.")
      stop("`type` must be chosen between 'upper' and 'lower'.")
    }
    levels <- sort(unique(levels))
  }

  if(is.null(x)) {
    stop("Must provide input for `x`.")
  }#else{
    #if(!(is.numeric(x)||is.character(x)) && !(is.vector(x) || is.array(x))) stop("`x` must be a numeric/character vector or numeric/character array.")
  #}

  #if(!is.null(y)) {
    #if(!is.numeric(y) && !(is.vector(y) || is.array(y))) stop("`y` should be a numeric vector or numeric array..")
  #}

  nd <- length(dim(SCB$scb_up))
  if (nd == 0L) {  # 1D
    if(!(length(SCB$scb_up) == length(SCB$scb_low))) {
      stop("Dimensions of `SCB$scb_up` and `SCB$scb_low` must match.")
    }
    if (any(SCB$scb_low > SCB$scb_up, na.rm = TRUE)) {
      warning("Found entries where `scb_low > scb_up`. Please check your inputs.")
    }
    # x must be provided as vector
    if (!is.null(dim(x)) && length(dim(x)) == 1) {
      x <- as.vector(x)
    }
    if (!(is.atomic(x) && is.null(dim(x)))) stop("For 1D, `x` must be a vector.")
    if(!(is.numeric(x)||is.character(x))){
      stop("`x` must be numeric/character.")
    }
    if (length(x) != length(SCB$scb_up)) {
      stop("For 1D, `length(x)` must match length of `SCB$scb_up/scb_low`.")
    }
    if(is.numeric(x)){
      x <- sort(x)
    }

    if(!is.null(mu_hat)){
      if(!is.numeric(mu_hat)) stop("Input values of `mu_hat` must be numeric.")
      if(!(length(SCB$scb_up) == length(mu_hat))) {
        stop("Dimensions of `SCB$scb_up`, `SCB$scb_low` and `mu_hat` must match.")
      }
      if (!((is.atomic(mu_hat) && is.null(dim(mu_hat))) || is.matrix(mu_hat) || is.array(mu_hat))) {
        stop("`mu_hat` must be a vector, array or matrix.")
      }
    }
    if(!is.null(mu_true)){
      if(!is.numeric(mu_true)) stop("Input values of `mu_true` must be numeric.")
      if(!(length(SCB$scb_up) == length(mu_true))){
        stop("Dimensions of `SCB$scb_up`, `SCB$scb_low` and `mu_true` must match.")
      }
      if (!((is.atomic(mu_true) && is.null(dim(mu_true))) || is.matrix(mu_true) || is.array(mu_true))) {
        stop("`mu_true` must be a vector, array or matrix.")
      }
    }
  } else if (nd == 2L) {  # 2D
    if(!(identical(dim(SCB$scb_up), dim(SCB$scb_low)))) {
      stop("Dimensions of `SCB$scb_up` and `SCB$scb_low` must match.")
    }
    if (any(SCB$scb_low > SCB$scb_up, na.rm = TRUE)) {
      warning("Found entries where `scb_low > scb_up`. Please check your inputs.")
    }
    if(!(is.numeric(x))){
      stop("`x` must be numeric.")
    }
    if (!is.null(dim(x)) && length(dim(x)) == 1) {
      x <- as.vector(x)
    }
    if (!(is.atomic(x) && is.null(dim(x)))) stop("For 2D, `x` must be a vector.")
    if(!(is.numeric(y))){
      stop("`y` must be numeric.")
    }
    if (!is.null(dim(y)) && length(dim(y)) == 1) {
      y <- as.vector(y)
    }
    if (!(is.atomic(y) && is.null(dim(y)))) stop("For 2D, `y` must be a vector.")
    if (length(x) != nrow(SCB$scb_up)) {
      stop("For 2D, `length(x)` must equal `nrow(SCB$scb_up)`.")
    }
    x <- sort(x)
    if (!is.null(y) && length(y) != ncol(SCB$scb_up)) {
      stop("For 2D, `length(y)` must equal `ncol(SCB$scb_up)`.")
    }
    y <- sort(y)
    if(!is.null(mu_hat)){
      if(!is.numeric(mu_hat)) stop("Input values of `mu_hat` must be numeric.")
      if(!identical(dim(SCB$scb_up), dim(mu_hat))) {
        stop("Dimensions of `SCB$scb_up`, `SCB$scb_low` and `mu_hat` must match.")
      }
      if (!(is.atomic(mu_hat) && is.null(dim(mu_hat)) || is.matrix(mu_hat) || is.array(mu_hat))) {
        stop("`mu_hat` must be a vector, array or matrix.")
      }
    }
    if(!is.null(mu_true)){
      if(!is.numeric(mu_true)) stop("Input values of `mu_true` must be numeric.")
      if(!identical(dim(SCB$scb_up), dim(mu_true))) {
        stop("Dimensions of `SCB$scb_up`, `SCB$scb_low` and `mu_true` must match.")
      }
      if (!(is.atomic(mu_true) && is.null(dim(mu_true)) || is.matrix(mu_true) || is.array(mu_true))) {
        stop("`mu_true` must be a vector, array or matrix.")
      }
    }
  }

  if(is.null(mu_hat) && is.null(mu_true)) stop("An input must be provided for either `mu_hat` or `mu_true`.")

  if(!is.null(xlab)){
    if (!is.character(xlab) || length(xlab) != 1L) stop("`xlab` must be a single string.")
  }
  if (!is.null(ylab)){
    if ((!is.character(ylab) || length(ylab) != 1L)) {
      stop("`ylab` must be NULL or a single string.")
    }
  }

  if (!is.character(palette) || length(palette) != 1) {
    stop("`palette` must be a single character string.")
  }

  if (!is.character(color_level_label) || length(color_level_label) != 1) {
    stop("`color_level_label` must be a single character string.")
  }

  dim = length(dim(SCB$scb_up))
  if(dim == 0){# 1 dimension case
    if(is.character(x)){# if we have discrete coordinates
      if(is.null(mu_true)){
        df_plot = data.frame(x = x, low = SCB$scb_low, up = SCB$scb_up,est_mean = mu_hat) %>%
          mutate(x = fct_reorder(x, desc(low)))
        p = df_plot %>% ggplot(aes(x = x))+ geom_errorbar(aes(ymin = low, ymax = up)) +
          geom_hline(yintercept = levels, linetype="dashed")
        p = p + geom_point(data = df_plot,aes(x = x, y = mu_hat),color = "black")
        #if(type == "interval"){ i <- 1 }
        for(l in levels){
          if(type == "upper"){
            df_plot_l = df_plot %>%
              mutate(l_in = ifelse(low >= l, l, NA),
                     l_est = ifelse(est_mean >= l & is.na(l_in), l, NA),
                     l_out = ifelse(up >= l & is.na(l_est) & is.na(l_in), l, NA)
              )
          }else if(type == "lower"){
            df_plot_l = df_plot %>%
              mutate(l_in = ifelse(up <= l, l, NA),
                     l_est = ifelse(est_mean <= l & is.na(l_in), l, NA),
                     l_out = ifelse(low <= l & is.na(l_est) & is.na(l_in), l, NA)
              )
          }else if(type == "interval"){
            # l_dim = dim(levels)
            #l = levels[i,]
            #df_plot_l = df_plot %>%
              #mutate(l_in = ifelse(low >= l$low & up <= l$up, l, NA),
                     #l_est = ifelse((est_mean >= l$low & est_mean <= l$up) & is.na(l_in), l$low, NA),
                     #l_out = ifelse((up >= l$low & low <= up) & is.na(l_est) & is.na(l_in), l$low, NA)
              #)
            #i = i + 1
            stop("'interval' is not avaliable for plotting, please choose between 'upper' and 'lower'.")
          }else if(type == "two-sided"){
            stop("'two-sided' is not avaliable for plotting, please choose between 'upper' and 'lower'.")
          }else{
            stop("`type` must be chosen between 'upper' and 'lower'.")
          }
          if(!all(is.na(df_plot_l$l_out))){
            p = p+ geom_point(aes(x, l_out), data = df_plot_l, color = "blue")
          }
          if(!all(is.na(df_plot_l$l_est))){
            p = p+ geom_point(aes(x, l_est), data = df_plot_l, color = "orange")
          }
          if(!all(is.na(df_plot_l$l_in))){
            p = p+ geom_point(aes(x, l_in), data = df_plot_l, color = "red")
          }
        }
      }else{
        df_plot = data.frame(x = x, low = SCB$scb_low, true_mean = mu_true, up = SCB$scb_up, est_mean = mu_hat)%>%
          mutate(x = fct_reorder(x, desc(low)))
        p = df_plot %>% ggplot(aes(x = x, y = true_mean))+ geom_point(color = "black")+ geom_errorbar(aes(ymin = low, ymax = up)) +
          geom_hline(yintercept = levels, linetype="dashed") +
          geom_text(data = data.frame(x = rep(levels(df_plot$x)[1], length(levels)), y = levels, labels = levels),
                    aes(x = x, y = y, label = labels,vjust = -0.5, hjust = -0.05))
        #scale_y_continuous(breaks =levels)
        #p = p + geom_point(data = df_plot,aes(x = x, y = mu_hat),color = "black")
        #if(type == "interval"){ i <- 1}
        for(l in levels){
          if(type == "upper"){
            df_plot_l = df_plot %>%
              mutate(l_in = ifelse(low >= l, l, NA),
                     l_est = ifelse(true_mean >= l & is.na(l_in), l, NA),
                     l_out = ifelse(up >= l & is.na(l_est) & is.na(l_in), l, NA)
              )
          }else if(type == "lower"){
            df_plot_l = df_plot %>%
              mutate(l_in = ifelse(up <= l, l, NA),
                     l_est = ifelse(true_mean <= l & is.na(l_in), l, NA),
                     l_out = ifelse(low <= l & is.na(l_est) & is.na(l_in), l, NA)
              )
          }#else if(type == "interval"){
            # l_dim = dim(levels)
            #l = levels[i,]
            #df_plot_l = df_plot %>%
              #mutate(l_in = ifelse(low >= l$low & up <= l$up, l$low, NA),
                     #l_est = ifelse((true_mean >= l$low & true_mean <= l$up) & is.na(l_in), l$low, NA),
                     #l_out = ifelse((up >= l$low & low <= up) & is.na(l_est) & is.na(l_in), l$low, NA)
              #)
            #i = i + 1
          #}
          if(!all(is.na(df_plot_l$l_out))){
            p = p+ geom_point(aes(x, l_out), data = df_plot_l, color = "blue")
          }
          if(!all(is.na(df_plot_l$l_true))){
            p = p+ geom_point(aes(x, l_true), data = df_plot_l, color = "green")
          }
          if(!all(is.na(df_plot_l$l_in))){
            p = p+ geom_point(aes(x, l_in), data = df_plot_l, color = "red")
          }
        }
      }
      p = p + ggtitle("Confidence regions") + theme_light() +
        theme(plot.title = element_text(face = "bold", hjust = 0.5,), plot.margin = unit(c(0.2,0.2,0.2,0.2), "mm"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        labs(x = "Coefficients", y = "Magnitude")
      return(p)
    }else{# Plotting for continuous x coordinate
      if(!is.numeric(x)) stop("`x` should be of type numeric for continuous variable.")
      if(is.null(mu_true)){
        df_plot = data.frame(x = x, low = SCB$scb_low, up = SCB$scb_up, est_mean = mu_hat)
        p = df_plot %>% ggplot(aes(x = x, y = est_mean)) +
          geom_ribbon(aes(ymin = low, ymax = up),alpha = 0.5) +
          geom_hline(yintercept = levels, linetype="dashed") +
          geom_text(data = data.frame(x = rep(min(x), length(levels)), y = levels, labels = levels),
                    aes(x = x, y = y, label = labels,vjust = -0.5))
        p = p + geom_line(data = df_plot,aes(x = x, y = mu_hat),color = "black")
        #if(type == "interval"){ i <- 1}
        for(l in levels){
          if(type == "upper"){
            df_plot_l = df_plot %>% mutate(l_in = ifelse(low >= l, l, NA),
                                           l_est = ifelse(est_mean >= l, l, NA),
                                           l_out = ifelse(up >= l , l, NA))
          }else if(type == "lower"){
            df_plot_l = df_plot %>% mutate(l_in = ifelse(up <= l, l, NA),
                                           l_est = ifelse(est_mean <= l, l, NA),
                                           l_out = ifelse(low <= l , l, NA))
          }else if(type == "interval"){
            # l_dim = dim(levels)
            #l = levels[i,]
            #df_plot_l = df_plot %>% mutate(l_in = ifelse(low >= l$low & up <= l$up, l$low, NA),
                                           #l_est = ifelse(est_mean >= l$low & est_mean <= l$up, l$low, NA),
                                           #l_out = ifelse(up >= l$low & low <= l$lup , l$low, NA))
            #i = i + 1
            stop("'interval' is not avaliable for plotting, please choose between 'upper' and 'lower'.")
          }else if(type == "two-sided"){
            stop("'two-sided' is not avaliable for plotting, please choose between 'upper' and 'lower'.")
          }else{
            stop("Type must be chosen between 'upper' and 'lower'.")
          }
          if(!all(is.na(df_plot_l$l_out))){
            p = p + geom_line(data = df_plot_l,aes(x = x, y = l_out),color = "blue",lwd=1.5)
          }
          if(!all(is.na(df_plot_l$l_est))){
            p = p + geom_line(data = df_plot_l,aes(x = x, y = l_est),color = "orange",lwd=1.5)
          }
          if(!all(is.na(df_plot_l$l_in))){
            p = p + geom_line(data = df_plot_l,aes(x = x, y = l_in),color = "red",lwd=1.5)
          }
        }
      }else{
        df_plot = data.frame(x = x, low = SCB$scb_low, true_mean = mu_true, up = SCB$scb_up,
                             est_mean = mu_hat)
        p = df_plot %>% ggplot(aes(x = x, y = true_mean))+ geom_line(color = "black")+
          geom_ribbon(aes(ymin = low, ymax = up),alpha = 0.3) +
          geom_hline(yintercept = levels, linetype="dashed")+
          geom_text(data = data.frame(x = rep(min(x), length(levels)), y = levels, labels = levels),
                    aes(x = x, y = y, label = labels, vjust = -0.5))
        #p = p + geom_line(data = df_plot,aes(x = x, y = mu_hat),color = "black")
        #if(type == "interval"){ i <- 1}
        for(l in levels){
          if(type == "upper"){
            df_plot_l = df_plot %>% mutate(l_in = ifelse(low >= l, l, NA),
                                           l_est = ifelse(true_mean >= l, l, NA),
                                           l_out = ifelse(up >= l , l, NA))
          }else if(type == "lower"){
            df_plot_l = df_plot %>% mutate(l_in = ifelse(up <= l, l, NA),
                                           l_est = ifelse(true_mean <= l, l, NA),
                                           l_out = ifelse(low <= l , l, NA))
          }#else if(type == "interval"){
            # l_dim = dim(levels)
            #l = levels[i,]
            #df_plot_l = df_plot %>% mutate(l_in = ifelse(low >= l$low & up <= l$up, l$low, NA),
                                           #l_est = ifelse(true_mean >= l$low & true_mean <= l$up, l$low, NA),
                                           #l_out = ifelse(up >= l$low & low <= l$up , l$low, NA))
            #i = i + 1
          #}
          if(!all(is.na(df_plot_l$l_out))){
            p = p + geom_line(data = df_plot_l,aes(x = x, y = l_out),color = "blue",lwd=1.5)
          }
          if(!all(is.na(df_plot_l$l_true))){
            p = p + geom_line(data = df_plot_l,aes(x = x, y = l_true),color = "green",lwd=1.5)
          }
          if(!all(is.na(df_plot_l$l_in))){
            p = p + geom_line(data = df_plot_l,aes(x = x, y = l_in),color = "red",lwd=1.5)
          }


          # p = p + geom_line(data = df_plot_l,aes(x = x, y = l_out),
          #                   arrow=arrow(ends = 'both',type = "closed", length = unit(1, "mm")), color = "blue") +
          #   geom_line(data = df_plot_l,aes(x = x, y = l_in),
          #             arrow=arrow(ends = 'both', type = "closed", length = unit(1, "mm")), color = "red")
        }
      }
      if(is.null(ylab)){
        p = p + ggtitle("Confidence regions") + theme_light() +
          theme(plot.title = element_text(face = "bold", hjust = 0.5,), plot.margin = unit(c(1,0.2,0.2,0.2), "mm"))+
          labs(x = xlab, y = "")
      }else{
        p = p + ggtitle("Confidence regions") + theme_light() +
          theme(plot.title = element_text(face = "bold", hjust = 0.5,), plot.margin = unit(c(1,0.2,0.2,0.2), "mm"))+
          labs(x = xlab, y = ylab)
      }
      return(p)
    }

  }else if(dim == 2){# 2 dimension case
    # remember to explain why 2D don't need to specify the type of set fitted
    # normalized quantity dataframe
    if(is.null(x)||is.null(y)){
      x = seq(0,1,dim(SCB$scb_up)[1])
      y = seq(0,1,dim(SCB$scb_up)[2])
    }
    if(!(is.numeric(x) || is.numeric(y))){
      stop("`x` and `y` should be of type numeric for 2D data.")
    }
    rownames(SCB$scb_up) = x
    colnames(SCB$scb_up) = y
    rownames(SCB$scb_low) = x
    colnames(SCB$scb_low) = y
    SCB$scb_up = suppressWarnings(melt(SCB$scb_up))
    SCB$scb_low = suppressWarnings(melt(SCB$scb_low))
    colnames(SCB$scb_up)=c("X1","X2", "scb_up")
    colnames(SCB$scb_low)=c("X1","X2", "scb_low")
    if(!is.null(mu_hat)){
      rownames(mu_hat) = x
      colnames(mu_hat) = y
      mu_hat = suppressWarnings(melt(mu_hat))
      colnames(mu_hat)=c("X1","X2", "est")
    }
    # Mean frame for plotting
    # TO GET RID OF THE WARNING REFER TO THIS LINK
    # https://stackoverflow.com/questions/50359647/plotting-2-dimensional-function-in-ggplot
    if(!is.null(mu_true)){
      mu = mu_true
      rownames(mu) = x
      colnames(mu) = y
      mu = suppressWarnings(melt(mu))
      colnames(mu)=c("X1","X2", "true_mean")
    }else if(!is.null(mu_hat)){
      mu = mu_hat
      #rownames(mu) = x
      #colnames(mu) = y
      #mu = suppressWarnings(melt(mu))
      colnames(mu)=c("X1","X2", "estimated_mean")
    }else{
      stop("Must provide 'mu_true' or 'mu_hat'.")
    }
    #Plotting]

    if (!is.logical(together) || length(together) != 1L) {
      stop("`together` must be a single logical value.")
    }

    if(together == T){
      colors = hcl.colors(length(levels),palette = palette, rev =T)
      if(!is.null(mu_true)){
        #browser()
        max_mu = round(max(mu$true_mean, na.rm = T)+0.2, digits = 1)
        min_mu = round(min(mu$true_mean, na.rm = T)-0.2, digits = 1)
        p_u <- ggplot()+
          geom_raster(aes(x=X1, y = X2, fill = true_mean), data = mu)+
          ggtitle("Outer Confidence Regions")+
          scale_fill_distiller(palette = "Spectral", direction = -1, limits = c(min_mu,max_mu), na.value = "transparent")+
          labs(fill="True mean", x = xlab, y = ylab)
        #scale_fill_gradientn(colours =
        #                       c("blue","green", "yellow", "orange", "red"))
        p_l <- ggplot()+
          geom_raster(aes(x=X1, y = X2, fill = true_mean),
                      data = mu,show.legend = F)+
          ggtitle("Inner confidence regions")+
          scale_fill_distiller(palette = "Spectral", direction = -1, na.value = "transparent")+
          labs(fill="True mean", x = xlab, y = ylab)
        p_est <- ggplot()+
          geom_raster(aes(x=X1, y = X2, fill = true_mean),
                      data = mu,show.legend = F)+
          ggtitle("Estimated regions")+
          scale_fill_distiller(palette = "Spectral", direction = -1, na.value = "transparent")+
          labs(fill="True mean", x = xlab, y = ylab)
      }else{
        max_mu = round(max(mu$estimated_mean, na.rm = T)+0.2, digits = 1)
        min_mu = round(min(mu$estimated_mean, na.rm = T)-0.2, digits = 1)
        #browser()
        p_u <- ggplot()+
          geom_raster(aes(x=X1, y = X2, fill = estimated_mean), data = mu)+
          ggtitle("Outer confidence regions")+
          scale_fill_distiller(palette = "Spectral", direction = -1, limits = c(min_mu,max_mu), na.value = "transparent")+
          labs(fill="Estimated \n mean", x = xlab, y = ylab)
        p_l <- ggplot()+
          geom_raster(aes(x=X1, y = X2, fill = estimated_mean),
                      data = mu,show.legend = F)+
          ggtitle("Inner confidence regions")+
          scale_fill_distiller(palette = "Spectral", direction = -1, na.value = "transparent")+
          labs(fill="Estimated \n mean", x = xlab, y = ylab)
        p_est <- ggplot()+
          geom_raster(aes(x=X1, y = X2, fill = estimated_mean),
                      data = mu,show.legend = F)+
          ggtitle("Estimated regions")+
          scale_fill_distiller(palette = "Spectral", direction = -1, na.value = "transparent")+
          labs(fill="Estimated \n mean", x = xlab, y = ylab)
      }
      for(i in 1:length(levels)){
        p_u <- p_u +
          stat_contour(aes(X1, X2, z= scb_up),data = SCB$scb_up,
                       breaks=levels[i],colour = colors[i],lwd=0.9)
        if(level_label){
          p_u <- p_u + geom_text_contour(aes(X1, X2, z= scb_up),label = levels[i],
                                         data = SCB$scb_up,
                                         breaks = levels[i],inherit.aes	= T, colour = color_level_label,
                                         position = position_jitter(), min.size = min.size)
        }
      }
      p_u = p_u +
        theme_light() +
        scale_x_continuous(expand = c(0.01,0.01)) + scale_y_continuous(expand = c(0.01,0.01))+
        theme(plot.title = element_text(face = "bold", hjust = 0.5,),
              plot.margin = unit(c(0.2,10,0.2,0.2), "mm"),
              aspect.ratio=1)

      for(i in 1:length(levels)){
        p_l <- p_l +
          stat_contour(aes(X1, X2, z= scb_low),data = SCB$scb_low,
                       breaks=levels[i],colour = colors[i],lwd=0.9)
        if(level_label){
          p_l <- p_l + geom_text_contour(aes(X1, X2, z= scb_low),label = levels[i],
                                         data = SCB$scb_low,
                                         breaks = levels[i], inherit.aes = T, colour = color_level_label,
                                         position = position_jitter(), min.size = min.size)
        }

        #theme(panel.background = element_blank(), plot.title = element_text(face = "bold"))
      }
      p_l = p_l +
        theme_light() +
        scale_x_continuous(expand = c(0.01,0.01)) + scale_y_continuous(expand = c(0.01,0.01))+
        theme(plot.title = element_text(face = "bold", hjust = 0.5,),
              plot.margin = unit(c(0.2,0.2,0.2,0.2), "mm"),
              aspect.ratio=1)

      for(i in 1:length(levels)){
        p_est <- p_est +
          stat_contour(aes(X1, X2, z= est),data = mu_hat,
                       breaks=levels[i],colour = colors[i],lwd=0.9)
        if(level_label){
          p_est <- p_est +geom_text_contour(aes(X1, X2, z= est),label = levels[i],
                                            data = mu_hat,
                                            breaks = levels[i], inherit.aes = T, colour = color_level_label,
                                            position = position_jitter(), min.size = min.size)
        }

        #theme(panel.background = element_blank(), plot.title = element_text(face = "bold"))
      }
      p_est = p_est +
        theme_light() +
        scale_x_continuous(expand = c(0.01,0.01)) + scale_y_continuous(expand = c(0.01,0.01))+
        theme(plot.title = element_text(face = "bold", hjust = 0.5,),
              plot.margin = unit(c(0.2,10,0.2,0.2), "mm"),
              aspect.ratio=1)
      #p = ggarrange(p_u,p_l, nrow = 1, common.legend = T, legend = "bottom")
      p = p_u + p_est  + p_l
      p[[1]] = p[[1]] + theme(
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank() )
      p[[2]] = p[[2]] + theme(
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank() )

      # Remove title from third subplot
      p[[3]] = p[[3]] + theme(
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),

        axis.ticks.x = element_blank(),
        axis.title.x = element_blank() )

      p = p + plot_layout(guides = "collect")
      return(p)
    }else{
      # Plotting the level one by one
      p <- lapply(1:length(levels), function(i){
        if(!is.null(mu_true)){
          max_mu = round(max(mu$true_mean, na.rm = T)+0.2, digits = 1)
          min_mu = round(min(mu$true_mean, na.rm = T)-0.2, digits = 1)
          temp = ggplot()+
            geom_raster(aes(x=X1, y = X2, fill = true_mean), data = mu)+
            labs(fill="True mean", x = xlab, y = ylab)+ theme_light()
        }else{
          max_mu = round(max(mu$estimated_mean, na.rm = T))
          min_mu = round(min(mu$estimated_mean, na.rm = T))
          temp = ggplot()+
            geom_raster(aes(x=X1, y = X2, fill = estimated_mean), data = mu)+
            labs(fill="Estimated \n mean", x = xlab, y = ylab)+ theme_light()
        }
        if(i == 1){
          temp = temp + theme(axis.ticks.x = element_blank(),
                              axis.title.x = element_blank(),
                              plot.title = element_text(face = "bold", hjust = 0.5,),
                              plot.margin = unit(c(0.2,10,0.2,0.2), "mm"),
                              aspect.ratio=1)
        }
        if(i == 2){
          temp = temp + theme(axis.ticks.y = element_blank(),
                              axis.title.y = element_blank(),
                              plot.title = element_text(face = "bold", hjust = 0.5,),
                              plot.margin = unit(c(0.2,10,0.2,0.2), "mm"),
                              aspect.ratio=1)
        }
        if(i == 3){
          temp = temp + theme(axis.ticks.x = element_blank(),
                              axis.title.x = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank(),
                              plot.title = element_text(face = "bold", hjust = 0.5,),
                              plot.margin = unit(c(0.2,0.2,0.2,0.2), "mm"),
                              aspect.ratio=1)
        }
        temp +
          scale_fill_distiller(palette = "Spectral", direction = -1, limits = c(min_mu,max_mu), na.value = "transparent") +
          stat_contour(aes(X1, X2, z= scb_up),data = SCB$scb_up,
                       breaks=levels[i],colour =  "blue",lwd=0.9)+
          stat_contour(aes(X1, X2, z= scb_low),data = SCB$scb_low,
                       breaks=levels[i],colour = "red",lwd=0.9)+
          stat_contour(aes(X1, X2, z= est),data = mu_hat,
                       breaks=levels[i],colour = "green",lwd=0.9)+
          ggtitle(paste("level =",  levels[i]))+
          scale_x_continuous(expand = c(0.01,0.01)) + scale_y_continuous(expand = c(0.01,0.01))
      })
      return(wrap_plots(p)+ plot_layout(guides = "collect"))
    }
  }else{
    stop("The dimension of `SCB$scb_up` and `SCB$scb_low` exceed 2.")
  }
}
# library(forcats)
# library(dplyr)
# Check whether the input should be 2 dimen
