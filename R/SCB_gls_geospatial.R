#' Construct Simultaneous Confidence Bands for a Spatial Generalized Least Squares Model
#'
#' @param sp_list A list containing the spatial coordinates and the observations.
#'  Should include the following components:
#'   \itemize{
#'     \item \code{x}: A numeric vector of x-coordinates (e.g., longitude).
#'     \item \code{y}: A numeric vector of y-coordinates (e.g., latitude).
#'     \item \code{obs}: A 3D array of observations with dimensions
#'     \code{length(x)} × \code{length(y)} × \code{n}.
#'   }
#' @param level A optional numeric threshold value used to test whether the
#' estimated mean surface significantly deviates from it. Default is NULL.
#' @param data_fit A design matrix used to fit the generalized least squares
#' (GLS) model. Each row corresponds to an observation, and each column to a covariate
#' to be included in the model. Outcome/observation should not be included.
#' The first column is typically an intercept column,
#' which will contain only 1s, if an intercept is included in the model.
#' Categorical variables in `data_fit` should be converted to dummy variables.
#' Default is `matrix(1, n, 1)` (only keep the intercept term)
#' @param w A numeric vector specifying the target function for constructing the SCB,
#' by giving a linear combination of the regression coefficients in the GLS model.
#' Default is `matrix(1, 1, 1)`, will only construct the SCB for the first regression
#' coefficient.
#' @param correlation A character string specifying the name of
#' the correlation structure (e.g., \code{"corAR1"}, \code{"corCompSymm"})
#' to be used in the GLS model. If \code{NULL}, no correlation structure is assumed.
#' @param corpar A list containing parameters to be passed to
#' the correlation structure function specified in \code{correlation}.
#' @param groups A vector of group identifiers used to define
#' the within-group correlation structure (e.g., repeated measures, time blocks).
#' If not specified, defaults to \code{rep(1, n)}, assuming all observations
#' belong to a single group.
#' @param V An optional array of known covariance matrices of
#' shape \code{[length(x), length(y), n, n]}, where each \code{V[i,j,,]}
#' corresponds to the covariance matrix for the observations at spatial location
#' \code{(x[i], y[j])}. If V is given, then the GLS model will be fitted based on V.
#' Otherwise, the GLS model will be fitted based on correlation structure.
#' If neither is provided, the model reduces to ordinary least squares (OLS) regression.
#' @param alpha A numerical value specifying the confidence level
#' for the Simultaneous Confidence Bands. Defalut is `0.1`.
#' @param N An integer specifying the number of bootstrap samples to
#' construct the Simultaneous Confidence Bands. Default is `1000`.
#' @param mask An optional logical matrix same dimensions as
#' \code{c(length(sp_list$x), length(sp_list$y))}, indicating spatial locations to
#' include in the SCB computation. Non-included locations (e.g., water areas)
#' should be set to 0 or \code{NA}.
#' Default is `array(1, dim = c(length(sp_list$x), length(sp_list$y)))`
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{\code{scb_up}}{A matrix of upper bounds for
#'   the simultaneous confidence bands at each spatial location
#'   corresponding to the target function specified by `w`.}
#'   \item{\code{scb_low}}{A matrix of lower bounds for
#'   the simultaneous confidence bands at each spatial location
#'   corresponding to the target function specified by `w`.}
#'   \item{\code{mu_hat}}{A matrix of estimated mean values at each spatial location
#'   corresponding to the target function specified by `w`.}
#'   \item{\code{norm_est}}{A matrix of standardized test statistics \code{(mu_hat - level) / SE}.}
#'   \item{\code{thres}}{The bootstrap threshold used to construct the confidence bands.}
#'   \item{\code{x}}{The vector of x-coordinates corresponding to the columns of the spatial grid.}
#'   \item{\code{y}}{The vector of y-coordinates corresponding to the rows of the spatial grid.}
#' }
#'
#' @references
#' Sommerfeld, M., Sain, S., & Schwartzman, A. (2018).
#' Confidence regions for spatial excursion sets from repeated random field observations, with an application to climate.
#' \emph{Journal of the American Statistical Association}, 113(523), 1327–1340.
#' \doi{10.1080/01621459.2017.1356318}
#'
#' Ren, J., Telschow, F. J. E., & Schwartzman, A. (2024).
#' Inverse set estimation and inversion of simultaneous confidence intervals.
#' \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, 73(4), 1082–1109.
#' \doi{10.1093/jrsssc/qlae027}
#'
#' @importFrom MASS lm.gls
#' @import nlme
#' @importFrom stats quantile formula sd
#' @importFrom Matrix bdiag
#'
#' @export
#'
#' @examples
#' data(climate_data)
#' # Construct confidence sets for the increase of the mean temperature (June-August)
#' # in North America between the 20th and 21st centuries
#' \dontrun{
#' temp = SCB_gls_geospatial(sp_list = climate_data$Z, level = 2, data_fit = climate_data$X,
#'                        w = c(1,0,0,0), correlation = climate_data$correlation,
#'                        mask = climate_data$mask, alpha = 0.1)
#' }
#'
#' example_list <- list(x = climate_data$Z$x[50:70], y = climate_data$Z$y[40:60],
#' obs = climate_data$Z$obs[50:70, 40:60,])
#' temp = SCB_gls_geospatial(sp_list = example_list, level = 2, data_fit = climate_data$X,
#'                        w = c(1,0,0,0), correlation = NULL,
#'                        mask = NULL, alpha = 0.1, N = 100)
#'
SCB_gls_geospatial =
  function (sp_list, level = NULL, data_fit = NULL, w = NULL,
            correlation = NULL, corpar = NULL, groups = NULL, V = NULL,
            alpha = 0.1, N = 1000, mask = NULL)
  {
    # require(nlme)
    if(is.null(sp_list)) stop("Must provide input for `sp_list`.")
    if(!is.list(sp_list)) stop("`sp_list` should be a list.")
    if(!all(c("x", "y", "obs") %in% names(sp_list))) {
      stop("`sp_list` must have elements named 'x', 'y' and 'obs'.")
    }

    if (!(is.numeric(sp_list$x) && length(dim(sp_list$x)) <= 1)) {
      stop("`sp_list$x` must be a 1D numeric vector.")
    }
    if (!(is.numeric(sp_list$y) && length(dim(sp_list$y)) <= 1)) {
      stop("`sp_list$y` must be a 1D numeric vector.")
    }
    if (!(is.array(sp_list$obs) && is.numeric(sp_list$obs))) {
      stop("`sp_list$obs` must be a numeric array.")
    }
    if (length(dim(sp_list$obs)) != 3) {
      stop("`sp_list$obs` must be a 3D array with dims [length(x), length(y), n].")
    }
    d <- dim(sp_list$obs)
    if (d[1] != length(sp_list$x) || d[2] != length(sp_list$y)) {
      stop(sprintf("`sp_list$obs` dims mismatch: got [%d, %d, %d], expected [%d, %d, n].",
                   d[1], d[2], d[3], length(sp_list$x), length(sp_list$y)))
    }
    if (d[3] < 1) {
      stop("`sp_list$obs` third dimension (n) must be >= 1.")
    }

    if(!is.null(level)) {
      if(!(is.numeric(level) && length(level) == 1)) stop("`level` must be a single numeric value.")
    }

    if (!is.numeric(N) || N <= 0 || N %% 1 != 0){
        stop("`N` must be a positive integer.")
    }

    if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1){
      stop("`alpha` must be in (0, 1).")
    }

    x = sp_list$x
    y = sp_list$y
    Y = sp_list$obs # observations
    n = dim(Y)[3]
    nloc <- length(x) * length(y)

    if (is.null(data_fit)) {
      data_fit <- matrix(1, n, 1) # design matrix
    }else{
      if(!(is.matrix(data_fit)||is.array(data_fit)) || !is.numeric(data_fit)){
        stop("`data_fit` should be a numeric matrix or array.")
      }
      if (nrow(data_fit) != n) {
        stop("The number of rows in `data_fit` must be equal to the third dimension of `sp_list$obs`.")
      }
    }

    if (is.null(w)) {
      w <- matrix(1, 1, 1) # covariate (linear combination)
    }else{
      if(!((is.atomic(w) && is.null(dim(w)))||is.array(w)||is.matrix(w)) || !is.numeric(w)){
        stop("`w` should be a numeric vector, matrix or array.")
      }
    }

    p <- ncol(data_fit)
    if(is.vector(w)){
      if (length(w) != p) {
        stop("The length of `w` must be equal to the number of columns of `data_fit`.")
      }
    }else{
      if (nrow(w) != p && length(dim(w)) != 1) {
        stop("Dimension of `w` should be 1, and the number of rows in `w` must
             be equal to the number of columns of `data_fit`.")
      }
    }

    M <- c(length(x), length(y))
    if (!is.null(mask)) {
      if (!is.array(mask)) {
        stop("`mask` must be a numeric matrix or array.")
      }
      if (length(dim(mask)) != 2) {
        stop("`mask` must be 2-dimensional with dimensions length(x) x length(y).")
      }
      if (!identical(dim(mask), M)) {
        stop("`mask` must have dimensions ",
             paste(M, collapse = " x "), ".")
      }
    }

    if (!is.null(V)) {
      if (!(is.array(V)) || (!is.numeric(V) || length(dim(V)) != 4)) {
        stop("`V` must be a numeric 4-dimensional array.")
      }

      dims <- dim(V)
      if (dims[3] != dims[4] || dims[3] != n) {
        stop("The last two dimensions of `V` must be equal (square matrices),
              and also be equal to the third dimension of `Z$obs`.")
      }
    }

    if(!is.null(corpar)){
      #if (!is.list(corpar)) {
        #stop("`corpar` must be a list of parameters for the correlation structure function.")
      #}
      #if (is.null(names(corpar)) || any(names(corpar) == "")) {
        #warning("Some elements in `corpar` are not named; this may cause errors in `do.call`.")
      #}
      #if (!is.numeric(corpar) || !is.vector(corpar)) {
        #stop("`corpar` should be a numeric vector.")
      #}
    }

    if (is.null(groups)) {
      groups <- rep(1, n)
    }
    if (!( (is.atomic(groups) && is.numeric(groups)) || is.factor(groups) )) {
      stop("`groups` must be either a numeric vector or a factor.")
    }
    if (length(groups) != n) {
      stop("`groups` must have length equal to the number of observations (", n, ")/
           the third dimension of `sp_list$obs`.")
    }

    invsqrtm <- function(A) {
      E <- eigen(A)
      U <- E$vectors
      D <- diag(E$values)
      U %*% diag(1/sqrt(E$values)) %*% t(U)
    }
    deR <- array(0, c(length(x), length(y), n)) # for bootstrap
    vabs <- matrix(0, length(x), length(y))
    norm_est <- matrix(0, length(x), length(y))
    mu_hat <- matrix(0, length(x), length(y))
    if (!is.null(correlation))
      correlation = do.call(get(correlation), c(corpar, form = ~1 |
                                                  groups)) # correlation structure between observations

    for (i in 1:length(x)) {
      for (j in 1:length(y)) {
        ytemp <- Y[i, j, ]
        if (sum(is.na(ytemp)) == length(ytemp)) {
          mu_hat[i, j] = NA
          norm_est[i, j] = NA
        }
        else {
          df <- data.frame(cbind(ytemp = ytemp, data_fit, groups = groups))
          df <- df[order(groups), ]
          groups <- sort(groups)
          fo <- paste(names(df)[1], "~", paste(names(df)[-c(1,
                                                            p + 2)], collapse = " + "), "-1") # include intercept in X
          if (is.null(V) && !is.null(correlation)) {
            model <- gls(formula(fo), data = df,
                               correlation = correlation)
          }
          else if(!is.null(V)) {
            model <- lm.gls(formula(fo), data = df,
                                  W = V[i, j, , ], inverse = TRUE)
          }else{
            #stop("Must provide one of 'correlation' and 'V'.")
            model <- gls(formula(fo), data = df)
          }

          mu_hat[i, j] <- t(w) %*% model$coefficients
          if (!is.null(correlation)) {
            cM <- corMatrix(model$modelStruct$corStruct,
                                  corr = F) # list of covariance matrix for all groups
            if (!is.list(cM))
              cM <- list(cM)
            invsqrtmOmega <- as.matrix(bdiag(cM)) # intergrate
            deR[i, j, ] <- invsqrtmOmega %*% model$residuals # decorrelation
            deR[i, j, ] <- deR[i, j, ]/sd(deR[i, j, ]) # standardization
          }
          else if (!is.null(V)) {
            deR[i, j, ] <- chol(solve(V[i, j, , ])) %*%
              model$residuals
            deR[i, j, ] <- deR[i, j, ]/sd(deR[i, j, ])
            #design_mat <- model.matrix(model)
            #attributes(design_mat) <- attributes(design_mat)[c("dim", "dimnames")]
            model$varBeta = solve(t(data_fit) %*% solve(V[i, j,
                                                   , ], data_fit))
          }
          else {
            deR[i, j, ] <- model$residuals
            deR[i, j, ] <- deR[i, j, ]/sd(deR[i, j, ])
          }
          vabs[i, j] <- sqrt(t(w) %*% model$varBeta %*%
                               w)
          norm_est[i, j] <- (mu_hat[i, j] - level)/vabs[i,
                                                        j]
        }
        #first_iter <- FALSE
      }
    }
    if (is.null(mask)) {
      mask = array(1, dim = c(length(x), length(y)))
    }
    else {
      mask[which(!is.na(mask), arr.ind = TRUE)] = 1
    }
    mu_hat <- mu_hat * mask
    norm_est <- norm_est * mask
    deR_mask = sweep(deR, 1:2, mask, FUN = "*" )
    a_MB = quantile(MB_(x = x, y = y, R = deR, N = N),
                    probs = 1 - alpha, type = 8)
    norm_est[i, j] <- (mu_hat[i, j] - level)/vabs[i,
                                                  j]
    scb_up = mu_hat + a_MB*vabs
    scb_low = mu_hat - a_MB*vabs
    # return index, scb_up, scb_low
    return(list(scb_up = scb_up, scb_low = scb_low, mu_hat = mu_hat,
                norm_est = norm_est, thres = a_MB, x = x, y = y))
  }

#' Multiplier Bootstrap for Simultaneous Confidence Band Threshold
#'
#' Internal function used to compute the threshold value for constructing simultaneous confidence bands via multiplier bootstrap.
#'
#' @param x A numeric vector of x-coordinates.
#' @param y A numeric vector of y-coordinates.
#' @param R A 3D array of standardized residuals with dimensions \code{[length(x), length(y), n]}, where \code{n} is the sample size.
#' @param N An integer specifying the number of bootstrap samples. Default is \code{1000}.
#'
#' @return A numeric vector of length \code{N}, containing the maximum standardized deviation across all spatial locations for each bootstrap sample.
#'   These can be used to compute the \code{(1 - alpha)} quantile as the SCB threshold.
#'
#' @importFrom stats rnorm
#' @keywords internal
#'
#' @examples
#' # Used internally by SCB_gls_geospatial
#'
MB_ = function (x, y, R, N = 1000)
{
  n = dim(R)[3]
  g = matrix(rnorm(n * N), n, N)
  apply(abs(matrix(R, ncol = n) %*% g), 2, max, na.rm = T)/sqrt(n -
                                                                  2)
}
