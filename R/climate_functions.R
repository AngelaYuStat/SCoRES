#' Construct Simultaneous Confidence Bands for a Generalized Least Square Model
#'
#' @param Z A list containing the spatial coordinates and the observations.
#'  Should include the following components:
#'   \itemize{
#'     \item \code{x}: A numeric vector of x-coordinates (e.g., longitude).
#'     \item \code{y}: A numeric vector of y-coordinates (e.g., latitude).
#'     \item \code{z}: A 3D array of observations with dimensions \code{length(x)} × \code{length(y)} × \code{n}.
#'   }
#' @param level A numeric threshold value used to test whether the estimated mean surface significantly deviates from it.
#' @param X A design matrix used in the generalized least squares (GLS) model. Each row corresponds to an observation, and each column to a covariate.
#' @param w A numeric vector specifying a linear combination of the regression coefficients.
#' @param correlation A character string specifying the name of the correlation structure (e.g., \code{"corAR1"}, \code{"corCompSymm"})
#'   to be used in the GLS model. If \code{NULL}, no correlation structure is assumed.
#' @param corpar A list of parameters to be passed to the correlation structure function specified in \code{correlation}.
#' @param groups A vector of group identifiers used to define the within-group correlation structure (e.g., repeated measures, time blocks).
#'   If not specified, defaults to \code{rep(1, n)}, assuming all observations belong to a single group.
#' @param V An optional array of known covariance matrices of shape \code{[length(x), length(y), n, n]},
#'   where each \code{V[i,j,,]} corresponds to the covariance matrix for the observations at spatial location \code{(x[i], y[j])}.
#'   If V is given, then the GLS model will be fitted based on V. Otherwise, the GLS model will be fitted based on correlation structure.
#' @param alpha A numerical value specifying the confidence level for the Simultaneous Confidence Bands.
#' @param N An integer specifying the number of bootstrap samples to construct the Simultaneous Confidence Bands.
#' @param mask An optional logical matrix same dimensions as \code{mu_hat}, indicating spatial locations to include in the SCB computation.
#'   Non-included locations (e.g., water areas) should be set to 0 or \code{NA}.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{\code{scb_up}}{A matrix of upper bounds for the simultaneous confidence bands at each spatial location.}
#'   \item{\code{scb_low}}{A matrix of lower bounds for the simultaneous confidence bands at each spatial location.}
#'   \item{\code{mu_hat}}{A matrix of estimated mean values at each spatial location.}
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
#' @importFrom MASS lm.gls
#' @importFrom nlme gls corMatrix
#' @importFrom stats quantile formula sd
#'
#' @export
#'
#' @examples
#' data(climate_data)
#' library(nlme)
#' # construct confidence sets for the increase of the mean temperature (June-August)
#' # in North America between the 20th and 21st centuries
#' temp = SCB_gls_climate(Z = climate_data$Z, level = 2, X = climate_data$X,
#'                        w = c(1,0,0,0), correlation = climate_data$correlation,
#'                        mask = climate_data$mask, alpha = 0.1)
SCB_gls_climate =
  function (Z, level, X = NULL, w = NULL, correlation = NULL, corpar = NULL,
            groups = NULL, V = NULL, alpha = 0.1, N = 1000,
            mask = NULL)
  {
    # require(nlme)
    x = Z$x
    y = Z$y
    Y = Z$z # observations
    n = dim(Y)[3]
    nloc <- length(x) * length(y)
    if (is.null(X)) {
      X <- matrix(1, n, 1) # design matrix
      w <- matrix(1, 1, 1) # covariate (linear combination)
    }
    p <- ncol(X)
    if (is.null(groups)) {
      groups <- rep(1, n)
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
          df <- data.frame(cbind(ytemp = ytemp, X, groups = groups))
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
            stop("Must provide one of 'correlation' and 'V'.")
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
            model$varBeta = solve(t(X) %*% solve(V[i, j,
                                                   , ], X))
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
    return(list(scb_up = scb_up, scb_low = scb_low, mu_hat = mu_hat, norm_est = norm_est, thres = a_MB,
                x = x, y = y))
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
#' # Used internally by SCB_gls_climate
#'
MB_ = function (x, y, R, N = 1000)
{
  n = dim(R)[3]
  g = matrix(rnorm(n * N), n, N)
  apply(abs(matrix(R, ncol = n) %*% g), 2, max, na.rm = T)/sqrt(n -
                                                                  2)
}

# Delete other functions
# Make Data Object
