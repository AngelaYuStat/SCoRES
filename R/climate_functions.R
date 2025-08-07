#' Construct Simultaneous Confidence Bands for a Generalized Least Square Model
#'
#' @param Z A list containing the spatial coordinates and the observations.
#'  Should include the following components:
#'   \itemize{
#'     \item \code{x}: A numeric vector of x-coordinates (e.g., longitude).
#'     \item \code{y}: A numeric vector of y-coordinates (e.g., latitude).
#'     \item \code{z}: A 3D array of observations with dimensions \code{length(x)} × \code{length(y)} × \code{n}.
#'   }
#' @param level A optional numeric threshold value used to test whether the estimated mean surface significantly deviates from it. Default is NULL.
#' @param X A design matrix used in the generalized least squares (GLS) model. Each row corresponds to an observation, and each column to a covariate.
#'          Default is `matrix(1, n, 1)` (only keep the intercept term)
#' @param w A numeric vector specifying a linear combination of the regression coefficients. Default is `matrix(1, 1, 1)`.
#' @param correlation A character string specifying the name of the correlation structure (e.g., \code{"corAR1"}, \code{"corCompSymm"})
#'   to be used in the GLS model. If \code{NULL}, no correlation structure is assumed.
#' @param corpar A list of parameters to be passed to the correlation structure function specified in \code{correlation}.
#' @param groups A vector of group identifiers used to define the within-group correlation structure (e.g., repeated measures, time blocks).
#'   If not specified, defaults to \code{rep(1, n)}, assuming all observations belong to a single group.
#' @param V An optional array of known covariance matrices of shape \code{[length(x), length(y), n, n]},
#'   where each \code{V[i,j,,]} corresponds to the covariance matrix for the observations at spatial location \code{(x[i], y[j])}.
#'   If V is given, then the GLS model will be fitted based on V. Otherwise, the GLS model will be fitted based on correlation structure.
#' @param alpha A numerical value specifying the confidence level for the Simultaneous Confidence Bands. Defalut is `0.1`.
#' @param N An integer specifying the number of bootstrap samples to construct the Simultaneous Confidence Bands. Default is `1000`.
#' @param mask An optional logical matrix same dimensions as \code{c(length(Z$x), length(Z$y))}, indicating spatial locations to include in the SCB computation.
#'   Non-included locations (e.g., water areas) should be set to 0 or \code{NA}. Default is `array(1, dim = c(length(Z$x), length(Z$y)))`
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
  function (Z, level = NULL, X = NULL, w = NULL, correlation = NULL, corpar = NULL,
            groups = NULL, V = NULL, alpha = 0.1, N = 1000,
            mask = NULL)
  {
    # require(nlme)
    if(!is.list(Z)) stop("`Z` should be a list.")
    if(!all(c("x", "y", "z") %in% names(Z))) {
      stop("`Z` must have elements named 'x', 'y' and 'z'.")
    }
    if(!is.numeric(Z$x)|| !is.numeric(Z$y)||!is.numeric(Z$z)) stop("All elements in `Z` should be numeric.")

    if(!is.null(level)) {
      if(!(is.numeric(level) && length(level) == 1)) stop("`level` must be a single numeric value.")
    }

    if (!is.numeric(N) || N <= 0 || N %% 1 != 0){
        stop("`nboot` must be a positive integer.")
    }

    if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1){
      stop("`alpha` must be in (0, 1).")
    }

    x = Z$x
    y = Z$y
    Y = Z$z # observations
    n = dim(Y)[3]
    nloc <- length(x) * length(y)
    if (is.null(X) || is.null(w)) {
      X <- matrix(1, n, 1) # design matrix
      w <- matrix(1, 1, 1) # covariate (linear combination)
    }else{
      if(!(is.matrix(X)||is.array(X)) || !is.numeric(X)){
        stop("`X` should be a numeric matrix or array.")
      }
      if(!(is.vector(w)||is.array(w)||is.matrix(w)) || !is.numeric(w)){
        stop("`w` should be a numeric vector, matrix or array.")
      }
      if (nrow(X) != n) {
        stop("The number of rows in `X` must be equal to the third dimension of `Z$z`.")
      }
    }

    p <- ncol(X)
    if(is.vector(w)){
      if (length(w) != p) {
        stop("The length of `w` must be equal to the number of rows of `X`.")
      }
    }else{
      if (nrow(w) != p && length(dim(w)) != 1) {
        stop("Dimension of `w` should be 1, and the number of rows in `w` must
             be equal to the number of rows of `X`.")
      }
    }

    M <- c(length(x), length(y))
    if (!is.null(mask)) {
      if(!identical(dim(mask), M)) stop("`mask` must have dimensions ", paste(M, collapse = " x "), ".")
      if (!(is.matrix(mask) || is.array(mask))) {
        stop("`mask` must be a matrix or array.")
      }
    }

    if (!is.null(V)) {
      if (!(is.matrix(V) || is.array(V)) || (!is.numeric(V) || length(dim(V)) != 4)) {
        stop("`V` must be a numeric 4-dimensional matrix or array.")
      }

      dims <- dim(V)
      if (dims[3] != dims[4] || dims[3] != n) {
        stop("The last two dimensions of `V` must be equal (square matrices),
              and also be equal to the third dimension of `Z$z`.")
      }
    }

    if(!is.null(corpar)){
      if (!is.list(corpar)) {
        stop("`corpar` must be a list of parameters for the correlation structure function.")
      }
      if (is.null(names(corpar)) || any(names(corpar) == "")) {
        warning("Some elements in `corpar` are not named; this may cause errors in `do.call`.")
      }
    }

    if (is.null(groups)) {
      groups <- rep(1, n)
    }
    if (!is.vector(groups) && !is.factor(groups)) {
      stop("`groups` must be a vector or factor.")
    }
    if (length(groups) != n) {
      stop("`groups` must have length equal to the number of observations (", n, ")/ the third dimension of `Z$z`.")
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

# Make Data Object
