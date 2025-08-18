#' Construct Simultaneous Confidence Bands for a Spatial Generalized Least Square Model
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
#' @param data_fit A named matrix or data frame used to fit the generalized least squares
#' (GLS) model. Each row corresponds to an observation, and each column to a covariate
#' to be included in the model. Outcome/boservation should not be included.
#' The first column is typically an intercept column,
#' which will contain only 1s, if an intercept is included in the model.
#' Default is `matrix(1, n, 1)` (only keep the intercept term)
#' @param subset An atomic character vector (e.g., c("X1 = 1"))
#' specified the target function for constructing the SCB.
#' Each element must be of the form <name> = <value>, where <name> is the name
#' of a covariate in data_fit and <value> is the desired value.
#' Whitespace is ignored. Default is NULL, will only construct the SCB for the
#' first covariate in `data_fit`.
#' @param correlation A character string specifying the name of
#' the correlation structure (e.g., \code{"corAR1"}, \code{"corCompSymm"})
#' to be used in the GLS model. If \code{NULL}, no correlation structure is assumed.
#' @param corpar A list of parameters to be passed to
#' the correlation structure function specified in \code{correlation}.
#' @param groups A vector of group identifiers used to define
#' the within-group correlation structure (e.g., repeated measures, time blocks).
#'   If not specified, defaults to \code{rep(1, n)},
#'   assuming all observations belong to a single group.
#' @param V An optional array of known covariance matrices of
#' shape \code{[length(x), length(y), n, n]}, where each \code{V[i,j,,]}
#' corresponds to the covariance matrix for the observations at spatial location
#' \code{(x[i], y[j])}. If V is given, then the GLS model will be fitted based on V.
#' Otherwise, the GLS model will be fitted based on correlation structure.
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
#' @importFrom MASS lm.gls
#' @importFrom nlme gls corMatrix
#' @importFrom stats quantile formula sd
#' @importFrom Matrix bdiag
#'
#' @export
#'
#' @examples
#' data(climate_data)
#' library(nlme)
#' # Construct confidence sets for the increase of the mean temperature (June-August)
#' # in North America between the 20th and 21st centuries
#' temp = SCB_gls_climate(sp_list = climate_data$Z, level = 2, data_fit = climate_data$X,
#'                        subset = c("X1 = 1"), correlation = climate_data$correlation,
#'                        mask = climate_data$mask, alpha = 0.1)
SCB_gls_climate =
  function (sp_list, level = NULL, data_fit = NULL, subset = NULL,
            correlation = NULL, corpar = NULL, groups = NULL, V = NULL,
            alpha = 0.1, N = 1000, mask = NULL)
  {
    # require(nlme)
    if(!is.list(sp_list)) stop("`sp_list` should be a list.")
    if(!all(c("x", "y", "obs") %in% names(sp_list))) {
      stop("`sp_list` must have elements named 'x', 'y' and 'obs'.")
    }
    if(!is.numeric(sp_list$x)|| !is.numeric(sp_list$y)||!is.numeric(sp_list$obs)){
      stop("All elements in `sp_list` should be numeric.")
    }

    if(!is.null(level)) {
      if(!(is.numeric(level) && length(level) == 1)) stop("`level` must be a single numeric value.")
    }

    if (!is.numeric(N) || N <= 0 || N %% 1 != 0){
        stop("`nboot` must be a positive integer.")
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
      w <- matrix(1, 1, 1) # covariate (linear combination)
    }else{
      if (!((is.data.frame(data_fit) && !is.null(colnames(data_fit))) ||
            (is.matrix(data_fit) && !is.null(colnames(data_fit))))) {
        stop("`data_fit` should be a named data frame or a matrix with colnames.")
      }

      # transform character variable to factor
      char_vars <- names(data_fit)[sapply(data_fit, is.character)]
      for (v in char_vars) {
        data_fit[[v]] <- factor(data_fit[[v]])
      }

      #if(!(is.vector(w)||is.array(w)||is.matrix(w)) || !is.numeric(w)){
        #stop("`w` should be a numeric vector, matrix or array with only one dimension.")
      #}
      if (nrow(data_fit) != n) {
        stop("The number of rows in `data_fit` must be equal to
             the third dimension of `sp_list$obs`.")
      }
    }

    p <- ncol(data_fit)
    #if(is.vector(w)){
      #if (length(w) != p) {
        #stop("The length of `w` must be equal to the number of rows of `data_fit`.")
      #}
    #}else{
      #if (nrow(w) != p && length(dim(w)) != 1) {
        #stop("Dimension of `w` should be 1, and the number of rows in `w` must
             #be equal to the number of rows of `data_fit`.")
      #}
    #}

    M <- c(length(x), length(y))
    if (!is.null(mask)) {
      if(!identical(dim(mask), M)) stop("`mask` must have dimensions ",
                                        paste(M, collapse = " x "), ".")
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
              and also be equal to the third dimension of `Z$obs`.")
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
    if (!( (is.vector(groups) && is.numeric(groups)) || is.factor(groups) )) {
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

    first_iter <- TRUE
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
            stop("Must provide one of 'correlation' and 'V'.")
          }
          if (!is.null(subset) && first_iter){
            if(is.null(data_fit)){
              stop("Must provide input for `data_fit` if `subset` is not NULL.")
            }
            w <- rep(0, length(model$coefficients))
            # Identify covariates
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

                if(is.data.frame(data_fit)){
                  if (!var %in% names(data_fit)) {
                    stop(paste0("Variable '", var, "' not found in `data_fit`."))
                  }
                  if (is.factor(df[[var]])) {
                    x <- data_fit[[var]]
                    lv <- levels(x)
                    if (!(val %in% lv)) {
                      stop(sprintf("Value '%s' not in levels of variable `%s`", val, var))
                    }
                    code <- which(lv == val)
                    col_idx <- match(var, colnames(data_fit))
                    w[col_idx+code-1] <- 1
                  }

                  if (is.numeric(data_fit[[var]])) { # numeric
                    #unique_vals <- unique(data_df[[var]])
                    #if (!(val %in% unique_vals)) {
                    #stop(paste0("Value '", val, "' not found in numeric variable '", var, "'. "))
                    #}
                    col_idx <- match(var, colnames(data_fit))
                    w[col_idx] <- val
                  }else{
                    stop(paste0("The variable '", var, "' is not of type numeric or factor. ",
                                "Please convert it to a numeric/factor variable."))
                  }
                }else{
                  if (!var %in% colnames(data_fit)) {
                    stop(paste0("Variable '", var, "' not found in `data_fit`."))
                  }
                  col_inx <- which(colnames(data_fit) == var)
                  if (is.factor(data_fit[,col_inx])) {
                    x <- data_fit[,col_inx]
                    lv <- levels(x)
                    if (!(val %in% lv)) {
                      stop(sprintf("Value '%s' not in levels of variable `%s`", val, var))
                    }
                    code <- which(lv == val)
                    w[col_inx+code-1] <- 1
                  }

                  if (is.numeric(data_fit[,col_inx])) { # numeric
                    #unique_vals <- unique(data_df[[var]])
                    #if (!(val %in% unique_vals)) {
                    #stop(paste0("Value '", val, "' not found in numeric variable '", var, "'. "))
                    #}
                    #col_idx <- match(var, colnames(data_fit))
                    w[col_inx] <- val
                  }else{
                    stop(paste0("The variable '", var, "' is not of type numeric or factor. ",
                                "Please convert it to a numeric/factor variable."))
                  }
                }
              }
            }
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
            design_mat <- model.matrix(model)
            attributes(design_mat) <- attributes(design_mat)[c("dim", "dimnames")]
            model$varBeta = solve(t(design_mat) %*% solve(V[i, j,
                                                   , ], design_mat))
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
        first_iter <- FALSE
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
#' # Used internally by SCB_gls_climate
#'
MB_ = function (x, y, R, N = 1000)
{
  n = dim(R)[3]
  g = matrix(rnorm(n * N), n, N)
  apply(abs(matrix(R, ncol = n) %*% g), 2, max, na.rm = T)/sqrt(n -
                                                                  2)
}
