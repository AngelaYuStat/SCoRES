#' Multiplier Bootstrap
#'
#' @param R Array of shape (..., N), where N is number of repetitions
#' @param Q Optional second sample array for two-sample SCB
#' @param alpha Significance level (default 0.05)
#' @param Mboots Number of bootstrap replications (default 5000)
#' @param method Method for SD estimation: "t" or "regular"
#' @param weights Multiplier type: "rademacher", "gaussian", or "mammen"
#'
#' @return A list with fields: `z` (distribution), `q` (threshold), and `samples`
#'
#' @references
#' Telschow, F. J. E., & Schwartzman, A. (2022).
#' Simultaneous confidence bands for functional data using the Gaussian Kinematic formula.
#' \emph{Journal of Statistical Planning and Inference}, 216, 70–94.
#' \doi{10.1016/j.jspi.2021.05.008}
#'
#' @keywords internal
#'
#' @examples
#' # Used internally by SCB_dense
#'
MultiplierBootstrap <- function( R,
                                 Q       = NULL,
                                 alpha   = 0.05,
                                 Mboots  = 5e3,
                                 method  = "t",
                                 weights = "rademacher"){

  #---- Check user input and put default values
  # Check R
  if( !is.array( R ) ){
    stop("'R' must be an array.")
  }
  # Check Q
  if( !( is.array( Q ) | is.null( Q ) ) ){
    stop("'Q' must be an array.")
  }

  #---- Check params and put defaults, if missing
  # Check Mboots
  if( is.numeric( Mboots ) ){
    if( Mboots %% 1 != 0 & Mboots <= 0 ){
      stop( "The input 'Mboots' needs to be a positiv natural number." )
    }
  }else{
    stop("The input 'Mboots' needs to be a positiv natural number.")
  }

  # Check alpha
  if( is.numeric( alpha ) ){
    if( alpha <= 0 | alpha >= 1 ){
      stop("The input 'alpha' needs to be a real number between 0 and 1.")
    }
  }else{
    stop("The input 'alpha' needs to be a real number between 0 and 1.")
  }

  # Check method
  if( is.character( method ) ){
    if( !( method %in% c( "t", "regular" ) ) ){
      stop( "Please, specify a valid choice for 'method'. Options are 't' and 'regular'.")
    }
  }else{
    stop( "Please, specify a valid choice for 'method'. Options are 't' and 'regular'.")
  }

  # Check weights
  if( is.character( weights ) ){
    if( !( weights %in% c( "gaussian", "rademacher", "mammen" ) ) ){
      stop( "Please, specify a valid choice for 'weights'. Options are 'gaussian', 'rademacher' and 'mammen'.")
    }
  }else{
    stop( "Please, specify a valid choice for 'weights'. Options are 'gaussian', 'rademacher' and 'mammen'.")
  }

  #----- One sample case
  if( is.null( Q ) ){
    #----- Precompute useful constants
    # dimension of input
    dimR = dim( R )

    # number of samples
    N    = dimR[ length( dimR ) ]
    D    = length( dimR ) - 1

    #----- Simulate multiplier weights
    if( weights == "gaussian" ){
      multiplier <- matrix( rnorm( N * Mboots ), N, Mboots )
    }else if( weights == "rademacher" ){
      multiplier <- matrix( sample( c( -1, 1 ), N * Mboots, replace = T ), N, Mboots )
    }else{
      multiplier <- matrix( sqrt( 5 ) *
                              rbinom( N * Mboots,
                                      1,
                                      ( sqrt( 5 ) - 1 ) / 2 / sqrt( 5 ) ) +
                              ( 1 - sqrt( 5 ) ) / 2,
                            N,
                            Mboots )
    }

    #----- Switch by dimension for computational speed
    if( D == 1 ){
      # Compute bootstrap means
      bootMeans <- R %*% multiplier / N

      # Estimate the variance from the sample
      if( method == "regular" ){
        data.sigma <- sqrt( matrixStats::rowVars( R ) )
      }else if( method == "t" ){
        bootSecMoments <- R^2 %*% multiplier^2 / N
        # We put an abs here to make sure that no NaNs are produced due to machine precision error.
        data.sigma <- sqrt( ( N / ( N - 1 ) ) * abs( bootSecMoments - bootMeans^2 ) )
      }

      # Compute bootstrap distribution of the maximum
      distVec <- sqrt( N ) * apply( abs( bootMeans / data.sigma ), 2, max )

    }else if( D > 1 ){
      # Compute bootstrap means
      bootMeans = array( matrix( R, prod( dimR[ -( D + 1 ) ] ), N ) %*% multiplier,
                         dim = c( dimR[ -( D + 1 ) ], Mboots ) ) / N

      # Estimate the variance from the sample
      if( method == "regular" ){
        data.sigma <- sqrt( apply( R, 1:D, var  ) )
        distVec    <- apply( array( as.vector( abs( bootMeans ) ) / as.vector( data.sigma ),
                                    dim = c( dimR[ -( D + 1 ) ], Mboots ) ), D + 1, max)  * sqrt( N )
      }else if( method == "t" ){
        bootSecMoments <- array( matrix( R^2, prod( dimR[ -( D + 1 ) ] ), N ) %*% multiplier^2,
                                 dim = c( dimR[ -( D + 1 ) ], Mboots ) ) / N
        # We put an abs here to make sure that no NaNs are produced due to machine precision error.
        data.sigma <- sqrt( ( N / ( N - 1 ) ) * abs( bootSecMoments - bootMeans^2 ) )
        distVec    <- apply( abs( bootMeans ) / data.sigma, D + 1, max ) * sqrt( N )
      }
    }
  }else{
    #----- Precompute useful constants
    # dimension of input
    dimR = dim( R )
    dimQ = dim( Q )
    # number of samples
    N    = dimR[ length( dimR ) ]
    M    = dimQ[ length( dimQ ) ]
    c    = N / M
    D    = length( dimR ) - 1

    #----- Obtain multiplier weights
    if( weights == "gaussian" ){
      multiplierR <- matrix( rnorm( N * Mboots ), N, Mboots )
      multiplierQ <- matrix( rnorm( M * Mboots ), M, Mboots )
    }else{
      multiplierR <- matrix( sample( 1:2, N * Mboots, replace = T ) * 2 - 3, N, Mboots )
      multiplierQ <- matrix( sample( 1:2, M * Mboots, replace = T ) * 2 - 3, M, Mboots )
    }

    #----- Switch by dimension for computational speed
    if( D == 1 ){
      #----- Compute bootstrap means
      bootMeansR <- R %*% multiplierR / N
      bootMeansQ <- Q %*% multiplierQ / M

      #----- Estimate the variance from the samples
      if( method == "regular" ){
        data.sigmaR <- matrixStats::rowVars( R )
        data.sigmaQ <- matrixStats::rowVars( Q )
      }else if( method == "t" ){
        bootSecMoments <- R^2 %*% multiplierR^2 / N
        data.varR      <- ( N / ( N - 1 ) ) * abs( bootSecMoments - bootMeansR^2 )

        bootSecMoments <- Q^2 %*% multiplierQ^2 / M
        data.varQ     <- ( M / ( M - 1 ) ) * abs( bootSecMoments - bootMeansQ^2 )
      }
      data.sigma  <- sqrt( ( 1 + 1 / c ) * data.varR + ( 1 + c ) * data.varQ )

      #----- Compute bootstrap distribution of the maximum
      distVec <- sqrt( N + M ) * apply( abs( ( bootMeansR + bootMeansQ ) / data.sigma ), 2, max )

    }else if( D > 1 ){
      #----- Compute bootstrap means
      bootMeansR = array( matrix( R, prod( dimR[ -( D+1 ) ] ), N ) %*%
                            multiplierR, dim = c( dimR[ -( D+1 ) ], Mboots ) ) / N
      bootMeansQ = array( matrix( Q, prod( dimR[ -( D+1 ) ] ), M ) %*%
                            multiplierQ, dim = c( dimQ[ -( D+1 ) ], Mboots ) ) / M

      #----- Estimate the variance from the sample
      if( method == "regular" ){
        data.varR <- apply( R, 1:D, var )
        data.varQ <- apply( Q, 1:D, var )
        data.sigma  <- sqrt( ( 1 + 1 / c ) * data.varR + ( 1 + c ) * data.varQ )

        distVec    <- apply( array( as.vector( abs( bootMeans ) ) /
                                      as.vector( data.sigma ),
                                    dim = c( dimR[ -( D+1 ) ], Mboots ) ), D + 1, max ) * sqrt( N + M - 2 )
      }else if( method == "t" ){
        #----- bootstrapped variance for R
        bootSecMoments <- array( matrix( R^2, prod( dimR[ -( D+1 ) ] ), N ) %*%
                                   multiplierR^2,
                                 dim = c( dimR[ -( D+1 ) ], Mboots ) ) / N
        data.varR <- ( N / ( N - 1 ) ) * abs( bootSecMoments - bootMeansR^2 )

        #----- bootstrapped variance for Q
        bootSecMoments <- array( matrix( Q^2, prod( dimQ[ -( D+1 ) ] ), N ) %*%
                                   multiplierQ^2,
                                 dim = c( dimQ[ -( D+1 ) ], Mboots ) ) / M
        data.varQ <- ( M / ( M - 1 ) ) * abs( bootSecMoments - bootMeansQ^2 )

        #----- bootstrapped variance of the asymptotic process
        data.sigma  <- sqrt( ( 1 + 1 / c ) * data.varR + ( 1 + c ) * data.varQ )

        #----- max statistic
        distVec <- apply( abs( bootMeansR + bootMeansQ )  / data.sigma,
                          D+1,
                          max )  * sqrt(N+M)
      }
    }
  }

  q = quantile( distVec, 1 - alpha, type = 8 )

  #----- Return quantile and bootstrap distribution
  return( list( z       = distVec,
                q       = q,
                samples = bootMeans / data.sigma * sqrt( N ) ) )
}
