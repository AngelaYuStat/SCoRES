data("climate_data")

test_that("Input validation: sp_list", {
  expect_error(SCB_gls_geospatial(), "sp_list")                      # not supplied
  expect_error(SCB_gls_geospatial(data.frame(x = 1)), "list")        # not a list
  expect_error(SCB_gls_geospatial(list(x=1, y=2, z=3)), "x.*y.*obs") # missing names
  expect_error(SCB_gls_geospatial(list(x = "a", y = 1:2, obs = array(0, dim = c(1,2,3)))),
    "numeric")# x is not numeric
  expect_error(SCB_gls_geospatial(list(x = 1:2, y = factor(c(1,2)), obs = array(0, dim = c(2,2,3)))),
    "numeric")# y is not numeric
  expect_error(SCB_gls_geospatial(list(x = matrix(1:2, ncol = 1), y = 1:2, obs = array(0, dim = c(2,2,3)))),
    "numeric vector")# x/y are matrix/array, not vector
  expect_error(SCB_gls_geospatial(list(x = 1:2, y = matrix(1:2, ncol = 1), obs = array(0, dim = c(2,2,3)))),
    "numeric vector")# x/y are matrix/array, not vector
  expect_error(SCB_gls_geospatial(list(x = 1:2, y = 1:2, obs = rnorm(8))),
    "numeric array")# obs not 3D
  expect_error(SCB_gls_geospatial(list(x = 1:3, y = 1:2, obs = array(0, dim = c(4,2,5)))),
    "dims mismatch|first two dims")# obs is 3D but miss match x/y
  expect_error(SCB_gls_geospatial(list(x = 1:2, y = 1:2, obs = array(0, dim = c(2,2,0)))),
    "third dimension.*>= 1")# obs has third dimension n = 0
})

test_that("Input validation: level/alpha/N", {
  expect_error(SCB_gls_geospatial(climate_data$Z, level = "a"), "level")
  expect_error(SCB_gls_geospatial(climate_data$Z, level = c(1,2)), "level")
  expect_error(SCB_gls_geospatial(climate_data$Z, N = -5), "positive")
  expect_error(SCB_gls_geospatial(climate_data$Z, alpha = 1.5), "alpha")
})

test_that("Input validation: data_fit", {
  sp <- list(x = 1:2, y = 1:2, obs = array(rnorm(2*2*5), c(2,2,5))) # 2x2x5
  df <- cbind(1, rnorm(5)) # 5x2
  w <- c(1,0) # 1x2
  groups <- rep(1,5)
  V  <- array(diag(5), dim = c(2,2,5,5)) # 2x2x5x5

  # data_fit not matrix/array
  expect_error(SCB_gls_geospatial(sp, data_fit = data.frame(x = 1)), "matrix or array")
  # data_fit rows mismatch
  expect_error(SCB_gls_geospatial(sp, data_fit = matrix(1, 4, 2)), "rows.*data_fit")
  # w mismatch
  expect_error(SCB_gls_geospatial(sp, data_fit = df, w=c(1,0,0)), "length of `w`")
  expect_error(SCB_gls_geospatial(sp, data_fit = df, w=matrix(c(1, 2, 3), ncol = 1)), "number of rows in `w`")
  # groups length mismatch
  expect_error(SCB_gls_geospatial(sp, data_fit = df, w=w, V=V, groups=rep(1,4)), "groups")

  # V last two dims not square / not n
  V_bad <- array(0, dim = c(2,2,5,4))
  expect_error(SCB_gls_geospatial(sp, data_fit = df, w=w, V=V_bad, groups=groups), "last two dimensions of `V`")
  # neither correlation nor V
  #expect_error(SCB_gls_geospatial(sp, data_fit = df, w=w, groups=groups), "provide one of 'correlation' and 'V'")
})

test_that("Input validation: mask", {
  x <- 1:3; y <- 1:2; n <- 4
  sp <- list(x = x, y = y, obs = array(rnorm(length(x)*length(y)*n),
                                       dim = c(length(x), length(y), n)))
  df <- cbind(1, rnorm(n))
  w <- c(1,0)
  groups <- rep(1, n)

  expect_error(SCB_gls_geospatial(sp, data_fit = df, w = w, groups = groups, mask = 1:6),
    "matrix or array") # mask is not an array/matrix
  mask_3d <- array(1, dim = c(3,2,2))
  expect_error(SCB_gls_geospatial(sp, data_fit = df, w = w, groups = groups, mask = mask_3d),
    "2-dimensional")# mask is an array but not 2D

  mask_wrong <- matrix(1, nrow = 4, ncol = 2)  # should be 3 x 2
  expect_error(SCB_gls_geospatial(sp, data_fit = df, w = w, groups = groups, mask = mask_wrong),
    "dimensions")# mask is 2D but wrong dimensions
})


test_that("Function works well", {
  res <- SCB_gls_geospatial(sp_list = climate_data$Z, level = 2, data_fit = climate_data$X,
                  w = c(1,0,0,0), correlation = climate_data$correlation,
                  mask = climate_data$mask, alpha = 0.1)

  expect_true(is.list(res))
  expect_equal(dim(res$scb_up), c(length(climate_data$Z$x), length(climate_data$Z$y)))
  expect_equal(dim(res$scb_low), c(length(climate_data$Z$x), length(climate_data$Z$y)))
  expect_true(all(res$scb_up[!is.na(res$scb_up)] >= res$scb_low[!is.na(res$scb_low)]))
})
