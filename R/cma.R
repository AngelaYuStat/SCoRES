cma = function(object, model, time, range, groups = NULL, subject = NULL, nboot = NULL){
  # object: a fitted functional outcome regression model object
  # model: character, specify the name of the fitted model associated with the object
  # time: a character specify the name of the time variable
  # range: a numerical vector represents the time region interested
  # group: a vector containing names of the variables corresponding to the groups of interest, null represents the reference group
  # subject: character, specify the name of the random effects for subject
  # nboot: number of samples

  mod_coef <- object$coefficients

  if(model %in% c("gam", "bam")){
    mod_cov <- vcov.gam(object) # containing all variance-covariance info for object
  }

  # get all names for terms
  coef_names <- names(mod_coef)

  # get index of the initial terms
  idx <- grep("^(Intercept)$|^s\\(seconds\\)\\.[0-9]+$", coef_names)

  # initialize dataframe
  predictors <- all.vars(formula(object))[-1]

  df_pred <- as_tibble(
    setNames(
      lapply(predictors, function(x) 0),
      predictors
    )
  )

  if(!is.null(group)){
    # for group interested
    group_idx <- unlist(lapply(group, function(g) grep(g, coef_names)))
    idx <- sort(unique(c(idx, group_idx)))
    df_pred[group] <- 1
  }

  if(!is.null(subject)){
    df_pred[subject] <- model.frame(object)[[subject]][1]
  }

  if(is.null(range)){
    s_pred <- unique(model.frame(object)[time])
  }else{
    s_pred <- sort(unique(range))
  }

  df_pred <- df_pred[rep(1, length(s_pred)), ]
  df_pred[[time]] <- s_pred

  # exract means and variance for group interested
  mod_coef <- mod_coef[idx]
  mod_cov <- mod_cov[idx, idx]

  # prepare design matrix
  lpmat <- predict(object, newdata=df_pred, se.fit=TRUE, type = "lpmatrix")

  # get mean response by groups with standard errors
  pred_df <- df_pred %>% mutate(mean = c(lpmat %*% mod_coef), se = c(sqrt(diag(lpmat %*% mod_cov %*% t(lpmat)))))

  lpmat <- lpmat[, idx]

  # Number of bootstrap samples (B)
  if(is.na(nboot)){
    nboot <- 1e4
  }

  # Set up container for bootstrap
  yhat_boot <- matrix(NA, nboot, length(s_pred))

  # Do the bootstrap

  for (i in 1:nboot) {
    beta_boot_i <- MASS::mvrnorm(n = 1, mu = mod_coef, Sigma = mod_cov)
    yhat_boot[i, ] <- lpmat %*% beta_boot_i
  }

  # Find the max statistic
  dvec_non <- apply(yhat_boot, 1, function(x) max(abs(x - pred_df_non$mean) / pred_df_non$se))

  # Set up container for bootstrap
  yhat_boot <- matrix(NA, nboot, length(s_pred))

  for (i in 1:nboot) {
    beta_boot_i <- MASS::mvrnorm(n = 1, mu = mod_coef_use, Sigma = mod_cov_use)
    yhat_boot[i, ] <- lpmat_use %*% beta_boot_i
  }

  # Find the max statistic
  dvec <- apply(yhat_boot, 1, function(x) max(abs(x - pred_df$mean) / pred_df$se)) # Substract mean estimate and divided by Df (element by element)

  # Get 95% global confidence band
  Z_global <- quantile(dvec, 0.95)
  y_hat_LB_global <- pred_df$mean - Z_global * pred_df$se
  y_hat_UB_global <- pred_df$mean + Z_global * pred_df$se

  global_df <- list(
    yhat = pred_df$mean,
    time = s_pred,
    se_hat = pred_df$se,
    scb_low = y_hat_LB_global,
    scb_uo = y_hat_UB_global,
    type = "Global Confidence Interval (CMA)"
  )

  return(global_df)
}
