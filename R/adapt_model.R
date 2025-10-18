# ============================================================
# R/adapt_model.R
# Lightweight wrapper to choose and apply a model adaptation method.
# ============================================================

# Expects:
# - method$type: "none", "recalibrate", or "revise"
# - Any other fields in `method` are passed through to the specific adapter.
# - data_so_far: data up to (and including) the current block (for fitting)
# - model: the current model object (whatever your code uses)
#
# Returns:
# - An updated model object (or the original if type == "none" or unknown)

adapt_model <- function(df, method, sim_out) {
  
  if (identical(method$type, "none")) {
    return(NULL)
  }
  
  if (identical(method$type, "recalibrate")) {
    return(adapt_model_recalibrate(
      df,
      sim_out,
      coef_original      = method$coef_original,
      shrinkage_factor   = method$shrinkage_factor
    ))
  }
  
  if (identical(method$type, "revise")) {
    return(adapt_model_revise(
      df,
      sim_out,
      coef_original      = method$coef_original,
      shrinkage_factor   = method$shrinkage_factor
    ))
  }
  
  # Unrecognized type: leave coefficients unchanged
  NULL
}

adapt_model_recalibrate <- function(df,
                                    sim_out,
                                    coef_original,
                                    shrinkage_factor = 5000) {
  
  
  # Pull & prep data
  y  <- sim_out$results$outcome
  z  <- as.integer(sim_out$results$treat)             
  id <- as.integer(sim_out$results$person_id)
  r  <- as.numeric(sim_out$risks_mat[id, 1, drop = FALSE]) # baseline risk

  # Shrinkage
  n <- length(y)
  shrinkage <- 1 / (1 + n / shrinkage_factor) 
  
    
  # Get linear predictor
  lp <- log(-log(1-r))
  
  # Group aligned to same ids
  g  <- df$group[id]
  
  # Starting coefficients for GLM
  start <- c(
    "(Intercept)" = 0,
    "lp"          = 1,
    "z"           = 0
     )
  
  # Fit GLM
  fit <- glm(y ~ ., family = binomial(link = "cloglog"), start = start,
             data = data.frame(y = y, lp = lp, z = z))
  
  # Coefficients from GLM fit
  co <- coef(fit)
  
  # Grab intercept from relevant GLM model
  gamma0 <- as.numeric(co["(Intercept)"])
  
  # Grab slope from relevant GLM model (first non-intercept term)
  gamma1 <- as.numeric(co["lp"])
  
  # Compute kappa = s + (1 - s) * slope
  kappa <- shrinkage + (1 - shrinkage) * gamma1
  
  # Initialize new coefficients
  coef_adjusted <- coef_original

  # Update group-specific coefficients
  for (grp in names(coef_original)) {
    
    # Overwrite betas
    coef_adjusted[[grp]]$betas <- coef_original[[grp]]$betas * kappa
    
    # Overwrite mean_x
    s0 <- coef_original[[grp]]$s0
    coef_adjusted[[grp]]$mean_x <- kappa * coef_original[[grp]]$mean_x -
                                (1 - shrinkage) * gamma0 -
                                kappa * log(-log(s0)) + log(-log(s0))
  }

  # Return weighted average of two thresholds
  return( coef_adjusted )
  
}


adapt_model_revise <- function(df,
                               sim_out,
                               coef_original,
                               shrinkage_factor = 5000) {
  
  # Pull and prepare data
  y  <- sim_out$results$outcome
  z  <- as.integer(sim_out$results$treat)
  id <- as.integer(sim_out$results$person_id)
  g    <- df$group[id]
  
  # Shrinkage
  n <- length(y)
  shrinkage <- 1 / (1 + n / shrinkage_factor) 
  
  # Covariate names
  cov_names <- names(coef_original[[1]]$betas)
  
  # Build design matrix
  X_fit <- data.frame(df[id, cov_names, drop = FALSE], z = z)   

  # Fit with cloglog link
  fit <- glm(y ~ ., family = binomial(link = "cloglog"),
             data = X_fit)

  # Recover new coefficients
  bhat <- coef(fit)
  
  # Initialize new coefficients
  coef_adjusted <- coef_original
  
  # Loop over groups
  for (grp in names(coef_adjusted)) {
    

    # Relevant indices of persons in group
    idx   <- which(g == grp)
    
    # Old coefficients
    bold   <- coef_original[[grp]]$betas[cov_names]
    b0_old <- log(-log(coef_original[[grp]]$s0)) - coef_original[[grp]]$mean_x
    
    # Align to available coefs (dropped cols get 0 pull)
    bhat_slopes <- setNames(rep(0, length(cov_names)), cov_names)
    common <- intersect(names(bhat_slopes), names(bhat))
    bhat_slopes[common] <- bhat[common]
    coef_adjusted[[grp]]$betas <-        shrinkage * coef_original[[grp]]$betas[cov_names] +
                                 (1 - shrinkage) * bhat_slopes
    
    # Shrink intercept to old values (and get implied mean_x)
    s0 <- coef_original[[grp]]$s0
    b0_new <- bhat["(Intercept)"]
    coef_adjusted[[grp]]$mean_x <-  shrinkage * coef_original[[grp]]$mean_x -
                                  (1-shrinkage) * b0_new -
                                  shrinkage * log(-log(s0)) + log(-log(s0))

  }
  
  coef_adjusted
}



