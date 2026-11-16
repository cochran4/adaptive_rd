# ============================================================
# R/adapt_model.R
# Lightweight wrapper to choose and apply a model adaptation method.
# ============================================================

# Select the requested model update using data observed so far.
adapt_model <- function(df, method, sim_out) {
  # NULL tells pce_predict() to use the original coefficients.
  if (identical(method$type, "none")) {
    return(NULL)
  }

  if (identical(method$type, "recalibrate")) {
    return(adapt_model_recalibrate(
      df, sim_out,
      coef_original = method$coef_original,
      recalibration_penalty = method$recalibration_penalty
    ))
  }

  if (identical(method$type, "revise")) {
    return(adapt_model_revise(
      df, sim_out,
      coef_original = method$coef_original,
      revision_penalty = method$revision_penalty
    ))
  }

  NULL
}

# This recalibration method assumes the PCE complementary log-log scale.
# Other model scales require a matching risk transformation, inverse link,
# and conversion of fitted coefficients to the prediction model format.
adapt_model_recalibrate <- function(
    df,
    sim_out,
    coef_original,
    recalibration_penalty = 4
) {
  # Outcomes, assignments, and original risks for all patients observed so far.
  y  <- sim_out$results$outcome
  z  <- as.integer(sim_out$results$treat)
  id <- as.integer(sim_out$results$person_id)
  r  <- as.numeric(sim_out$risks_mat[id, 1, drop = FALSE]) # baseline risk

  n <- length(y)

  # Avoid infinite values when transforming risk to the complementary log-log
  # scale.
  r <- pmin(pmax(r, .Machine$double.eps), 1 - .Machine$double.eps)
  lp <- log(-log1p(-r))

  # Minimize the summed negative log-likelihood plus a fixed penalty.
  objective <- function(gamma) {
    # gamma: calibration intercept, calibration slope, assignment coefficient.
    eta <- gamma[1] + gamma[2] * lp + gamma[3] * z

    # Convert the linear predictor to event probabilities.
    p <- 1 - exp(-exp(eta))

    # Keep log-likelihood calculations finite at numerical endpoints.
    p <- pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
    negative_log_likelihood <- -sum(stats::dbinom(y, 1, p, log = TRUE))

    # Shrink the intercept toward 0 and slope toward 1; do not penalize
    # assignment.
    penalty <- recalibration_penalty / 2 *
      (gamma[1]^2 + (gamma[2] - 1)^2)
    negative_log_likelihood + penalty
  }

  # Start at the original model with no assignment effect.
  initial <- c(gamma0 = 0, gamma1 = 1, gammaA = 0)
  fit <- stats::optim(
    par = initial,
    fn = objective,
    method = "BFGS",
    hessian = TRUE,
    control = list(maxit = 1000)
  )

  gamma0 <- unname(fit$par["gamma0"])
  gamma1 <- unname(fit$par["gamma1"])

  # Record the fitted assignment coefficient; exclude it from risk predictions.
  gammaA <- unname(fit$par["gammaA"])

  # Express the recalibrated predictions in the format used by pce_predict().
  coef_adjusted <- coef_original

  # Apply the same calibration intercept and slope to every PCE group.
  for (grp in names(coef_original)) {
    # Scale the original predictor coefficients by the calibration slope.
    coef_adjusted[[grp]]$betas <- coef_original[[grp]]$betas * gamma1

    # Adjust the centering constant to incorporate the calibration intercept.
    s0 <- coef_original[[grp]]$s0
    coef_adjusted[[grp]]$mean_x <- gamma1 * coef_original[[grp]]$mean_x -
      gamma0 -
      gamma1 * log(-log(s0)) + log(-log(s0))
  }

  # Save fitted parameters and convergence information for reviewing the
  # updates.
  attr(coef_adjusted, "recalibration") <- list(
    n = n,
    gamma0 = gamma0,
    gamma1 = gamma1,
    gammaA = gammaA,
    recalibration_penalty = recalibration_penalty,
    objective = fit$value,
    convergence = fit$convergence
  )

  coef_adjusted
}

# Start every racial group with the White PCE equation for the corresponding
# sex.
# Keep the group names expected by pce_predict().
white_pce_for_all <- function(coef_original) {
  coef_initial <- coef_original
  coef_initial$Black_F <- coef_original$White_F
  coef_initial$Other_F <- coef_original$White_F
  coef_initial$Black_M <- coef_original$White_M
  coef_initial$Other_M <- coef_original$White_M
  coef_initial
}

# This revision method assumes the PCE complementary log-log scale.
# Other model scales require different starting coefficients, an inverse link,
# and a conversion of fitted coefficients to the prediction model format.
adapt_model_revise <- function(
    df,
    sim_out,
    coef_original,
    revision_penalty = 3000
) {
  # Collect outcomes, assignments, and patient characteristics observed so far.
  y  <- sim_out$results$outcome
  z  <- as.integer(sim_out$results$treat)
  id <- as.integer(sim_out$results$person_id)

  # Extract sex from the group label; race is not used in the revised model.
  sex <- ifelse(grepl("_F$", df$group[id]), "Female", "Male")
  n <- length(y)

  # Build the predictor matrix, including an intercept.
  # PCE groups share predictor names; each sex receives its own fitted
  # coefficients.
  cov_names <- names(coef_original$White_F$betas)
  x <- as.matrix(df[id, cov_names, drop = FALSE])
  storage.mode(x) <- "double"
  x <- cbind(`(Intercept)` = 1, x)

  # Convert the starting White PCE equations to the complementary log-log scale.
  # The intercept incorporates the original survival and centering constants.
  white_theta <- function(entry) {
    c(
      `(Intercept)` = log(-log(entry$s0)) - entry$mean_x,
      entry$betas[cov_names]
    )
  }

  # Keep these starting coefficients fixed as the penalty targets at every
  # update.
  theta0_f <- white_theta(coef_original$White_F)
  theta0_m <- white_theta(coef_original$White_M)
  # Number of coefficients per sex, including the intercept.
  p <- length(theta0_f)

  # Minimize the summed negative log-likelihood plus a Gaussian penalty.
  objective <- function(par) {
    # Parameter order: female coefficients, male coefficients, assignment
    # effect.
    theta_f <- par[seq_len(p)]
    theta_m <- par[p + seq_len(p)]
    gamma_a <- par[2 * p + 1L]

    # Apply the coefficients for each patient's sex, then add assignment.
    eta <- numeric(n)
    female <- sex == "Female"
    male <- sex == "Male"
    eta[female] <- as.numeric(x[female, , drop = FALSE] %*% theta_f)
    eta[male] <- as.numeric(x[male, , drop = FALSE] %*% theta_m)
    eta <- eta + gamma_a * z

    # Convert the complementary log-log predictor to an event probability.
    event_probability <- 1 - exp(-exp(eta))

    # Keep probabilities away from numerical endpoints during optimization.
    event_probability <- pmin(
      pmax(event_probability, .Machine$double.eps),
      1 - .Machine$double.eps
    )
    # Sum the Bernoulli log-likelihood contributions from all observed patients.
    negative_log_likelihood <- -sum(
      stats::dbinom(y, 1, event_probability, log = TRUE)
    )

    # Penalize departures from starting coefficients, including the intercepts.
    # Leave the shared assignment coefficient unpenalized.
    penalty <- revision_penalty / 2 * (
      sum((theta_f - theta0_f)^2) +
        sum((theta_m - theta0_m)^2)
    )
    negative_log_likelihood + penalty
  }

  # Start optimization at the White equations with assignment effect zero.
  initial <- c(theta0_f, theta0_m, gammaA = 0)
  fit <- stats::optim(
    par = initial,
    fn = objective,
    method = "BFGS"
  )

  # Extract the fitted coefficients in the same order used by the objective.
  theta_f <- fit$par[seq_len(p)]
  theta_m <- fit$par[p + seq_len(p)]
  gamma_a <- unname(fit$par[2 * p + 1L])

  # Convert fitted intercepts and slopes back to the format used by
  # pce_predict().
  revised_entry <- function(theta, reference_entry) {
    entry <- reference_entry
    # Separate the predictor coefficients from the fitted intercept.
    entry$betas <- stats::setNames(theta[-1L], cov_names)
    # Hold s0 fixed and encode the fitted intercept through mean_x.
    entry$mean_x <- log(-log(entry$s0)) - unname(theta[1L])
    entry
  }

  # Use the same revised equation for all racial groups within each sex.
  female_entry <- revised_entry(theta_f, coef_original$White_F)
  male_entry <- revised_entry(theta_m, coef_original$White_M)
  coef_adjusted <- white_pce_for_all(coef_original)
  coef_adjusted[c("White_F", "Black_F", "Other_F")] <- list(female_entry)
  coef_adjusted[c("White_M", "Black_M", "Other_M")] <- list(male_entry)

  # Save update diagnostics. The returned risk equations exclude gamma_a.
  attr(coef_adjusted, "revision") <- list(
    n = n,
    gammaA = gamma_a,
    revision_penalty = revision_penalty,
    # Euclidean distances from the starting female and male coefficients.
    female_distance = sqrt(sum((theta_f - theta0_f)^2)),
    male_distance = sqrt(sum((theta_m - theta0_m)^2)),
    objective = fit$value,
    convergence = fit$convergence
  )

  coef_adjusted
}
