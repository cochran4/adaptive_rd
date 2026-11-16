# Parameters
BANDWIDTH = .02

#------------------------------------------------------------------------------#
# estimate_spline: one GLM with separate risk splines and shared PC effects.
# - Assignment groups get separate natural spline bases (knots/boundaries
#   learned from each group's risk_last distribution).
# - Each basis is then evaluated on the full sample and masked by group so the
#   columns only contribute within their group.
# - History is residualized on the complementary log-log scale; its PCs
#   are added linearly as shared adjustments. Splines use probability risk.
#
# Example:
#   fit_out <- estimate_spline(sim_out, var_explained = 0.90, df_ns = 3,
#                              family = binomial("logit"))
#   fit_out <- estimate_spline(sim_out, family = gaussian("identity"))
#------------------------------------------------------------------------------#
estimate_spline <- function(sim_out, family, var_explained = 0.90, df_ns = 2) {
  # Prepare outcomes, assignments, and risk histories for sampled patients.
  y  <- sim_out$results$outcome
  z  <- as.integer(sim_out$results$treat)  # Assigned = 1; unassigned = 0.
  id <- as.integer(sim_out$results$person_id)

  # Rows follow simulated patients; columns follow deployed model versions.
  R  <- sim_out$risks_mat[id, , drop = FALSE]

  # Keep final probability risk for outcome splines and kernel averaging.
  risk_last <- R[, ncol(R)]

  # Use the PCE linear predictor scale for residualization, not the margins
  # R - C. The observed threshold history is common to all patient profiles.
  # Clip probabilities away from 0 and 1 to keep the transformation finite.
  eta <- log(-log1p(-pmin(
    pmax(R, .Machine$double.eps),
    1 - .Machine$double.eps
  )))

  eta_last <- eta[, ncol(eta)]
  eta_hist <- if (ncol(eta) > 1) eta[, -ncol(eta), drop = FALSE] else NULL

  # Summarize history variation remaining after adjustment for final risk.
  pcs <- NULL
  pca <- NULL
  K <- 0
  cum_var_total <- 1.0

  if (!is.null(eta_hist)) {
    # Regress every earlier transformed score on an intercept and final score.
    # Affine recalibrations leave no residual variation beyond numerical noise.
    Xd   <- cbind(1, eta_last)
    B    <- qr.solve(crossprod(Xd), crossprod(Xd, eta_hist))
    Yres <- eta_hist - Xd %*% B

    # Remove columns with negligible residual variance.
    keep <- apply(Yres, 2, var) > 1e-12
    Yres <- Yres[, keep, drop = FALSE]

    if (ncol(Yres) > 0) {
      # Center residuals without standardizing their variances.
      pca      <- prcomp(Yres, center = TRUE, scale. = FALSE)
      eig_res  <- pca$sdev^2
      var_r    <- var(eta_last)
      total    <- var_r + sum(eig_res)

      # Count final transformed risk first, then retain the fewest PCs needed
      # to reach var_explained of the combined variance.
      cum_prop <- cumsum(c(var_r, eig_res) / total)
      need     <- which(cum_prop >= var_explained)[1]
      K        <- max(0, need - 1L)

      if (K > 0)
        pcs <- setNames(
          as.data.frame(pca$x[, seq_len(K), drop = FALSE]),
          paste0("PC", seq_len(K))
        )

      cum_var_total <- cum_prop[1 + K]
    }
  }

  # Construct separate natural spline bases for the assignment groups.
  x0 <- risk_last[z == 0]
  x1 <- risk_last[z == 1]

  # Learn knots and boundaries from observed final risks within each group.
  b0_tmp <- splines::ns(x0, df = df_ns)
  k0     <- attr(b0_tmp, "knots")
  bk0    <- attr(b0_tmp, "Boundary.knots")

  b1_tmp <- splines::ns(x1, df = df_ns)
  k1     <- attr(b1_tmp, "knots")
  bk1    <- attr(b1_tmp, "Boundary.knots")

  # Evaluate both bases for all patients, using their respective group knots.
  B0_full <- splines::ns(risk_last, knots = k0, Boundary.knots = bk0)
  B1_full <- splines::ns(risk_last, knots = k1, Boundary.knots = bk1)
  colnames(B0_full) <- paste0("B0_", seq_len(ncol(B0_full)))
  colnames(B1_full) <- paste0("B1_", seq_len(ncol(B1_full)))

  # Use numeric assignment so the formula can mask the two spline bases.
  t <- as.integer(z)

  # Assemble predictors for the joint outcome model.
  dat <- data.frame(
    y = y,
    t = t,
    risk_last = risk_last,  # Retained for prediction and kernel averaging.
    B0_full,
    B1_full
  )

  if (!is.null(pcs)) dat <- cbind(dat, pcs)

  # Each spline contributes only in its assignment group; PCs are shared.
  #   y ~ t + I((1 - t) * B0_1) + ... + I((1 - t) * B0_K)
  #         + I(t * B1_1) + ... + I(t * B1_K)
  #         + PCs
  b0_terms <- paste0("I((1 - t) * ", colnames(B0_full), ")")
  b1_terms <- paste0("I(t * ",         colnames(B1_full), ")")
  pc_terms <- if (!is.null(pcs)) names(pcs) else NULL

  rhs <- c("t", b0_terms, b1_terms, pc_terms)
  fml <- stats::as.formula(paste("y ~", paste(rhs, collapse = " + ")))

  # Fit all observed outcomes without kernel weights.
  fit_joint <- glm(fml, family = family, data = dat)

  # Save the fit and spline definitions for prediction under either assignment.
  list(
    fit       = fit_joint,
    basis_df  = df_ns,
    knots0    = k0,
    bknots0   = bk0,
    knots1    = k1,
    bknots1   = bk1,
    pca       = pca,
    K         = K,
    cum_var_total = cum_var_total,
    data      = dat
  )
}

#------------------------------------------------------------------------------#
# conditional_curves: kernel-averaged mean predictions under each assignment.
# Predict both arms at observed risks and history PCs; each target x
# determines only the kernel weights used to average those predictions.
#
# Args:
#   fit_out  : fitted model and prediction data returned by estimate_spline()
#   x_seq    : numeric grid of x values
#   level    : CI level (default 0.95)
#   Uses the shared BANDWIDTH for local averaging over final probability risk.
#
# Returns:
#   data.frame(risk_last, mean, l, u, grp) with grp in {"Control","Treated"}
#------------------------------------------------------------------------------#
conditional_curves <- function(fit_out, x_seq, level = 0.95) {

  # Stabilized kernel function
  Kg <- function(u) {
    ex <- -0.5 * u^2
    ex <- ex - max(ex)
    exp(ex)
  }

  # Retrieve the joint model and observed patient predictors.
  fit   <- fit_out$fit
  dat   <- fit_out$data

  fam   <- fit$family
  beta  <- stats::coef(fit)
  Vb    <- 0.5 * (stats::vcov(fit) + t(stats::vcov(fit)))     # symmetrize
  trm   <- stats::delete.response(stats::terms(fit))
  zcrit <- stats::qnorm((1 + level)/2)

  # Use the same final-risk kernel weights under both assignments.
  Xobs <- dat$risk_last
  h    <- BANDWIDTH

  # Set assignment only. Keep each patient's observed risk spline bases and PCs.
  mu_provider <- function(tval) {
    nd <- dat
    nd$t <- as.numeric(tval)

    Xmat  <- stats::model.matrix(trm, nd)
    eta_i <- drop(Xmat %*% beta)
    mu_i  <- fam$linkinv(eta_i)
    list(mu_i = mu_i, eta_i = eta_i, Xmat = Xmat)
  }

  # Delta-method interval for the weighted mean, constructed on the link scale.
  ci_g_link <- function(w, prov) {
    mu_i  <- prov$mu_i
    eta_i <- prov$eta_i
    Xmat  <- prov$Xmat

    # Average outcomes, then calculate the link-scale gradient.
    mhat   <- sum(w * mu_i)
    m_c    <- mhat
    eta_m  <- fam$linkfun(m_c)              # g(mhat)
    gprime <- 1 / fam$mu.eta(eta_m)         # inverse-function theorem

    dmu_i  <- fam$mu.eta(eta_i)
    grad_m <- drop(t(Xmat) %*% (w * dmu_i))
    grad_g <- gprime * grad_m

    se_g <- sqrt(pmax(drop(t(grad_g) %*% Vb %*% grad_g), 0))
    c(lo = fam$linkinv(eta_m - zcrit * se_g),
      hi = fam$linkinv(eta_m + zcrit * se_g),
      mhat = mhat)
  }

  # Predict once per assignment; vary only kernel weights over target risks.
  eval_group <- function(tval, label) {
    prov <- mu_provider(tval)
    one_x <- function(x) {
      w <- Kg((Xobs - x) / h)
      w <- w / sum(w)
      ci   <- ci_g_link(w, prov)
      c(mean = ci[["mhat"]], l = ci[["lo"]], u = ci[["hi"]])
    }
    out <- t(vapply(x_seq, one_x, numeric(3)))
    data.frame(risk_last = x_seq,
               mean = out[, "mean"], l = out[, "l"], u = out[, "u"],
               grp = label, row.names = NULL)
  }

  rbind(
    eval_group(0, "Control"),
    eval_group(1, "Treated")
  )
}

# -- Gradient of the kernel-averaged arm mean, with observed risks and PCs.
.grad_at <- function(fit_out, x, tval) {
  fit <- fit_out$fit
  dat <- fit_out$data
  fam <- fit$family
  beta <- stats::coef(fit)
  trm <- stats::delete.response(stats::terms(fit))
  # kernel weights on full X
  Kg <- function(u) {
    ex <- -0.5 * u^2
    ex <- ex - max(ex)
    exp(ex)
  }
  Xobs <- dat$risk_last
  h    <- BANDWIDTH
  w <- Kg((Xobs - x) / h)
  w <- w / sum(w)

  # Use the same observed-risk predictions as conditional_curves().
  nd <- dat
  nd$t <- as.numeric(tval)

  Xmat <- stats::model.matrix(trm, nd)
  eta  <- as.numeric(Xmat %*% beta)
  dmu  <- fit$family$mu.eta(eta)
  drop(t(Xmat) %*% (w * dmu))   # ∑ w_i μ'(η_i) X_i
}

# -- ATE at cutoff using conditional_curves() for means; joint delta for CI
ate_at_threshold <- function(fit_out, sim_out, ci = 0.95) {
  # Evaluate at the final threshold using the joint coefficient covariance.
  c_star <- tail(sim_out$results$threshold_used, 1)
  Vb <- 0.5 * (stats::vcov(fit_out$fit) + t(stats::vcov(fit_out$fit)))
  z  <- stats::qnorm((1 + ci)/2)

  # Obtain the two kernel-averaged means.
  cc <- conditional_curves(fit_out, x_seq = c_star, level = ci)
  m0 <- subset(cc, grp == "Control")$mean
  m1 <- subset(cc, grp == "Treated")$mean

  # Differentiate the contrast jointly to include covariance between arm means.
  g0 <- .grad_at(fit_out, c_star, tval = 0)
  g1 <- .grad_at(fit_out, c_star, tval = 1)
  g  <- g1 - g0
  se <- sqrt(pmax(drop(t(g) %*% Vb %*% g), 0))

  ate <- m1 - m0
  ci_lo <- ate - z * se
  ci_hi <- ate + z * se

  list(cutoff = c_star, ate = ate, se = se, ci = c(lo = ci_lo, hi = ci_hi),
       treated = m1, control = m0)
}

#------------------------------------------------------------------------------#
# conditional_curves_truth
#
# Purpose:
#   Estimates the *true* conditional mean curves
#   E[Y(a) | R = x] for a ∈ {0,1} over a grid x_seq, by plugging in the
#   **actual (truth) outcome mean model** supplied via `mean_fn`.
#   No model fitting is done here; we purely evaluate the DGP mean.
#
# Arguments:
#   cohort_df : data.frame of individuals (rows)
#   risk_fn   : function(df) -> numeric vector R in [0,1]; risk for each row
#   mean_fn   : function(cohort_df, idx, risk0, treat) -> numeric mean vector
#               Must follow your required outcome_fn signature and implement
#               the TRUE mean model of Y given (R, treat) for these data.
#   x_seq     : grid of target final risks for local averaging
#   h         : positive bandwidth for kernel smoothing over R
#
# Returns:
#   data.frame with columns:
#     risk  : grid value x
#     mean  : estimated E[Y(a) | R = x] using truth mean model
#     group : "Control" (a=0) or "Treated" (a=1)
#------------------------------------------------------------------------------#
conditional_curves_truth <- function(
    cohort_df,
    risk_mat,
    mean_fn,
    x_seq
) {

  # Compute first and last risks
  R_first <- as.numeric(risk_mat[,1])
  R_last  <- as.numeric(risk_mat[,ncol(risk_mat)])
  n       <- NROW(cohort_df)
  idx     <- seq_len(n)

  # Set bandwidth
  h    <- BANDWIDTH

  # Kernel helper function (kernel around last risk)
  Kg <- function(u) {
    ex <- -0.5 * u^2
    ex <- ex - max(ex)
    exp(ex)
  }
  kweights <- function(x) {
    w <- Kg((R_last - x)/h)
    sw <- sum(w)
    w / sum(w)
  }
  #    Evaluate the TRUE mean outcome for all rows under treatment a = tval.
  #    Note: we pass `risk0 = R_obs` to match the required signature; the mean
  #    function is expected to implement the truth (no fitting).
  eval_group <- function(tval, label) {
    mu_true_i <- mean_fn(
      cohort_df,
      idx,                # use all rows
      R_first,            # baseline risk per row (same length as idx)
      rep_len(tval, n)    # constant treatment assignment
    )

    # For each target x, take the kernel-weighted average of mu_true_i
    m_x <- vapply(
      x_seq,
      function(x) sum(kweights(x) * mu_true_i),
      numeric(1)
    )

    data.frame(
      risk_last  = x_seq,
      mean  = m_x,
      group = label,
      row.names = NULL
    )
  }
  # Return both curves: control (a=0) and treated (a=1)
  rbind(
    eval_group(0, "Control"),
    eval_group(1, "Treated")
  )
}

# -- ATE at cutoff using conditional_curves_truth() for means
ate_at_threshold_truth <- function(sim_out, cohort_df, mean_fn) {

  # Last threshold used
  c_star <- tail(sim_out$results$threshold_used, 1)

  # Obtain the two kernel-averaged means.
  cc <- conditional_curves_truth(cohort_df,
                                 sim_out$risks_mat,
                                 mean_fn,
                                 c_star)

  # Grab the relevant values
  m0 <- subset(cc, group == "Control")$mean
  m1 <- subset(cc, group == "Treated")$mean

  # Local ATE (ground truth)
  ate <- m1 - m0

  list(cutoff = c_star, ate = ate,
       treated = m1, control = m0)
}

ate_naive <- function(sim_out, family = gaussian()) {
  y <- sim_out$results$outcome
  z <- sim_out$results$treat
  z <- as.integer(z)  # in case it's logical

  fam_name <- family$family

  if (fam_name == "binomial") {
    p1 <- mean(y[z == 1])
    p0 <- mean(y[z == 0])
    n1 <- sum(z == 1)
    n0 <- sum(z == 0)

    est <- p1 - p0
    se  <- sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
    ci  <- est + c(-1, 1) * qnorm(0.975) * se

  } else if (fam_name == "gaussian") {
    test <- t.test(y ~ z, var.equal = FALSE)
    est  <- diff(test$estimate)
    ci   <- test$conf.int

  } else {
    stop("Unsupported family: ", fam_name)
  }

  list(ate = est, ci = c(lo = ci[1], hi = ci[2]))
}

ate_outreg <- function(sim_out, cohort_df, family, separate_arms = FALSE) {

  # Recover raw covariates from the PCE predictor columns.
  raw_df <- cohort_df %>%
    mutate(
      # Undo the logs
      age         = exp(ln_age),
      total_chol  = exp(ln_total_cholest),
      hdlC        = exp(ln_hdlC),

      # Reconstruct systolic_bp using whichever BP variable is non-zero
      systolic_bp = exp(ln_treated_BP + ln_untreated_BP ),

      # Reconstruct treated_bp: TRUE if treated version was used
      treated_bp = ln_treated_BP > 0,

      # Keep smoker and diabetes as-is (already binary)
      smoker   = smoker,
      diabetes = diabetes
    ) %>%
    select(age, total_chol, hdlC, systolic_bp, treated_bp, smoker, diabetes)

  # Extract relevant data
  y  <- sim_out$results$outcome
  z  <- as.integer(sim_out$results$treat)
  id <- as.integer(sim_out$results$person_id)
  X  <- raw_df[id, , drop = FALSE]
  risk_last <- sim_out$risks_mat[id, ncol(sim_out$risks_mat), drop = TRUE]
  c_star    <- tail(sim_out$results$threshold_used, 1)

  # Gaussian kernel weights centered at c_star
  kernel <- function(x) exp(-0.5 * (x / BANDWIDTH)^2)
  wts <- kernel(risk_last - c_star)

  df_fit <- data.frame(y = y, z = z, X)
  if (separate_arms) {
    # Each arm has its own intercept and covariate effects.
    fit1 <- glm(y ~ ., data = data.frame(y = y[z == 1], X[z == 1, ]),
                family = family)
    fit0 <- glm(y ~ ., data = data.frame(y = y[z == 0], X[z == 0, ]),
                family = family)

    # Predict both outcomes for the same patients before local averaging.
    mu1 <- predict(fit1, newdata = X, type = "response")
    mu0 <- predict(fit0, newdata = X, type = "response")
  } else {
    # Original pooled model with shared covariate effects.
    fit <- glm(y ~ z + ., data = df_fit, family = family)

    df_cf1 <- df_fit
    df_cf1$z <- 1
    mu1 <- predict(fit, newdata = df_cf1, type = "response")
    df_cf0 <- df_fit
    df_cf0$z <- 0
    mu0 <- predict(fit, newdata = df_cf0, type = "response")
  }

  # Weighted average difference at c_star
  est <- weighted.mean(mu1 - mu0, wts)

  # Comparator confidence intervals were not calculated.
  ci <- c(lo = NaN, hi = NaN)

  list(ate = est, ci = ci)
}

# Separate-arm outcome regression, using the same local ATE weights.
ate_outreg_separate <- function(sim_out, cohort_df, family) {
  ate_outreg(sim_out, cohort_df, family, separate_arms = TRUE)
}

ate_ipw <- function(sim_out, cohort_df, family) {

  # Back-transform raw variables
  raw_df <- cohort_df %>%
    mutate(
      age         = exp(ln_age),
      total_chol  = exp(ln_total_cholest),
      hdlC        = exp(ln_hdlC),
      systolic_bp = exp(ln_treated_BP + ln_untreated_BP),
      treated_bp  = ln_treated_BP > 0,
      smoker      = smoker,
      diabetes    = diabetes
    ) %>%
    select(age, total_chol, hdlC, systolic_bp, treated_bp, smoker, diabetes)

  # Extract sim results
  y  <- sim_out$results$outcome
  z  <- as.integer(sim_out$results$treat)
  id <- as.integer(sim_out$results$person_id)
  X  <- raw_df[id, , drop = FALSE]
  risk_last <- sim_out$risks_mat[id, ncol(sim_out$risks_mat), drop = TRUE]
  c_star    <- tail(sim_out$results$threshold_used, 1)

  # Compute kernel weights centered at c_star
  kernel <- function(x) exp(-0.5 * (x / BANDWIDTH)^2)
  wts <- kernel(risk_last - c_star)

  # Estimate propensity score (Pr[z = 1 | X]) using logistic regression
  df_ps <- data.frame(z = z, X)
  ps_model <- glm(z ~ ., data = df_ps, family = binomial())
  ps <- predict(ps_model, type = "response")

  # Compute the difference between separately normalized weighted means.
  ipw_ate <- function(y, z, ps, wts) {
    treated   <- (z == 1)
    control   <- (z == 0)

    # Normalize inverse-probability and kernel weights within each group.
    mu1 <- sum(wts * y * treated / ps) / sum(wts * treated / ps)
    mu0 <- sum(wts * y * control / (1 - ps)) /
      sum(wts * control / (1 - ps))

    mu1 - mu0
  }

  est <- ipw_ate(y, z, ps, wts)

  # Comparator confidence intervals were not calculated.
  list(ate = est, ci = c(lo = NaN, hi = NaN))
}

ate_aipw <- function(sim_out, cohort_df, family) {

  # Transform covariates back to raw scale
  raw_df <- cohort_df %>%
    mutate(
      age         = exp(ln_age),
      total_chol  = exp(ln_total_cholest),
      hdlC        = exp(ln_hdlC),
      systolic_bp = exp(ln_treated_BP + ln_untreated_BP),
      treated_bp  = ln_treated_BP > 0,
      smoker      = smoker,
      diabetes    = diabetes
    ) %>%
    select(age, total_chol, hdlC, systolic_bp, treated_bp, smoker, diabetes)

  # Extract relevant variables
  y  <- sim_out$results$outcome
  z  <- as.integer(sim_out$results$treat)
  id <- as.integer(sim_out$results$person_id)
  X  <- raw_df[id, , drop = FALSE]
  risk_last <- sim_out$risks_mat[id, ncol(sim_out$risks_mat), drop = TRUE]
  c_star    <- tail(sim_out$results$threshold_used, 1)

  # Compute weights around c_star
  kernel <- function(x) exp(-0.5 * (x / BANDWIDTH)^2)
  wts <- kernel(risk_last - c_star)

  # Fit the assignment model without kernel weights.
  df_ps <- data.frame(z = z, X)
  ps_model <- glm(z ~ ., data = df_ps, family = binomial())
  ps <- predict(ps_model, type = "response")

  # Fit the pooled outcome model without kernel weights.
  df_fit <- data.frame(y = y, z = z, X)
  fit <- glm(y ~ z + ., data = df_fit, family = family)

  df_cf1 <- df_fit
  df_cf1$z <- 1
  mu1 <- predict(fit, newdata = df_cf1, type = "response")

  df_cf0 <- df_fit
  df_cf0$z <- 0
  mu0 <- predict(fit, newdata = df_cf0, type = "response")

  # AIPW components
  treated   <- (z == 1)
  control   <- (z == 0)

  # Predicted contrasts plus inverse-probability-weighted outcome residuals.
  eif <- ((treated * (y - mu1)) / ps + mu1) -
    ((control * (y - mu0)) / (1 - ps) + mu0)

  # AIPW estimate
  est <- weighted.mean(eif, wts)

  # Comparator confidence intervals were not calculated.
  list(ate = est, ci = c(lo = NaN, hi = NaN))
}
