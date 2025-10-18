# Parameters
BANDWIDTH = .02


#------------------------------------------------------------------------------#
# estimate_spline (single GLM, group-specific ns() bases; no interactions)
# - Control and treated get their OWN natural spline bases (knots/boundaries
#   learned from each group's risk_last distribution).
# - Each basis is then evaluated on the full sample and masked by group so the
#   columns only contribute within their group.
# - PCs from residualized history are added linearly (shared adjustment).
#
# Example:
#   fit_out <- estimate_spline(sim_out, var_explained = 0.90, df_ns = 3,
#                              family = binomial("logit"))
#   fit_out <- estimate_spline(sim_out, family = gaussian("identity"))
#------------------------------------------------------------------------------#
estimate_spline <- function(sim_out, family, var_explained = 0.90, df_ns = 2) {
  #--------------------------- Pull & prep data ------------------------------#
  y  <- sim_out$results$outcome
  z  <- as.integer(sim_out$results$treat)             # 0/1
  id <- as.integer(sim_out$results$person_id)
  R  <- sim_out$risks_mat[id, , drop = FALSE]
  
  risk_last  <- R[, ncol(R)]                          # running variable
  risks_hist <- if (ncol(R) > 1) R[, -ncol(R), drop = FALSE] else NULL
  
  #------------------- PCA on residualized history ---------------------------#
  pcs <- NULL; pca <- NULL; K <- 0; cum_var_total <- 1.0
  if (!is.null(risks_hist)) {
    # residualize history on intercept + risk_last
    Xd   <- cbind(1, risk_last)
    B    <- qr.solve(crossprod(Xd), crossprod(Xd, risks_hist))
    Yres <- risks_hist - Xd %*% B
    
    # drop ~constant columns
    keep <- apply(Yres, 2, var) > 1e-12
    Yres <- Yres[, keep, drop = FALSE]
    
    if (ncol(Yres) > 0) {
      pca      <- prcomp(Yres, center = TRUE, scale. = TRUE)
      eig_res  <- pca$sdev^2
      var_r    <- var(risk_last)
      total    <- var_r + sum(eig_res)
      
      # choose K: risk_last counts first, then PCs until >= var_explained
      cum_prop <- cumsum(c(var_r, eig_res) / total)
      need     <- which(cum_prop >= var_explained)[1]
      K        <- max(0, need - 1L)
      if (K > 0)
        pcs <- setNames(as.data.frame(pca$x[, seq_len(K), drop = FALSE]),
                        paste0("PC", seq_len(K)))
      cum_var_total <- cum_prop[1 + K]
    }
  }
  
  #---------------- Group-specific ns() bases (own knots & boundaries) -------#
  x0 <- risk_last[z == 0]
  x1 <- risk_last[z == 1]
  
  # learn knots/boundaries from each group
  b0_tmp <- splines::ns(x0, df = df_ns)
  k0     <- attr(b0_tmp, "knots")
  bk0    <- attr(b0_tmp, "Boundary.knots")
  
  b1_tmp <- splines::ns(x1, df = df_ns)
  k1     <- attr(b1_tmp, "knots")
  bk1    <- attr(b1_tmp, "Boundary.knots")
  
  # evaluate each basis on the FULL sample using its group's parameters
  B0_full <- splines::ns(risk_last, knots = k0, Boundary.knots = bk0)
  B1_full <- splines::ns(risk_last, knots = k1, Boundary.knots = bk1)
  colnames(B0_full) <- paste0("B0_", seq_len(ncol(B0_full)))
  colnames(B1_full) <- paste0("B1_", seq_len(ncol(B1_full)))
  
  # numeric treatment (0/1) so arithmetic in the formula works
  t <- as.integer(z)   
  
  #-------------------- Assemble design & build formula ----------------------#
  dat <- data.frame(
    y = y,
    t = t,                                # numeric 0/1 used inside I(...)
    risk_last = risk_last,                 # for reference/plots if needed
    B0_full,
    B1_full
  )
  if (!is.null(pcs)) dat <- cbind(dat, pcs)
  
  # Formula:
  #   y ~ t + I((1 - t) * B0_1) + ... + I((1 - t) * B0_K)
  #         + I(t * B1_1) + ... + I(t * B1_K)
  #         + PCs
  b0_terms <- paste0("I((1 - t) * ", colnames(B0_full), ")")
  b1_terms <- paste0("I(t * ",         colnames(B1_full), ")")
  pc_terms <- if (!is.null(pcs)) names(pcs) else NULL
  
  rhs <- c("t", b0_terms, b1_terms, pc_terms)
  fml <- stats::as.formula(paste("y ~", paste(rhs, collapse = " + ")))
  
  fit_joint <- glm(fml, family = family, data = dat)
  
  #----------------------------- Return ---------------------------------------#
  list(
    fit       = fit_joint,               # one GLM w/ group-specific spline columns
    basis_df  = df_ns,
    knots0    = k0,  bknots0 = bk0,      # save parameters for prediction/reuse
    knots1    = k1,  bknots1 = bk1,
    pca       = pca, K = K, cum_var_total = cum_var_total,
    data      = dat
  )
}

#------------------------------------------------------------------------------#
# conditional_curves (joint GLM; both groups in one call, with I((1-t)*B0_*) & I(t*B1_*))
# Uses knots saved in fit_out to recompute bases at each x, then sets t := 0/1.
#
# Args:
#   fit_out  : list returned by estimate_spline() containing fit, data, knots0/1, bknots0/1
#   x_seq    : numeric grid of x values
#   level    : CI level (default 0.95)
#   bandwidth: optional kernel bandwidth on X; default = bw.nrd0 on FULL X
#
# Returns:
#   data.frame(risk_last, mean, l, u, grp) with grp in {"Control","Treated"}
#------------------------------------------------------------------------------#
conditional_curves <- function(fit_out, x_seq, level = 0.95) {

  # Stabilized kernel function
  Kg <- function(u) { ex <- -0.5 * u^2; ex <- ex - max(ex); exp(ex) }

  # pieces from fit_out
  fit   <- fit_out$fit
  dat   <- fit_out$data
  k0    <- fit_out$knots0;   bk0 <- fit_out$bknots0
  k1    <- fit_out$knots1;   bk1 <- fit_out$bknots1
  
  fam   <- fit$family
  beta  <- stats::coef(fit)
  Vb    <- 0.5 * (stats::vcov(fit) + t(stats::vcov(fit)))     # symmetrize
  trm   <- stats::delete.response(stats::terms(fit))
  zcrit <- stats::qnorm((1 + level)/2)
  
  # names of basis columns in the training data
  b0_names <- grep("^B0_", names(dat), value = TRUE)
  b1_names <- grep("^B1_", names(dat), value = TRUE)
  
  # kernel weights on FULL X (same for both groups)
  Xobs <- dat$risk_last
  h    <- BANDWIDTH
 
  # provider: at x and group tval, recompute bases with saved knots & set t := tval
  mu_provider <- function(x, tval) {
    n  <- nrow(dat)
    nd <- dat
    nd$risk_last <- x
    nd$t         <- as.numeric(tval)  # numeric 0/1 used by I((1 - t)*...), I(t*...)
    
    # recompute bases at 'x' using saved knots/boundaries; assign to matching cols
    B0x <- splines::ns(rep(x, n), knots = k0, Boundary.knots = bk0)
    B1x <- splines::ns(rep(x, n), knots = k1, Boundary.knots = bk1)
    colnames(B0x) <- b0_names
    colnames(B1x) <- b1_names
    nd[, b0_names] <- B0x
    nd[, b1_names] <- B1x
    
    Xmat  <- stats::model.matrix(trm, nd)
    eta_i <- drop(Xmat %*% beta)
    mu_i  <- fam$linkinv(eta_i)
    list(mu_i = mu_i, eta_i = eta_i, Xmat = Xmat)
  }
  
  # delta CI for g(E[μ]) at x for group t
  ci_g_link <- function(w, prov) {
    mu_i  <- prov$mu_i
    eta_i <- prov$eta_i
    Xmat  <- prov$Xmat
    
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
  
  # evaluate for Control (t=0) and Treated (t=1)
  eval_group <- function(tval, label) {
    one_x <- function(x) {
      w <- Kg((Xobs - x)/h); w <- w / sum(w)
      prov <- mu_provider(x, tval)
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


# -- helper: gradient of m_t(x) w.r.t. beta at a single x and t in {0,1}
.grad_at <- function(fit_out, x, tval) {
  fit <- fit_out$fit; dat <- fit_out$data
  fam <- fit$family; beta <- stats::coef(fit)
  trm <- stats::delete.response(stats::terms(fit))
  k0  <- fit_out$knots0; bk0 <- fit_out$bknots0
  k1  <- fit_out$knots1; bk1 <- fit_out$bknots1
  
  # kernel weights on full X
  Kg <- function(u) { ex <- -0.5*u^2; ex <- ex - max(ex); exp(ex) }
  Xobs <- dat$risk_last
  h    <- BANDWIDTH
  w <- Kg((Xobs - x)/h); w <- w / sum(w)
  
  # recompute bases at x and set t := tval for all rows
  n <- nrow(dat)
  nd <- dat
  nd$risk_last <- x
  nd$t <- as.numeric(tval)
  b0_names <- grep("^B0_", names(dat), value = TRUE)
  b1_names <- grep("^B1_", names(dat), value = TRUE)
  B0x <- splines::ns(rep(x, n), knots = k0, Boundary.knots = bk0)
  B1x <- splines::ns(rep(x, n), knots = k1, Boundary.knots = bk1)
  colnames(B0x) <- b0_names; nd[, b0_names] <- B0x
  colnames(B1x) <- b1_names; nd[, b1_names] <- B1x
  
  Xmat <- stats::model.matrix(trm, nd)
  eta  <- as.numeric(Xmat %*% beta)
  dmu  <- fit$family$mu.eta(eta)
  drop(t(Xmat) %*% (w * dmu))   # ∑ w_i μ'(η_i) X_i
}

# -- ATE at cutoff using conditional_curves() for means; joint delta for CI
ate_at_threshold <- function(fit_out, sim_out, ci = 0.95) {
  # cutoff and covariance
  c_star <- tail(sim_out$results$threshold_used, 1)
  Vb <- 0.5 * (stats::vcov(fit_out$fit) + t(stats::vcov(fit_out$fit)))
  z  <- stats::qnorm((1 + ci)/2)
  
  # point means from conditional_curves()
  cc <- conditional_curves(fit_out, x_seq = c_star, level = ci)
  m0 <- subset(cc, grp == "Control")$mean
  m1 <- subset(cc, grp == "Treated")$mean
  
  # joint gradient difference and SE
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
#   x_seq     : numeric grid of target risk values at which to compute E[Y(a)|R=x]
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
  Kg <- function(u) { ex <- -0.5*u^2; ex <- ex - max(ex); exp(ex) }
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

  # point means from conditional_curves()
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


ate_outreg <- function(sim_out, cohort_df, family) {
  
  # Transform back to raw avlues
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
  
  # Fit unweighted model
  df_fit <- data.frame(y = y, z = z, X)
  fit <- glm(y ~ z + ., data = df_fit, family = family)
  
  # Predict counterfactual outcomes for each person
  df_cf1 <- df_fit
  df_cf1$z <- 1
  mu1 <- predict(fit, newdata = df_cf1, type = "response")
  
  df_cf0 <- df_fit
  df_cf0$z <- 0
  mu0 <- predict(fit, newdata = df_cf0, type = "response")
  
  # Weighted average difference at c_star
  est <- weighted.mean(mu1 - mu0, wts)
  
  # TO-DO: Add confidence intervals
  ci <- c(lo = NaN, hi = NaN)
  
  list(ate = est, ci = ci)
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
  
  # IPW estimating equation
  ipw_ate <- function(y, z, ps, wts) {
    treated   <- (z == 1)
    control   <- (z == 0)
    
 # IPW estimates of E[Y(1)] and E[Y(0)]
  mu1 <- sum(wts * y * treated / ps) / sum(wts * treated / ps)
  mu0 <- sum(wts * y * control / (1 - ps)) / sum(wts * control / (1 - ps))
    
    mu1 - mu0
  }
  
  est <- ipw_ate(y, z, ps, wts)
  
  list(ate = est, ci = c(lo = NaN, hi = NaN))  # Placeholder CI
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
  
  # Estimate propensity scores
  df_ps <- data.frame(z = z, X)
  ps_model <- glm(z ~ ., data = df_ps, family = binomial())
  ps <- predict(ps_model, type = "response")
  
  # Estimate outcome regression
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
  
  # Efficient influence function (EIF)
  eif <- ((treated * (y - mu1)) / ps + mu1) -
    ((control * (y - mu0)) / (1 - ps) + mu0)
  
  # AIPW estimate
  est <- weighted.mean(eif, wts)
  
  # Placeholder CI
  list(ate = est, ci = c(lo = NaN, hi = NaN))
}