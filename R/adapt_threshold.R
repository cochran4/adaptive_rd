# R/adapt_threshold.R
# Utilities to adapt a binary treatment threshold from continuous risk.

# ============================================================
# Quantile-based threshold
# ============================================================

#' Adapt a risk threshold based on a desired treatment rate
#'
#' Finds the risk threshold c such that a specified proportion
#' (`desired_treat_rate`) of individuals have risk >= c.
#'
#' @param data Data frame containing the analytic dataset.
#' @param rows Optional integer vector of row indices to use. If NULL,
#'             all rows in `data` are used.
#' @param risk_fn Function taking (data, coef_adjust) and returning a numeric
#'                vector of risk estimates in [0, 1].
#' @param coef_adjust Optional object passed to `risk_fn` (e.g., PCE overrides or
#'                    LP rescaling controls). Only needed if `risk_fn` uses it.
#' @param desired_treat_rate Target fraction of individuals to treat
#'                           (e.g., 0.50 for 50% treated).
#'
#' @return Numeric risk threshold value.
adapt_threshold_quantile <- function(
    risks,
    desired_treat_rate = 0.50
) {
  
  # threshold so that desired_treat_rate fraction are >= threshold
  threshold <- stats::quantile(
    risks,
    probs = 1 - desired_treat_rate,  # higher treat rate => lower threshold
    na.rm = TRUE
  )
  as.numeric(threshold)
}

# ============================================================
# Number-needed-to-treat (NNT)-based threshold (pointwise)
# ============================================================
#
# Chooses the risk threshold c such that the pointwise treatment effect
# at risk = c yields an NNT closest to `desired_nnt`. Continuous outcomes only.
# Uses:
#   - estimate_spline(sim_out, "gaussian")
#   - conditional_curves(fit_out, x_seq = risk_grid, level)
#   - pooled SD from sim_out$results$outcome by treat
#
# Returns: list with
#   - threshold : chosen risk threshold
#   - grid      : data.frame(risk_threshold, effect, d, nnt)
#
adapt_threshold_nnt <- function(sim_out,
                                risks,
                                desired_nnt,
                                min_treat_rate = 0,
                                max_treat_rate = 1, 
                                shrinkage = .5) {
  
  # Previous threshold
  prev_threshold <- tail(sim_out$results$threshold_used, 1)
  
  # Compute target d from desired NNT
  target_d <- sqrt(2) * qnorm(0.5 * (1 + 1 / desired_nnt))
  
  # Minimum risk based on max possible treat rate
  min_risk <- adapt_threshold_quantile(risks, max_treat_rate)

  # Maximum risk based on min possible treat rate
  max_risk <- adapt_threshold_quantile(risks, min_treat_rate)
  
  # Possible risk thresholds
  risk_grid = ((min_risk*1000):(max_risk*1000))/1000
  
  #---------------
  # Get estimates
  #---------------
  
  # Only works for Gaussian
  family <- gaussian()

  # Fit spline to simulation output
  fit_out <- estimate_spline(sim_out, family)
  
  # curves: long DF with columns risk_last, mean, grp ∈ {Control, Treated}
  curves <- conditional_curves(fit_out, risk_grid)
  control <- curves[curves$grp == "Control", ]
  treated <- curves[curves$grp == "Treated", ]
  
  # pointwise effect at each risk_last (assumes identical order)
  effect <- treated$mean - control$mean
  
  # Pooled standard deviations (assumes continuous outcome)
  df <- sim_out$results
  y0 <- df$outcome[df$treat == 0]
  y1 <- df$outcome[df$treat == 1]
  n0 <- length(y0); n1 <- length(y1)
  s0 <- sd(y0, na.rm = TRUE); s1 <- sd(y1, na.rm = TRUE)
  sp <- sqrt(((n0 - 1) * s0^2 + (n1 - 1) * s1^2) / (n0 + n1 - 2))
  
  # Cohen's d for each LATE
  cohen_d <- - effect / sp

  # Pick risk threshold where d is closest to target_d ---
  opt_idx <- which.min(abs(cohen_d - target_d))
  
  # New threshold
  new_threshold <- risk_grid[opt_idx]
  
  # Return weighted average of two thresholds
  return( shrinkage*prev_threshold + (1-shrinkage)*new_threshold )

}

# ============================================================
# R/adapt_threshold.R
# Lightweight wrapper to choose and apply a threshold adaptation method.
# ============================================================

# Expects:
# - method$type: "none", "quantile", or "nnt"
# - method$initial_threshold (used when type == "none")
# - method$desired_treat_rate (for "quantile")
# - method$treatment_effect_fn, method$desired_nnt (for "nnt")
# - optional: method$min_treat_rate, method$max_treat_rate (for "nnt")

adapt_threshold <- function(method, risks, sim_out) {
  
  if (identical(method$type, "none")) {
    return(method$initial_threshold)
  }
  
  if (identical(method$type, "quantile")) {
    return(adapt_threshold_quantile(
      risks              = risks,
      desired_treat_rate = method$desired_treat_rate
    ))
  }
  
  if (identical(method$type, "nnt")) {
    return(adapt_threshold_nnt(sim_out,
                               risks,
                               method$desired_nnt))
  }
  
  # If type is unrecognized, just keep current threshold.
  method$initial_threshold
}

