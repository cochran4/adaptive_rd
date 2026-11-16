# R/adapt_threshold.R
# Utilities to adapt a binary treatment threshold from continuous risk.

# ============================================================
# Quantile-based threshold
# ============================================================

#' Adapt a risk threshold based on a desired treatment rate
#'
#' Finds the risk threshold c such that a specified proportion
#' (`target_treat_rate`) of individuals have risk >= c.
#'
#' @param risks Risk predictions for all patients observed so far.
#' @param target_treat_rate Target assignment rate for the next block.
#'
#' @return Numeric risk threshold value.
adapt_threshold_quantile <- function(
    risks,
    target_treat_rate = 0.50
) {
  if (target_treat_rate <= 0) return(Inf)
  if (target_treat_rate >= 1) return(-Inf)

  # Threshold so that target_treat_rate fraction are >= threshold
  threshold <- stats::quantile(
    risks,
    probs = 1 - target_treat_rate,  # higher treat rate => lower threshold
    na.rm = TRUE
  )
  as.numeric(threshold)
}

# Target assignment rate for the next block after update_number updates.
quantile_target_rate <- function(method, update_number) {
  method$target_treat_rate *
    min(update_number / method$blocks_to_target_rate, 1)
}

# ============================================================
# Benefit-based threshold: local increase/decrease/hold rule
# ============================================================
# Assumes that assignment benefit increases with risk. Adapt the local
# target-seeking principle of Ivanova and Kim (2009), Biometrics 65:307-315,
# DOI: 10.1111/j.1541-0420.2008.01045.x, to the estimated assignment contrast.
adapt_threshold_benefit <- function(sim_out,
                                    target_reduction,
                                    step_size = 0.01,
                                    confidence_level = 2 * pnorm(1) - 1) {
  # Use the proposed estimator and its joint interval at the current cutoff.
  fit_out <- estimate_spline(sim_out, family = gaussian())
  estimate <- ate_at_threshold(fit_out, sim_out, ci = confidence_level)
  current_threshold <- estimate$cutoff

  # The estimator reports assigned minus unassigned LDL change. Reverse the
  # signs and endpoints to express the interval as an additional LDL reduction.
  reduction_lower <- -estimate$ci["hi"]
  reduction_upper <- -estimate$ci["lo"]

  # Larger benefits require a higher risk threshold under this scenario's
  # monotonicity assumption. Hold when the interval includes the target.
  if (reduction_upper < target_reduction) {
    current_threshold <- current_threshold + step_size
  } else if (reduction_lower > target_reduction) {
    current_threshold <- current_threshold - step_size
  }

  as.numeric(pmin(pmax(current_threshold, 0), 1))
}

# ============================================================
# R/adapt_threshold.R
# Lightweight wrapper to choose and apply a threshold adaptation method.
# ============================================================

# Expects:
# - method$type: "none", "quantile", or "benefit"
# - method$initial_threshold (used when type == "none")
# - method$target_treat_rate and method$blocks_to_target_rate (for "quantile")
# - method$initial_threshold, method$target_reduction, method$step_size, and
#   method$confidence_level (for "benefit")

adapt_threshold <- function(method, risks, sim_out, update_number) {

  if (identical(method$type, "none")) {
    return(method$initial_threshold)
  }

  if (identical(method$type, "quantile")) {
    current_target_rate <- quantile_target_rate(method, update_number)
    return(adapt_threshold_quantile(
      risks             = risks,
      target_treat_rate = current_target_rate
    ))
  }

  if (identical(method$type, "benefit")) {
    return(adapt_threshold_benefit(
      sim_out          = sim_out,
      target_reduction = method$target_reduction,
      step_size        = method$step_size,
      confidence_level = method$confidence_level
    ))
  }

  # If type is unrecognized, just keep current threshold.
  method$initial_threshold
}
