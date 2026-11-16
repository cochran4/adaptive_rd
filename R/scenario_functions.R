## ---------------------Threshold adaptation-------------------------

# No adaptation; the caller supplies the fixed threshold.
threshold_method_none <- list(type = "none")

# Start with no assignments, then ramp the capacity target to 30%.
threshold_method_quantile <- list(
  type                  = "quantile",
  initial_threshold     = Inf,
  target_treat_rate     = 0.30,
  blocks_to_target_rate = 6
)

# Benefit adaptation: target an additional 30 mg/dL reduction in LDL
threshold_method_benefit <- list(
  type              = "benefit",
  initial_threshold = 0.10,
  target_reduction  = 30,
  step_size         = 0.01,
  confidence_level  = 2 * pnorm(1) - 1  # Estimate plus/minus one SE.
)

## ---------------------Model adaptation--------------------------

# No adaptation (fixed model)
model_method_none <- list(type = "none")

# Recalibrate original PCE risks; penalize intercept/slope changes.
model_method_recalibrate <- list(
  type                  = "recalibrate",
  coef_original         = pce_coefs,
  recalibration_penalty = 4
)

# Start with White PCE equations for each sex, then revise without race.
model_method_revise <- list(
  type             = "revise",
  coef_original    = pce_coefs,
  coef_initial     = white_pce_for_all(pce_coefs),
  revision_penalty = 3000
)

## ---------------------Outcomes--------------------------------

# ---- Outcome 1: Attendance at prevention program ----

# Coefficients approximate Gidlow et al. (2019), Figure 2a, not a trial fit.
attend_params <- list(
  intercept         = -0.58,
  risk_slope        = -3,
  assignment_effect = 0.78
)

# Logistic approximation to Gidlow et al. (2019), Figure 2a:
# standard letters versus telephone invitations.
attendance_prob <- function(
    cohort_df, idx, risk_baseline, treat, params = attend_params
) {
  stats::plogis(params$intercept + params$risk_slope * risk_baseline +
                  params$assignment_effect * treat)
}

# Draw attendance using the probability under the assigned invitation.
outcome_attend <- function(
    cohort_df, idx, risk_baseline, treat, params = attend_params
) {
  p <- attendance_prob(cohort_df, idx, risk_baseline, treat, params)
  stats::rbinom(n = length(idx), size = 1, prob = p)
}

# ---- Outcome 2: Change in LDL cholesterol ----

# Risk/LDL pairings are illustrative. Percent changes come from pooled
# CLEAR Harmony/Wisdom results (EMA Nustendi report, Figure 16).
# The residual SD approximates extension-study change SDs (Table 23).
cholest_params <- list(
  risk_anchors       = c(0, 0.15, 0.30),
  ldl_anchors        = c(80, 120, 160),
  control_percent    = c(6.137, -0.803, -5.766),
  medication_percent = c(-12.486, -18.435, -25.219),
  residual_sd        = 25
)

# Natural splines interpolate the trial-informed means at illustrative risks.
# Beyond the boundary risks, the natural splines continue linearly.
cholest_mean <- function(
    cohort_df, idx, risk_baseline, treat, params = cholest_params
) {
  basis <- function(r) cbind(1, splines::ns(
    r, knots = params$risk_anchors[2],
    Boundary.knots = params$risk_anchors[c(1, 3)]
  ))
  # Convert percentage changes to mg/dL and solve for spline coefficients.
  control <- solve(
    basis(params$risk_anchors),
    params$ldl_anchors * params$control_percent / 100
  )
  medication <- solve(
    basis(params$risk_anchors),
    params$ldl_anchors * params$medication_percent / 100
  )

  # Evaluate each patient's mean under the assigned arm.
  b <- basis(risk_baseline)
  as.numeric((1 - treat) * (b %*% control) + treat * (b %*% medication))
}

# Add independent normal variation around the assigned-arm LDL mean.
outcome_cholest <- function(
    cohort_df, idx, risk_baseline, treat, params = cholest_params
) {
  mu <- cholest_mean(cohort_df, idx, risk_baseline, treat, params)
  stats::rnorm(length(idx), mu, params$residual_sd)
}

# ---- Outcome 3: Cardiovascular event ----

# CLEAR Outcomes primary-prevention analysis: three-component endpoint.
# Use the displayed 54-month horizon and reported hazard ratio of 0.64.
cvd_params <- list(
  follow_up_months  = 54,
  reference_months  = 120,
  assignment_hr     = 0.64
)

# ---- Simulated 54-month cardiovascular-event probability ----
cvd_prob <- function(
    cohort_df, idx, risk_baseline, treat, params = cvd_params
) {
  # Avoid infinite values in the complementary log-log transformation.
  risk_baseline <- pmin(
    pmax(risk_baseline, .Machine$double.eps),
    1 - .Machine$double.eps
  )

  # Assume a constant baseline hazard to shorten the horizon; assignment
  # multiplies that hazard by assignment_hr.
  eta <- log(params$follow_up_months / params$reference_months) +
    log(-log1p(-risk_baseline)) +
    log(params$assignment_hr) * treat

  1 - exp(-exp(eta))
}

# Draw events using the probability under the assigned intervention.
outcome_cvd <- function(
    cohort_df, idx, risk_baseline, treat, params = cvd_params
) {
  p <- cvd_prob(cohort_df, idx, risk_baseline, treat, params)
  stats::rbinom(n = length(idx), size = 1, prob = p)
}
