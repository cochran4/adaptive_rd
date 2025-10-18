## ---------------------Threshold adaptation-------------------------

# Quantile-based adaptation: set threshold so that a fraction is treated
threshold_method_quantile <- list(
  type               = "quantile",
  initial_threshold  = 0.10, # initial_threshold,
  desired_treat_rate = 0.30  # 30% referral
)

# NNT adaptation: set threshold so that we target a specific NNT
target_nnt <- 3 # 3 persons
threshold_method_nnt <- list(
  type               = "nnt",
  initial_threshold  = 0.10,       # initial_threshold,
  desired_nnt        = target_nnt   # Target NNT
) 

## ---------------------Model adaptation--------------------------

# No adaptation (fixed model)
model_method_none <- list(type = "none")

# Recalibrate model
model_method_recalibrate <- list(
  type              = "recalibrate",
  coef_original     = pce_coefs,
  coef_adjusted     = NULL,
  shrinkage_factor  = 3000
)

# Revise model
model_method_revise <- list(
  type               = "revise",
  coef_original      = pce_coefs,
  coef_adjusted      = NULL,
  shrinkage_factor   = 3000
)

## ---------------------Outcomes--------------------------------

# ---- Outcome 1: Attendance to prevention program ----

# Parameters
attend_params <- list(
  treatment_effect_intercept = 0.25,
  reference_intercept        = 0.00,
  reference_slope            = 0.10
)

# ---- Attendance probability ----
attendance_prob <- function(cohort_df, idx, risk_baseline, treat, params = attend_params) {
  
  # Probability formula
  p <- params$reference_intercept +
    params$reference_slope * risk_baseline +
    treat * params$treatment_effect_intercept*(risk_baseline+1/2)*(1-risk_baseline)*16/9
  
}

outcome_attend <- function(cohort_df, idx, risk_baseline, treat, params = attend_params) {
  p <- attendance_prob(cohort_df, idx, risk_baseline, treat, params)
  stats::rbinom(n = length(idx), size = 1, prob = p)
}


# ---- Outcome 2: Change in cholesterol ----


# Parameters
cholest_params <- list(
  treatment_effect_intercept =   0.00,  # shift for treated
  treatment_effect_slope     = -10.00,  # effect scales with risk
  reference_intercept        =  2.00,   # baseline drift
  reference_slope            =  0.00,   # dependence on baseline_sbp
  within_sd                  =  5.00    # residual SD
)

# Mean model: expected change in cholesterol
cholest_mean <- function(cohort_df, idx, risk_baseline, treat, params = cholest_params) {
  
  mu <- params$reference_intercept +
    params$reference_slope * risk_baseline +
    treat * (params$treatment_effect_intercept +
               params$treatment_effect_slope * risk_baseline)
  mu
}

# Outcome generator: draws cholesterol change given the mean model
outcome_cholest <- function(cohort_df, idx, risk_baseline, treat, params = cholest_params) {
  mu <- cholest_mean(cohort_df, idx, risk_baseline, treat, params = params)
  stats::rnorm(n = length(idx), mean = mu, sd = params$within_sd)
}


# ---- Outcome 3: Cardiovascular event ----

# Parameters
cvd_params <- list(
  treatment_effect_intercept =  0.4,    # shift for treated
  baseline_intercept         =  0.1,    # baseline drift
  baseline_slope             =  0.9     # dependence on baseline_sbp
)

# ---- Simulated probability (pce model is inaccurate) ----
cvd_prob <- function(cohort_df, idx, risk_baseline, treat, params = cvd_params) {
  
  # Copy df for local edits
  df <- cohort_df[idx, ]
  
  # Recompute baseline risk with updated labels
  risk_new <- pce_predict(df, coefs = pce_coefs, coef_adjust = NULL)
  
  # Get on linear predictor scale
  lp_new <- log(-log(1-risk_new))
  
  # Apply affine transformation
  lp_new <- params$baseline_intercept + 
    params$baseline_slope * lp_new + 
    params$treatment_effect_intercept * treat
  
  # Transform back to response scale
  p <- 1 - exp(-exp(lp_new))
  
  # Return probability
  return(p)
}

outcome_cvd <- function(cohort_df, idx, risk_baseline, treat, params = cvd_params) {
  p <- cvd_prob(cohort_df, idx, risk_baseline, treat, params)
  stats::rbinom(n = length(idx), size = 1, prob = p)
}
