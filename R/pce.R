# R/pce.R
# Compute 10-year ASCVD risk using the Pooled Cohort Equations (PCE).
# Input data columns expected:
#   person_id, age, sex, race, sbp, smoker_current, diabetes, bp_treated, hdl_c, tot_chol

# --------------------------------------------------------------------
# Compute a single row's linear predictor given that row's terms and a
# coefficient entry (from pce_coefs or overrides).
#
# coef_entry must contain:
#   $betas  : named numeric vector (names must be present in terms_row)
#   $s0     : baseline survival at 10 years (not used here, but upstream)
#   $mean_x : the mean LP used in derivation (not used here directly)
# --------------------------------------------------------------------
.lp_row <- function(terms_row, coef_entry) {

  # Coefficients
  betas  <- coef_entry$betas
  
  # Identify which betas correspond to columns we actually built
  common <- intersect(names(betas), names(terms_row))
  
  # Weighted sum = linear predictor contribution
  sum(as.numeric(betas[common])*as.numeric(terms_row[common]))
}

# --------------------------------------------------------------------
# Select the active coefficient bundle:
#   - use `coef_adjust$coefs` if provided,
#   - otherwise fall back to the default `coefs` (usually `pce_coefs`).
# --------------------------------------------------------------------
.get_active_coefs <- function(default_coefs, coef_adjust) {
  if (!is.null(coef_adjust)) return(coef_adjust)
  default_coefs
}


# --------------------------------------------------------------------
# Compute 10-year ASCVD risk using PCE, with optional coefficient adjustments
# --------------------------------------------------------------------
pce_predict <- function(df, coefs = pce_coefs, coef_adjust = NULL) {
  
  # Select the coefficient bundle (override if provided)
  active_coefs <- .get_active_coefs(coefs, coef_adjust)
  
  # Preallocate output vectors for efficiency
  n  <- nrow(df)
  lp <- numeric(n)  # linear predictor (before optional rescale)
  s0 <- numeric(n)  # baseline survival at 10 years (group-specific)
  mean_x <- numeric(n)  # mean LP used in derivation (group-specific)

  # Compute per-group LP and pull s0/mean_x
  groups_unique <- unique(df$group)
  for (g in groups_unique) {
    idx <- which(df$group == g)  # rows in this group
    ce  <- active_coefs[[g]]        # coefficient entry for this group
    
    # Row-wise LP = sum(betas * terms) for each person in the group
    lp[idx] <- apply(df[idx, , drop = FALSE], 1L, .lp_row, coef_entry = ce)
    
    # Baseline survival and mean LP are group constants pulled from coef bundle
    s0[idx]     <- ce$s0
    mean_x[idx] <- ce$mean_x
  }

  # Convert LP to 10-year risk per the PCE formula:
  risk10 <- 1 - (s0 ^ exp(lp - mean_x))
  
  # Numerical guard: clamp to (0, 1)
  risk10 <- pmin(pmax(risk10, 1e-6), 1 - 1e-6)
  
  # Return risk prediction
  risk10
}


# --------------------------------------------------------------------
# Coefficients for the ACC/AHA Pooled Cohort Equations (PCE).
# Groups: White_F (W_f), Black_F (AF_f), White_M (W_m), Black_M (AF_m)
# --------------------------------------------------------------------
pce_coefs <- list(
  White_F = list(
    betas = c(
      ln_age                  = -29.799,
      ln_age_sqrd             =   4.884,
      ln_total_cholest        =  13.540,   # Corrected: Was 3.540 in All of Us Paper
      ln_age_totcholest       =  -3.114,
      ln_hdlC                 = -13.578,
      ln_age_hdlC             =   3.149,
      ln_treated_BP           =   2.019,
      ln_age_BP               =   0.000,   # no term for WF per your table
      ln_untreated_BP         =   1.957,
      ln_age_ln_untreated_BP  =   0.000,   # no term for WF per your table
      smoker                  =   7.574,
      ln_age_smoker           =  -1.665,
      diabetes                =   0.661
    ),
    s0     = 0.9665,
    mean_x = -29.18
  ),
  Black_F = list(
    betas = c(
      ln_age                  =  17.114,
      ln_age_sqrd             =   0.000,
      ln_total_cholest        =   0.940,
      ln_age_totcholest       =   0.000,
      ln_hdlC                 = -18.920,
      ln_age_hdlC             =   4.475,
      ln_treated_BP           =  29.291,
      ln_age_BP               =  -6.432,   # this is ln_age * ln_sbp_treated
      ln_untreated_BP         =  27.820,
      ln_age_ln_untreated_BP  =  -6.087,   # this is ln_age * ln_sbp_untreated
      smoker                  =   0.691,
      ln_age_smoker           =   0.000,
      diabetes                =   0.874
    ),
    s0     = 0.9533,
    mean_x = 86.61
  ),
  White_M = list(
    betas = c(
      ln_age                  =  12.344,
      ln_age_sqrd             =   0.000,
      ln_total_cholest        =  11.853,
      ln_age_totcholest       =  -2.664,
      ln_hdlC                 =  -7.990,
      ln_age_hdlC             =   1.769,
      ln_treated_BP           =   1.797,
      ln_age_BP               =   0.000,
      ln_untreated_BP         =   1.764,
      ln_age_ln_untreated_BP  =   0.000,
      smoker                  =   7.837,
      ln_age_smoker           =  -1.795,
      diabetes                =   0.658
    ),
    s0     = 0.9144,
    mean_x = 61.18
  ),
  Black_M = list(
    betas = c(
      ln_age                  =   2.469,
      ln_age_sqrd             =   0.000,
      ln_total_cholest        =   0.302,
      ln_age_totcholest       =   0.000,
      ln_hdlC                 =  -0.307,
      ln_age_hdlC             =   0.000,
      ln_treated_BP           =   1.916,
      ln_age_BP               =   0.000,
      ln_untreated_BP         =   1.809,
      ln_age_ln_untreated_BP  =   0.000,
      smoker                  =   0.549,
      ln_age_smoker           =   0.000,
      diabetes                =   0.645
    ),
    s0     = 0.8954,
    mean_x = 19.54
  )
)

# Use White coefficients for individuals labeled as "Other"
pce_coefs$Other_F <- pce_coefs$White_F
pce_coefs$Other_M <- pce_coefs$White_M

