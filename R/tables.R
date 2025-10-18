# tables.R

# Control variables 
DIGITS = 3 # Number of display digits

#-------------------------------------------------------------------------------
# compare_ate_at_threshold
#-------------------------------------------------------------------------------
# Creates a comparison table of estimated vs. true ATE at the last threshold.
#
# Arguments:
#   fit_out   : object returned by model fitting
#   sim_out   : object returned by simulation
#   cohort_df : data frame (cohort)
#   risk_fn   : function to compute risks from cohort_df
#   mean_fn   : function to compute mean outcome (truth)
#
# Returns:
#   A knitr::kable table
#-------------------------------------------------------------------------------

compare_ate_at_threshold <- function(fit_out, sim_out, cohort_df, risk_fn, mean_fn) {
  
  # Estimated local ATE
  est <- ate_at_threshold(fit_out, sim_out)
  
  # Actual local ATE (based on DGM)
  tru <- ate_at_threshold_truth(sim_out, cohort_df, mean_fn)
  
  # Recover values
  est_ate <- est$ate
  ci_lo   <- est$ci["lo"]
  ci_hi   <- est$ci["hi"]
  tru_ate <- tru$ate
  
  # Comparison
  bias    <- est_ate - tru_ate
  covered <- if (tru_ate >= ci_lo && tru_ate <= ci_hi) "Yes" else "No"
  
  tab <- data.frame(
    Metric = c("Estimated ATE", "95% CI", "Truth", "Bias", "Coverage"),
    Value  = c(
      sprintf(paste0("%.", DIGITS, "f"), est_ate),
      sprintf(paste0("[%.", DIGITS, "f, %.", DIGITS, "f]"), ci_lo, ci_hi),
      sprintf(paste0("%.", DIGITS, "f"), tru_ate),
      sprintf(paste0("%.", DIGITS, "f"), bias),
      covered
    ),
    stringsAsFactors = FALSE
  )
  
  knitr::kable(tab, align = c("l", "c"))
}