## ------------------------------------------------------------------
##                             SETUP
## ------------------------------------------------------------------

## ---------------------Initialize-----------------------------------

# Core package imports
source("R/imports.R")           # attaches dplyr, tibble, mice, etc.

# Setting seed for reproducibility
set.seed(1234)

## ---------------------Build datasets--------------------------------

data_dir     <- 'data_raw'
nh_list      <- read_nhanes_2017_2018(data_dir)
cohort_bundle <- build_cohort(nh_list, impute = TRUE, summarize = TRUE)
cohort_df     <- cohort_bundle$data
rm(list = c("data_dir", "nh_list"))



## ---------------------Scenarios--------------------------------

sim_context <- list(
  simulate_design = simulate_design,
  cohort_df     = cohort_df,
  pce_predict   = pce_predict,
  threshold_method_quantile = threshold_method_quantile,
  threshold_method_nnt      = threshold_method_nnt,
  model_method_none         = model_method_none,
  model_method_revise       = model_method_revise,
  model_method_recalibrate  = model_method_recalibrate,
  outcome_attend   = outcome_attend,
  outcome_cholest  = outcome_cholest,
  outcome_cvd      = outcome_cvd,
  attendance_prob  = attendance_prob,
  cholest_mean     = cholest_mean,
  cvd_prob         = cvd_prob
)

scenarios <- list(
  S1 = function(ctx) list(
    sim_out = ctx$simulate_design(
      ctx$cohort_df, ctx$pce_predict, 400, 100, 26,
      ctx$threshold_method_quantile, ctx$model_method_none, ctx$outcome_attend
    ),
    family  = binomial(link = "logit"),
    mean_fn = ctx$attendance_prob
  ),
  
  S2 = function(ctx) list(
    sim_out = ctx$simulate_design(
      ctx$cohort_df, ctx$pce_predict, 400, 100, 26,
      ctx$threshold_method_quantile, ctx$model_method_none, ctx$outcome_cholest
    ),
    family  = gaussian(),
    mean_fn = ctx$cholest_mean
  ),
  
  S3 = function(ctx) list(
    sim_out = ctx$simulate_design(
      ctx$cohort_df, ctx$pce_predict, 400, 100, 26,
      ctx$threshold_method_nnt, ctx$model_method_none, ctx$outcome_cholest
    ),
    family  = gaussian(),
    mean_fn = ctx$cholest_mean
  ),
  
  S4 = function(ctx) list(
    sim_out = ctx$simulate_design(
      ctx$cohort_df, ctx$pce_predict, 400, 100, 26,
      ctx$threshold_method_quantile, ctx$model_method_recalibrate, ctx$outcome_cvd
    ),
    family  = binomial(link = "logit"),
    mean_fn = ctx$cvd_prob
  ),
  
  S5 = function(ctx) list(
    sim_out = ctx$simulate_design(
      ctx$cohort_df, ctx$pce_predict, 400, 100, 26,
      ctx$threshold_method_quantile, ctx$model_method_revise, ctx$outcome_cvd
    ),
    family  = binomial(link = "logit"),
    mean_fn = ctx$cvd_prob
  )
)


## ------------------------------------------------------------------
##                      Monte Carlo Simulation
## ------------------------------------------------------------------

# One run of the simulation for a single scenario
run_once <- function(scenario_fn, ctx) {
  
  source("R/imports.R") 
  source("R/scenario_functions.R") 
  
  out     <- scenario_fn(ctx)
  sim_out <- out$sim_out
  fam     <- out$family
  mean_fn <- out$mean_fn
  
  # Estimation
  fit_out     <- suppressWarnings(estimate_spline(sim_out, family = fam))
  est_b       <- ate_at_threshold(fit_out, sim_out)
  
  # Ground truth ATE at decision threshold
  tru_b      <- ate_at_threshold_truth(sim_out, ctx$cohort_df, mean_fn)
  
  # Estimators for comparison
  est_naive  <- ate_naive(sim_out, family = fam)        # Unadjusted difference in means
  est_outreg <- ate_outreg(sim_out, cohort_df, family = fam)  # Outcome regression
  est_ipw    <- ate_ipw(sim_out, cohort_df, family = fam)     # Inverse probability weighting
  est_aipw   <- ate_aipw(sim_out, cohort_df, family = fam)    # Augmented IPW (doubly robust)
  
        
  # Return estimates
  list(
    truth = list(
      ate = tru_b$ate
    ),
    spline = list(
      ate = est_b$ate,
      ci  = est_b$ci
    ),
    naive = list(
      ate = est_naive$ate,
      ci  = est_naive$ci
    ),
    outreg = list(
      ate = est_outreg$ate,
      ci  = est_outreg$ci
    ),
    ipw   = list(
      ate = est_ipw$ate,
      ci  = est_ipw$ci
    ),
    aipw   = list(
      ate = est_aipw$ate,
      ci  = est_aipw$ci
    )
  )
}


run_mc_ate <- function(scenarios, sim_context, B = 200) {
  out <- vector("list", length(scenarios))
  names(out) <- names(scenarios)
  
  future::plan("multisession")
  handlers(global = TRUE)
  
  for (nm in names(scenarios)) {
    
    progressr::with_progress({
      results <- future_lapply(
        seq_len(B),
        FUN = function(b, scenario_fn, ctx) run_once(scenario_fn, ctx),
        scenario_fn = scenarios[[nm]],
        ctx = sim_context,
        future.seed = TRUE
      )
    })
    
    methods <- c("spline", "ipw", "outreg", "naive", "aipw")  # extend if needed
    
    df_list <- lapply(methods, function(m) {
      est  <- vapply(results, function(r) r[[m]]$ate, numeric(1))
      tru  <- vapply(results, function(r) r$truth$ate, numeric(1))
      
      df <- data.frame(
        method   = m,
        scenario = nm,
        est      = est,
        truth    = tru
      )
      
      if (m == "spline") {
        df$ci_lo <- vapply(results, function(r) r[[m]]$ci["lo"], numeric(1))
        df$ci_hi <- vapply(results, function(r) r[[m]]$ci["hi"], numeric(1))
      } else {
        df$ci_lo <- NA
        df$ci_hi <- NA
      }
      
      df
    })
    
    out[[nm]] <- do.call(rbind, df_list)
  }
  
  do.call(rbind, out)
}

mc <- run_mc_ate(
  scenarios   = scenarios,
  sim_context = sim_context,
  B           = 2000
)
mc
saveRDS(mc, file = "mc.rds")

plot_comparison(
  mc,
  method_labels = c("spline" = "Adaptive RD", "naive" = "Naive", "outreg" = "Outcome Regression", "ipw" = "Inverse Probability Weighting",  "aipw" = "Augmented Inverse Probability Weighting"),
  scenario_labels = c(
    "S1" = "Scenario 1",
    "S2" = "Scenario 2",
    "S3" = "Scenario 3",
    "S4" = "Scenario 4",
    "S5" = "Scenario 5"
  ),
  save_path = "figures/comparison.rds",
  save_path2 = "figures/fig_comparison.pdf"
)