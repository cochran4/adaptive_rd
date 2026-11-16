## ------------------------------------------------------------------
##                             SETUP
## ------------------------------------------------------------------

## ---------------------Initialize-----------------------------------

# Load packages, shared functions, and scenario settings.
source("R/imports.R")

# Set the initial random seed; cohort imputation also uses its own fixed seed.
set.seed(1234)

## ---------------------Build datasets--------------------------------

data_dir      <- "data_raw"
nh_list       <- read_nhanes_2017_2018(data_dir)
cohort_bundle <- build_cohort(nh_list, impute = TRUE, summarize = TRUE)
cohort_df      <- cohort_bundle$data
rm(list = c("data_dir", "nh_list"))

## ---------------------Scenarios--------------------------------

# Use the selected adaptation settings for the final estimator comparison.
comparison_benefit <- threshold_method_benefit
comparison_revision <- model_method_revise
comparison_revision$revision_penalty <- 3000

sim_context <- list(
  simulate_design = simulate_design,
  cohort_df     = cohort_df,
  pce_predict   = pce_predict,
  threshold_method_quantile = threshold_method_quantile,
  threshold_method_benefit  = comparison_benefit,
  model_method_none         = model_method_none,
  model_method_revise       = comparison_revision,
  model_method_recalibrate  = model_method_recalibrate,
  outcome_attend   = outcome_attend,
  outcome_cholest  = outcome_cholest,
  outcome_cvd      = outcome_cvd,
  attendance_prob  = attendance_prob,
  cholest_mean     = cholest_mean,
  cvd_prob         = cvd_prob
)

# Each scenario has 30 blocks of 100 patients and returns its outcome family
# and known mean function alongside the simulated data.
scenarios <- list(
  S1 = function(ctx) list(
    sim_out = ctx$simulate_design(
      ctx$cohort_df, ctx$pce_predict, 100, 30,
      ctx$threshold_method_quantile, ctx$model_method_none, ctx$outcome_attend
    ),
    family  = binomial(link = "logit"),
    mean_fn = ctx$attendance_prob
  ),

  S2 = function(ctx) list(
    sim_out = ctx$simulate_design(
      ctx$cohort_df, ctx$pce_predict, 100, 30,
      ctx$threshold_method_benefit, ctx$model_method_none, ctx$outcome_cholest
    ),
    family  = gaussian(),
    mean_fn = ctx$cholest_mean
  ),

  S3 = function(ctx) list(
    sim_out = ctx$simulate_design(
      ctx$cohort_df, ctx$pce_predict, 100, 30,
      ctx$threshold_method_quantile, ctx$model_method_recalibrate,
      ctx$outcome_cvd
    ),
    family  = binomial(link = "logit"),
    mean_fn = ctx$cvd_prob
  ),

  S4 = function(ctx) list(
    sim_out = ctx$simulate_design(
      ctx$cohort_df, ctx$pce_predict, 100, 30,
      ctx$threshold_method_quantile, ctx$model_method_revise, ctx$outcome_cvd
    ),
    family  = binomial(link = "logit"),
    mean_fn = ctx$cvd_prob
  )
)

## ------------------------------------------------------------------
##                      Monte Carlo Simulation
## ------------------------------------------------------------------

# Simulate one dataset, then evaluate all estimators on that same dataset.
run_once <- function(scenario_fn, ctx) {

  # Load dependencies and shared definitions in the parallel worker.
  source("R/imports.R")
  source("R/scenario_functions.R")

  out     <- scenario_fn(ctx)
  sim_out <- out$sim_out
  fam     <- out$family
  mean_fn <- out$mean_fn

  # Fit the proposed estimator and evaluate it at the final threshold.
  fit_out     <- suppressWarnings(estimate_spline(sim_out, family = fam))
  est_b       <- ate_at_threshold(fit_out, sim_out)

  # Evaluate the known outcome model near the final decision threshold.
  # Outcomes use original race-specific PCE risk, even when revision
  # starts from a different deployed model (Scenario 4).
  original_risk <- ctx$pce_predict(ctx$cohort_df)
  truth_mean <- function(cohort_df, idx, risk_baseline, treat) {
    mean_fn(cohort_df, idx, original_risk[idx], treat)
  }
  tru_b <- ate_at_threshold_truth(sim_out, ctx$cohort_df, truth_mean)

  # Apply each comparator to the same simulated patients and assignments.
  est_naive <- ate_naive(sim_out, family = fam)
  est_outreg <- ate_outreg(sim_out, ctx$cohort_df, family = fam)
  est_outreg_separate <- ate_outreg_separate(
    sim_out, ctx$cohort_df, family = fam
  )
  est_ipw <- ate_ipw(sim_out, ctx$cohort_df, family = fam)
  est_aipw <- ate_aipw(sim_out, ctx$cohort_df, family = fam)

  # Return the shared truth and each method's estimate and available interval.
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
    outreg_separate = list(
      ate = est_outreg_separate$ate,
      ci  = est_outreg_separate$ci
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

# Run B replications per selected scenario and combine results into one table.
run_mc_ate <- function(scenarios, sim_context, B = 200) {
  out <- vector("list", length(scenarios))
  names(out) <- names(scenarios)

  # Use separate R sessions for parallel replications.
  future::plan("multisession")
  handlers(global = TRUE)

  for (nm in names(scenarios)) {

    progressr::with_progress({
      results <- future_lapply(
        seq_len(B),
        FUN = function(b, scenario_fn, ctx) run_once(scenario_fn, ctx),
        scenario_fn = scenarios[[nm]],
        ctx = sim_context,
        # Generate reproducible random-number streams for parallel workers.
        future.seed = TRUE
      )
    })

    # Store estimates from all methods; intervals are retained only for ours.
    methods <- c(
      "spline", "ipw", "outreg", "outreg_separate", "naive", "aipw"
    )

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

# Rerun selected scenarios, retaining completed results for the others.
selected_scenarios <- getOption("comparison_scenarios", names(scenarios))
results_dir <- getOption("comparison_results_dir", "results/comparison")
results_file <- file.path(results_dir, "comparison_results.rds")

mc <- run_mc_ate(
  scenarios   = scenarios[selected_scenarios],
  sim_context = sim_context,
  B           = getOption("comparison_reps", 2000L)
)
# A partial rerun replaces only the selected scenarios in the saved results.
if (length(selected_scenarios) < length(scenarios)) {
  previous <- readRDS(results_file)
  mc <- rbind(previous[!previous$scenario %in% selected_scenarios, ], mc)
  mc <- mc[order(match(mc$scenario, names(scenarios))), ]
}

# Save the combined paper results and rebuild the comparison figure.
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
saveRDS(mc, results_file)
invisible(plot_comparison_paper(mc))
