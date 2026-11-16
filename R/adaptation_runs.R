# Adaptation simulation experiments
#
# Run this file from the project root after sourcing R/imports.R. The Monte
# Carlo results are saved independently of the Quarto document so rendering
# the document never reruns the simulations.

# Settings used in the paper: threshold rules in S1 and revision penalties in
# S4.
adaptation_default_grid <- function() {
  list(
    S1 = tibble::tibble(
      setting = c("fixed", "ramp_3", "ramp_6"),
      blocks_to_target_rate = c(0, 3, 6)
    ),
    S4 = tibble::tibble(
      setting = paste0("revision_penalty_", c(1000, 3000, 10000)),
      revision_penalty = c(1000, 3000, 10000)
    )
  )
}

# Build the source cohort used for every setting; imputation has its own seed.
prepare_adaptation_cohort <- function(data_dir = "data_raw", seed = 1234) {
  set.seed(seed)
  nhanes <- read_nhanes_2017_2018(data_dir)
  build_cohort(nhanes, impute = TRUE, summarize = FALSE)
}

# Convert one row of the settings grid into simulation methods and an outcome.
adaptation_scenario_spec <- function(scenario, setting_row, cohort_df) {
  # Copy shared settings so this grid row changes only the selected feature.
  quantile_method <- threshold_method_quantile
  revise_method <- model_method_revise

  if (scenario == "S1") {
    if (identical(setting_row$setting[[1]], "fixed")) {
      # Hold the source cohort's 70th-percentile risk fixed from block 1.
      quantile_method <- threshold_method_none
      quantile_method$initial_threshold <-
        unname(stats::quantile(pce_predict(cohort_df), 0.70))
    } else {
      # Keep the capacity target and vary only the number of ramp updates.
      quantile_method$blocks_to_target_rate <-
        setting_row$blocks_to_target_rate[[1]]
    }
    return(list(
      threshold_method = quantile_method,
      model_method = model_method_none,
      outcome_fn = outcome_attend,
      mean_fn = attendance_prob,
      family = stats::binomial(link = "logit")
    ))
  }

  if (scenario == "S4") {
    # Keep the six-update capacity ramp and vary only the revision penalty.
    revise_method$revision_penalty <- setting_row$revision_penalty[[1]]
    return(list(
      threshold_method = quantile_method,
      model_method = revise_method,
      outcome_fn = outcome_cvd,
      mean_fn = cvd_prob,
      family = stats::binomial(link = "logit")
    ))
  }
}

# Cumulative variance of the residual-history PCs, excluding final-risk
# variance.
# Leave entries as NA when fewer than n_components PCs are available.
adaptation_pca_summary <- function(fit_out, n_components = 5L) {
  cumulative <- rep(NA_real_, n_components)
  if (!is.null(fit_out$pca)) {
    # Use all residual eigenvalues as the denominator, not just retained PCs.
    variance <- fit_out$pca$sdev^2
    available <- min(length(variance), n_components)
    cumulative[seq_len(available)] <-
      cumsum(variance)[seq_len(available)] / sum(variance)
  }

  tibble::tibble(
    component = seq_len(n_components),
    cumulative_variance = cumulative
  )
}

# One row per block: deployed threshold, target rate, and observed assignment
# rate.
adaptation_threshold_summary <- function(sim_out) {
  sim_out$results |>
    dplyr::group_by(block) |>
    dplyr::summarise(
      threshold = dplyr::first(threshold_used),
      target_assignment_rate = dplyr::first(target_treat_rate),
      observed_assignment_rate = mean(treat),
      .groups = "drop"
    )
}

# Actual sampled assignments, summarized within risk bins in each run.
# Use common absolute risk bins across replications.
adaptation_assignment_summary <- function(sim_out, bin_width = 0.01) {
  # Save the final cutoff so the plot can center absolute-risk bins later.
  cutoff <- tail(sim_out$results$threshold_used, 1)

  sim_out$results |>
    dplyr::mutate(
      # Group actual assignments by risk at the time of assignment.
      risk_bin = floor(risk / bin_width) * bin_width + bin_width / 2
    ) |>
    dplyr::group_by(risk_bin) |>
    dplyr::summarise(
      mean_risk = mean(risk),
      final_threshold = cutoff,
      n_patients = dplyr::n(),
      n_assigned = sum(treat),
      assignment_proportion = mean(treat),
      .groups = "drop"
    )
}

# Simulate and summarize one replication of one setting.
# Returns performance and the summaries needed for the adaptation figure.
run_adaptation_once <- function(
    scenario,
    setting_row,
    cohort_df,
    replication,
    seed,
    block_size = 100L,
    n_blocks = 30L
) {
  # Reuse the same seed for a replication across settings. This gives each
  # setting the same random-number stream at the start of the simulation.
  set.seed(seed + replication)
  spec <- adaptation_scenario_spec(scenario, setting_row, cohort_df)

  sim_out <- simulate_design(
    cohort_df = cohort_df,
    risk_fn = pce_predict,
    block_size = block_size,
    n_blocks = n_blocks,
    threshold_adapt_method = spec$threshold_method,
    model_adapt_method = spec$model_method,
    outcome_fn = spec$outcome_fn
  )

  # Evaluate the proposed estimator at this replication's final cutoff.
  fit_out <- suppressWarnings(estimate_spline(sim_out, family = spec$family))
  estimate_result <- ate_at_threshold(fit_out, sim_out)

  # The outcome DGM uses original race-specific PCE risk, even when the
  # deployed starting model differs (Scenario 4).
  original_risk <- pce_predict(cohort_df)
  truth_mean <- function(cohort_df, idx, risk_baseline, treat) {
    spec$mean_fn(cohort_df, idx, original_risk[idx], treat)
  }
  truth_result <- ate_at_threshold_truth(sim_out, cohort_df, truth_mean)

  # Store errors and interval coverage relative to the same local truth.
  performance <- tibble::tibble(
    replication = replication,
    estimate = estimate_result$ate,
    standard_error = estimate_result$se,
    ci_lower = unname(estimate_result$ci["lo"]),
    ci_upper = unname(estimate_result$ci["hi"]),
    truth = truth_result$ate,
    error = estimate_result$ate - truth_result$ate,
    squared_error = (estimate_result$ate - truth_result$ate)^2,
    covered = estimate_result$ci["lo"] <= truth_result$ate &
      truth_result$ate <= estimate_result$ci["hi"],
    ci_width = unname(estimate_result$ci["hi"] - estimate_result$ci["lo"]),
    final_threshold = tail(sim_out$results$threshold_used, 1),
    overall_assignment_rate = mean(sim_out$results$treat)
  )

  # Keep the summaries used by the figure and numerical-results tables.
  list(
    performance = performance,
    thresholds = adaptation_threshold_summary(sim_out) |>
      dplyr::mutate(replication = replication, .before = 1),
    assignments = if (scenario == "S1") {
      adaptation_assignment_summary(sim_out) |>
        dplyr::mutate(replication = replication, .before = 1)
    } else tibble::tibble(),
    pca = adaptation_pca_summary(fit_out) |>
      dplyr::mutate(replication = replication, .before = 1),
    # Compare updated risks with the deployed starting model on the same cohort.
    risk_changes = if (scenario == "S4") {
      # Subtract each cohort member's starting risk from every later risk.
      differences <- sweep(sim_out$risks_mat, 1, sim_out$risks_mat[, 1], "-")
      tibble::tibble(
        replication = replication,
        block = seq_len(n_blocks),
        mean_absolute_change = colMeans(abs(differences))
      )
    } else tibble::tibble(),
    model_updates = sim_out$model_updates |>
      dplyr::mutate(replication = replication, .before = 1)
  )
}

# Stack one type of summary across replications and attach its setting labels.
bind_adaptation_component <- function(results, component, scenario, setting) {
  pieces <- lapply(results, `[[`, component)
  # Some summaries apply only to S1 or S4; omit their empty tables.
  pieces <- pieces[vapply(pieces, nrow, integer(1)) > 0L]
  if (length(pieces) == 0L) return(tibble::tibble())

  dplyr::bind_rows(pieces) |>
    dplyr::mutate(
      scenario = scenario,
      setting = setting,
      .before = 1
    )
}

# Run each setting in one scenario and collect its replication summaries.
run_adaptation_scenario <- function(
    scenario,
    grid,
    cohort_df,
    n_replications = 20L,
    seed = 20260930L,
    block_size = 100L,
    n_blocks = 30L
) {
  # Collect each type of result separately.
  output <- list(
    performance = list(),
    thresholds = list(),
    assignments = list(),
    pca = list(),
    model_updates = list(),
    risk_changes = list()
  )

  # Run each setting in the scenario's grid.
  for (setting_index in seq_len(nrow(grid))) {
    setting_row <- grid[setting_index, , drop = FALSE]
    setting_name <- setting_row$setting[[1]]
    message(scenario, ": ", setting_name)

    # Generate the requested replications for this setting.
    runs <- vector("list", n_replications)
    for (replication in seq_len(n_replications)) {
      runs[[replication]] <- run_adaptation_once(
        scenario = scenario,
        setting_row = setting_row,
        cohort_df = cohort_df,
        replication = replication,
        seed = seed,
        block_size = block_size,
        n_blocks = n_blocks
      )
    }

    # Combine replications and label them with the scenario and setting.
    for (component in c("performance", "thresholds", "assignments",
                        "pca", "model_updates", "risk_changes")) {
      output[[component]][[setting_name]] <- bind_adaptation_component(
        runs, component, scenario, setting_name
      )
    }
  }

  # Combine settings into one table for each type of result.
  lapply(output, dplyr::bind_rows)
}

# Run the requested scenarios, reusing complete settings when requested.
# Save scenario files and a combined file for plotting and Quarto tables.
run_adaptation <- function(
    n_replications = 500L,
    scenarios = c("S1", "S4"),
    grids = adaptation_default_grid(),
    cohort_df = NULL,
    results_dir = "results/adaptation_paper",
    seed = 20260930L,
    block_size = 100L,
    n_blocks = 30L,
    reuse_saved = TRUE
) {
  # Prepare the original cohort and output folder.
  if (is.null(cohort_df)) cohort_df <- prepare_adaptation_cohort(seed = seed)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  combined <- list()
  for (scenario in scenarios) {
    # Load this scenario's saved results, if reusing them.
    scenario_file <- file.path(
      results_dir, paste0("adaptation_", tolower(scenario), ".rds")
    )
    old <- if (reuse_saved && file.exists(scenario_file))
      readRDS(scenario_file) else NULL
    grid <- grids[[scenario]]

    # Identify settings with all requested replications already saved.
    complete <- vapply(grid$setting, function(setting) {
      if (is.null(old)) return(FALSE)
      rows <- old$performance[
        old$performance$setting == setting, , drop = FALSE
      ]
      all(seq_len(n_replications) %in% rows$replication)
    }, logical(1))
    if (any(complete)) message(
      "Reusing ", scenario, ": ",
      paste(grid$setting[complete], collapse = ", ")
    )

    # Run only incomplete settings; each receives the full replication count.
    missing <- grid[!complete, , drop = FALSE]
    fresh <- if (nrow(missing)) run_adaptation_scenario(
      scenario = scenario, grid = missing, cohort_df = cohort_df,
      n_replications = n_replications, seed = seed,
      block_size = block_size, n_blocks = n_blocks
    ) else old

    # Replace rerun settings and retain the other saved settings.
    scenario_result <- fresh
    if (!is.null(old) && nrow(missing)) {
      scenario_result <- lapply(names(fresh), function(component) {
        previous <- dplyr::bind_rows(old[[component]])
        if (nrow(previous))
          previous <- dplyr::filter(previous, !setting %in% missing$setting)
        dplyr::bind_rows(previous, fresh[[component]])
      })
      names(scenario_result) <- names(fresh)
    }

    # Save the scenario, then select the paper's settings and replications.
    saveRDS(scenario_result, scenario_file)
    message("Saved ", scenario_file)
    combined[[scenario]] <- lapply(scenario_result, function(component) {
      if (!nrow(component)) return(component)
      dplyr::filter(component, setting %in% grid$setting,
                    replication <= n_replications)
    })
  }

  # Combine scenarios and record the settings used for this analysis.
  combined_result <- list(
    config = list(
      n_replications = n_replications,
      scenarios = scenarios,
      grids = grids[scenarios],
      seed = seed,
      block_size = block_size,
      n_blocks = n_blocks
    ),
    performance = dplyr::bind_rows(lapply(combined, `[[`, "performance")),
    thresholds = dplyr::bind_rows(lapply(combined, `[[`, "thresholds")),
    assignments = dplyr::bind_rows(lapply(combined, `[[`, "assignments")),
    pca = dplyr::bind_rows(lapply(combined, `[[`, "pca")),
    model_updates = dplyr::bind_rows(lapply(combined, `[[`, "model_updates")),
    risk_changes = dplyr::bind_rows(lapply(combined, `[[`, "risk_changes"))
  )

  # Save the figure's input data and return it without printing large tables.
  saveRDS(combined_result, file.path(results_dir, "adaptation_results.rds"))
  invisible(combined_result)
}
