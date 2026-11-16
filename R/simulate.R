# R/simulate.R
# Streamed simulation over iid bootstrap samples from a cohort,
# with optional threshold adaptation (quantile- or benefit-based).

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

# iid sample of row indices (with replacement)
sample_idx <- function(n_rows, size)
  sample.int(n_rows, size = size, replace = TRUE)

# -------------------------------------------------------------------
# Main API
# -------------------------------------------------------------------

# simulate_design()
# - Draws n_blocks of size block_size, sampled iid (with replacement) from
#   cohort_df.
# - Assigns treatment if predicted risk meets or exceeds a threshold, which can
# adapt
#   over time depending on threshold_adapt_method.
# - Optionally updates the prediction model according to model_adapt_method.
#
# Arguments:
#   cohort_df            : data frame, the sampling frame for iid draws
# risk_fn              : function(df_block, coef_adjust?) -> vector of risks in
# [0,1]
#   block_size           : number of individuals per block
#   n_blocks             : total number of blocks
#   threshold_adapt_method : list describing threshold adaptation, with fields:
#                            $type = "none", "quantile", or "benefit"
# $initial_threshold = starting threshold for the first block
# $target_treat_rate = final target rate if type=="quantile"
# $blocks_to_target_rate = ramp length if type=="quantile"
# $target_reduction = target outcome reduction if type=="benefit"
#   model_adapt_method   : list describing model adaptation, with fields:
#                            $type = "none", "recalibrate", or "revise"
# $coef_initial = optional starting prediction coefficients
#                            plus any method-specific options
# outcome_fn           : function(cohort_df, idx, risk0, treat) -> simulated
# outcome
#
# Returns:
#   list(
# results = tibble(person_id, block, risk, threshold_used, treat, outcome, ...),
#     risks_mat = cohort risks under each deployed model,
#     model_updates = fitted model parameters after each update
#   )
simulate_design <- function(
    # Base cohort data; rows are sampled with replacement
    cohort_df,
    # Function(data, coef_adjust) -> predicted risks in [0,1]
    risk_fn,
    block_size,              # Sample size for each block
    n_blocks,                # Total number of blocks
    threshold_adapt_method,   # List with $type and threshold settings
    model_adapt_method,      # List with $type and optional adaptation logic
    outcome_fn               # Function(cohort_df, idx, risk0, treat) -> outcome
) {
  N <- nrow(cohort_df)

  # Start with the threshold specified by the scenario.
  current_threshold <- threshold_adapt_method$initial_threshold

  # Record the planned rate for quantile rules; other rules have no rate target.
  current_target_treat_rate <- if (
    identical(threshold_adapt_method$type, "quantile")
  ) {
    0
  } else {
    NA_real_
  }

  # Store one result per block and at most one model update between blocks.
  results <- vector("list", length = n_blocks)
  model_updates <- vector("list", length = max(n_blocks - 1L, 0L))

  # Retain all sampled cohort indices, including repeated draws.
  idx_hist <- integer(0)

  # Rows are source-cohort individuals; columns are deployed block models.
  risks_mat <- matrix(NA_real_, nrow = N, ncol = n_blocks)
  colnames(risks_mat) <- paste0("block_", seq_len(n_blocks))

  # Preserve risk under the original model for the outcome-generating mechanism.
  risk_baseline <- risk_fn(cohort_df, coef_adjust = NULL)

  # Precompute assignment risk for block 1. Model-revision scenarios can supply
  # their revised starting model before any outcomes have accumulated.
  initial_coef_adjust <- model_adapt_method$coef_initial
  risks_mat[, 1L] <- if (is.null(initial_coef_adjust)) {
    risk_baseline
  } else {
    risk_fn(cohort_df, coef_adjust = initial_coef_adjust)
  }

  # Process n_blocks equal-sized blocks; there is no separate warm-up block.
  for (block in seq_len(n_blocks)) {

    # Sample a set of individual indices with replacement
    idx <- sample_idx(N, block_size)

    # Pull risk values for the sampled individuals at the current block
    risk <- risks_mat[idx, block]

    # Assign the intervention when risk meets or exceeds the current threshold.
    treat <- as.integer(risk >= current_threshold)

    # Draw outcomes using risk under the original model so that changing the
    # deployed prediction model does not change the outcome-generating
    # mechanism.
    outcome <- outcome_fn(cohort_df, idx, risk_baseline[idx], treat)

    # Store results from this block
    results[[block]] <- tibble::tibble(
      person_id        = idx,
      block            = block,
      risk             = risk,
      threshold_used   = current_threshold,
      target_treat_rate = current_target_treat_rate,
      treat            = treat,
      outcome          = outcome
    )

    # Accumulate index history across blocks (used for adaptive updates)
    idx_hist <- c(idx_hist, idx)

    # Update the model and threshold for the next block; skip after the last
    # block.
    if (block < n_blocks) {

      # Fit the next model using only outcomes observed through this block.
      model_adapt_method$coef_adjusted <- adapt_model(
        df      = cohort_df,
        method  = model_adapt_method,
        sim_out       = list(
          results   = dplyr::bind_rows(results[1:block]),
          risks_mat = risks_mat[, 1:block, drop = FALSE]
        )
      )

      # Save any fitted parameters attached by the model updater.
      recalibration <- attr(model_adapt_method$coef_adjusted, "recalibration")
      if (!is.null(recalibration)) {
        model_updates[[block]] <- tibble::as_tibble(recalibration)
        model_updates[[block]]$after_block <- block
      }

      revision <- attr(model_adapt_method$coef_adjusted, "revision")
      if (!is.null(revision)) {
        model_updates[[block]] <- tibble::as_tibble(revision)
        model_updates[[block]]$after_block <- block
      }

      # Use observed patients' risks under the just-completed block's model.
      # The newly fitted model supplies predictions only for the next block.
      risks_cumulative <- risks_mat[idx_hist, block]
      current_threshold <- adapt_threshold(
        method        = threshold_adapt_method,
        risks         = risks_cumulative,
        update_number = block,
        sim_out       = list(
          results   = dplyr::bind_rows(results[1:block]),
          risks_mat = risks_mat[, 1:block, drop = FALSE]
        )
      )

      # Record the ramped assignment target for the next block.
      if (identical(threshold_adapt_method$type, "quantile")) {
        current_target_treat_rate <- quantile_target_rate(
          threshold_adapt_method, block
        )
      }

      # Evaluate the next model for the full source cohort.
      risks_mat[, block + 1L] <- risk_fn(
        cohort_df, coef_adjust = model_adapt_method$coef_adjusted
      )

    }
  }

  # Combine block results and update logs; retain the full risk history.
  list(
    results       = dplyr::bind_rows(results),
    risks_mat     = risks_mat,
    model_updates = dplyr::bind_rows(model_updates)
  )
}
