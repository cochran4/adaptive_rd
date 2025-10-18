# R/simulate.R
# Streamed simulation over iid bootstrap samples from a cohort,
# with optional threshold adaptation (quantile- or NNT-based).

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

# iid sample of row indices (with replacement)
sample_idx <- function(n_rows, size) sample.int(n_rows, size = size, replace = TRUE)

# -------------------------------------------------------------------
# Main API
# -------------------------------------------------------------------

# simulate_design()
# - Draws an initial warm-up block of individuals, followed by n_blocks
#   blocks of size block_size, sampled iid (with replacement) from cohort_df.
# - Assigns treatment if predicted risk exceeds a threshold, which can adapt
#   over time depending on threshold_adapt_method.
# - Optionally updates the outcome model itself according to model_adapt_method.
#
# Arguments:
#   cohort_df            : data frame, the sampling frame for iid draws
#   risk_fn              : function(df_block, coef_adjust?) -> vector of risks in [0,1]
#   initial_block_size   : number of individuals in the initial warm-up block
#   block_size           : number of individuals per subsequent block
#   n_blocks             : number of post–warm-up blocks
#   threshold_adapt_method : list describing threshold adaptation, with fields:
#                            $type = "none", "quantile", or "nnt"
#                            $initial_threshold = starting risk threshold in [0,1]
#                            $adapt_every = update frequency (in blocks)
#                            $desired_treat_rate = target treat rate if type=="quantile"
#                            $desired_nnt = target NNT if type=="nnt"
#                            $min_threshold, $max_threshold, $grid_points = search grid
#   model_adapt_method   : list describing model adaptation, with fields:
#                            $type = "none", "refit", or other strategies
#                            plus any method-specific options
#   outcome_fn           : function(cohort_df, idx, risk0, treat) -> outcome probability
#
# Returns:
#   list(
#     data = tibble(person_id, block, risk, threshold_used, treat, outcome),
#     threshold_log = tibble(block, threshold, adapt_reason)
#   )
simulate_design <- function(
    cohort_df,               # Base cohort data; rows are sampled with replacement
    risk_fn,                 # Function(df_block) -> vector of predicted risks in [0,1]
    initial_block_size,      # Sample size for the initial warm-up block
    block_size,              # Sample size for each subsequent block
    n_blocks,                # Number of post–warm-up blocks
    threshold_adapt_method,  # List with $type, $initial_threshold, etc.
    model_adapt_method,      # List with $type and optional adaptation logic
    outcome_fn               # Function(cohort_df, idx, risk0, treat) -> outcome 
) {
  N <- nrow(cohort_df)
  
  # Start with the user-specified initial threshold
  current_threshold <- threshold_adapt_method$initial_threshold
  
  # Storage for simulation output (1 warm-up block + n_blocks)
  results <- vector("list", length = 1L + n_blocks)
  
  # Index history across blocks (used for adaptation)
  idx_hist <- integer(0)
  
  # Matrix to store risk values for *all* cohort individuals over time
  risks_mat <- matrix(NA_real_, nrow = N, ncol = n_blocks + 1L)
  colnames(risks_mat) <- paste0("block_", 0:n_blocks)
  
  # Precompute risk scores for block 0
  risks_mat[, 1L] <- risk_fn(cohort_df, coef_adjust = NULL)
  
  # Iterate over each block (block 1 is warm-up)
  for (block in 1:(1 + n_blocks)) {
    
    # Set sample size for this block
    current_block_size <- if (block == 1L) initial_block_size else block_size
    
    # Sample a set of individual indices with replacement
    idx <- sample_idx(N, current_block_size)
    
    # Pull risk values for the sampled individuals at the current block
    risk <- risks_mat[idx, block]
    
    # Apply treatment rule: treat if above the current threshold
    treat <- as.integer(risk >= current_threshold)
    
    # Draw outcomes using the user-specified outcome function
    # Uses block-0 risk for consistency across treatment decisions
    outcome <- outcome_fn(cohort_df, idx, risks_mat[idx, 1L], treat)
    
    # Store results from this block
    results[[block]] <- tibble::tibble(
      person_id      = idx,
      block          = block,
      risk           = risk,
      threshold_used = current_threshold,
      treat          = treat,
      outcome        = outcome
    )
    
    # Accumulate index history across blocks (used for adaptive updates)
    idx_hist <- c(idx_hist, idx)

    # ----------- Adapt Threshold for Next Block -----------
    if (block < 1 + n_blocks) {
      
      # Update model using adaptation method
      model_adapt_method$coef_adjusted <- adapt_model(
        df      = cohort_df,
        method  = model_adapt_method,
        sim_out = list(
          results   = dplyr::bind_rows(results[1:block]),
          risks_mat = risks_mat[ , 1:block, drop = FALSE]
        )
      )
      
      # Update threshold using adaptation method
      risks_cumulative <- risks_mat[idx_hist, block]
      current_threshold <- adapt_threshold(
        method  = threshold_adapt_method,
        risks   = risks_cumulative,
        sim_out = list(
          results   = dplyr::bind_rows(results[1:block]),
          risks_mat = risks_mat[ , 1:block, drop = FALSE]
        )
      )
      
      # Pre-compute risk for next block
      risks_mat[, block + 1L] <- risk_fn(cohort_df, coef_adjust = model_adapt_method$coef_adjusted)
      
    }
  }
  
  # Return full results and risk matrix
  list(
    results   = dplyr::bind_rows(results),
    risks_mat = risks_mat
  )
}
  