
# =========================
# Styling helpers (reusable)
# =========================

# Okabe–Ito palette (color-blind friendly)
cb_palette <- c(
  "Control" = "#666666", # Dark gray, strong but neutral
  "Treated" = "#0072B2", # Blue (instead of green, closer luminance)
  "Accent"  = "#D55E00"  # Vermillion/orange-red
)

# Some controls
POINT_ALPHA = .1
ZOOM_FRAC   = .2
N_POINTS    = 200
BASE_SIZE   = 14
BASE_FAMILY = "sans"

# Modern, clean theme (no global side effects)
theme_modern <- function() {
  ggplot2::theme_minimal(base_size = BASE_SIZE, base_family = BASE_FAMILY) + 
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold"),
      plot.subtitle   = ggplot2::element_text(color = "gray30"),
      axis.title      = ggplot2::element_text(face = "bold"),
      axis.text       = ggplot2::element_text(color = "gray20"),
      panel.grid.minor= ggplot2::element_blank(),
      panel.grid.major= ggplot2::element_line(color = "gray90", linewidth = 0.4),
      axis.ticks      = ggplot2::element_blank(),
      legend.position  = "none",
      legend.title     = ggplot2::element_text(size = BASE_SIZE, face = "bold"),
      legend.text      = ggplot2::element_text(size = BASE_SIZE),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

# Convenience scales for our two groups
scale_colors_cb2 <- function() ggplot2::scale_color_manual(values = cb_palette[c("Control","Treated")], name = "")
scale_fills_cb2  <- function() ggplot2::scale_fill_manual(values  = cb_palette[c("Control","Treated")], name = "")


# ============================
# Plot threshold adaptation trajectory
# ============================
# This function takes simulation output and returns a ggplot showing how
# the referral threshold changed over time as patients were seen.
# Optionally, the plot can be saved as an RDS object for later use.

plot_threshold_trajectory <- function(sim_out, save_path = NULL) {
  
  # ----------------------------------
  # Prepare threshold data
  # ----------------------------------
  thr <- sim_out$results %>%
    dplyr::arrange(block) %>%  # sort by block if relevant
    dplyr::mutate(
      patient_num = dplyr::row_number(),                   # running patient index
      jump = threshold_used != dplyr::lag(threshold_used), # flag threshold jumps
      jump = tidyr::replace_na(jump, TRUE)                 # first point is always a jump
    )
  
  # ----------------------------------
  # Create plot object
  # ----------------------------------
  p <- ggplot2::ggplot(thr, ggplot2::aes(x = patient_num, y = threshold_used)) +
    ggplot2::geom_step(linewidth = 1.1, color = cb_palette["Accent"]) +           # step plot of threshold
    ggplot2::geom_point(                                                           # mark jumps with open circles
      data  = dplyr::filter(thr, jump),
      size  = 2.6, shape = 21, stroke = 0.9, fill = "white", color = cb_palette["Accent"]
    ) +
    ggplot2::labs(
      x = "Patient number (order seen)",
      y = "Referral threshold"
    ) +
    theme_modern()                                                                 # custom modern theme
  
  # ----------------------------------
  # Save plot if a path is provided
  # ----------------------------------
  if (!is.null(save_path)) {
    saveRDS(p, save_path)
  }
  
  # ----------------------------------
  # Return the plot object
  # ----------------------------------
  return(p)
}

#------------------------------------------------------------------------------#
# plot_spline_curves
#------------------------------------------------------------------------------#
# Plots the conditional mean of a binary outcome (e.g., attendance) as a function 
# of predicted risk. Shows fitted spline curves, ribbons (CIs), and ground truth.
# Optionally annotates the average treatment effect (ATE) at the threshold.
#
# Inputs:
#   - fit_out: model output with fitted splines and input data
#   - sim_out: simulation output used to get last threshold
#   - truth_risk_fn: ground truth risk function
#   - truth_mean_fn: ground truth mean function for E[Y|R]
#   - cohort_df: original simulated data used to get true response curves
#   - ci: width of confidence interval for fitted curves
#   - y_lab: label for y-axis
#   - labels: vector with control and treated group labels
#   - save_path: optional file path to save the plot:
#         * If ends with ".rds" → saves the ggplot object
#         * If ends with ".png" → saves a PNG image using ggsave()
#
# Output:
#   - ggplot object (also saved if save_path provided)

plot_spline_curves <- function(
    fit_out, sim_out, truth_risk_fn, truth_mean_fn, cohort_df,
    ci = 0.95,
    y_lab = "Probability of attendance",
    labels = c("Not referred", "Referred"),
    save_path = NULL   # <--- new argument
) {
  
  # ----------------------------------------
  # Setup: label definitions and data prep
  # ----------------------------------------
  lbl_ctrl <- labels[1]
  lbl_trt  <- labels[2]
  fit      <- fit_out$fit
  dat      <- fit_out$data
  n_patients <- nrow(dat)
  
  # Get the last threshold used (for reference line)
  cutoff <- tail(sim_out$results$threshold_used, 1)
  
  # Compute plotting window
  x_min <- min(dat$risk_last, na.rm = TRUE)
  x_max <- max(dat$risk_last, na.rm = TRUE)
  span  <- x_max - x_min
  #zoom_min <- max(x_min, cutoff - ZOOM_FRAC * span)
  #zoom_max <- min(x_max, cutoff + ZOOM_FRAC * span)
  zoom_min <- .00
  zoom_max <- .40
  
  # Generate sequence of x-values for smooth plotting
  x_seq <- seq(zoom_min, zoom_max, length.out = N_POINTS)
  
  # ----------------------------------------
  # Ground truth curves
  # ----------------------------------------
  cc_truth <- conditional_curves_truth(
    cohort_df, sim_out$risks_mat, truth_mean_fn, x_seq
  )
  cc_truth$grp <- ifelse(cc_truth$group == "Control", lbl_ctrl, lbl_trt)
  
  # ----------------------------------------
  # Fitted curves with CIs
  # ----------------------------------------
  cc <- conditional_curves(fit_out, x_seq, level = ci)
  cc$grp <- ifelse(cc$grp == "Control", lbl_ctrl, lbl_trt)
  
  # Split curves by side of threshold for alpha transparency
  s_ctrl <- subset(cc, grp == lbl_ctrl)
  s_trt  <- subset(cc, grp == lbl_trt)
  s_ctrl_left  <- subset(s_ctrl, risk_last <= cutoff); s_ctrl_left$alpha <- 1
  s_ctrl_right <- subset(s_ctrl, risk_last  > cutoff); s_ctrl_right$alpha <- 0.7
  s_trt_left   <- subset(s_trt,  risk_last  < cutoff); s_trt_left$alpha  <- 0.7
  s_trt_right  <- subset(s_trt,  risk_last >= cutoff); s_trt_right$alpha <- 1
  ss <- rbind(s_ctrl_left, s_ctrl_right, s_trt_left, s_trt_right)
  
  # ----------------------------------------
  # Observed data within plotting window
  # ----------------------------------------
  obs <- subset(dat, risk_last >= zoom_min & risk_last <= zoom_max)
  if ("t" %in% names(obs)) {
    obs$grp <- ifelse(obs$t == 0, lbl_ctrl, lbl_trt)
  } else {
    obs$grp <- ifelse(as.integer(as.character(obs$treat)) == 0, lbl_ctrl, lbl_trt)
  }
  
  # Set consistent factor levels
  lvl <- c(lbl_ctrl, lbl_trt)
  obs$grp      <- factor(obs$grp, levels = lvl)
  ss$grp       <- factor(ss$grp,  levels = lvl)
  cc_truth$grp <- factor(cc_truth$grp, levels = lvl)
  
  # ----------------------------------------
  # Annotate ATE at threshold
  # ----------------------------------------
  ate_res <- ate_at_threshold(fit_out, sim_out)
  y0 <- as.numeric(ate_res$control)
  y1 <- as.numeric(ate_res$treated)
  y_mid <- (y0 + y1) / 2
  
  lab <- bquote(beta[i] == .(sprintf("%.2f (%.2f, %.2f)",
                                     ate_res$ate, ate_res$ci["lo"], ate_res$ci["hi"])))
  
  dx   <- 0.0125 * span / 3
  x_br <- if (cutoff + 3 * dx <= zoom_max) cutoff + 2 * dx else cutoff - 2 * dx
  
  # Dummy data for ATE star in legend
  legend_df <- data.frame(x = zoom_min, y = min(ss$l, na.rm = TRUE), key = "ATE")
  
  # Get accent color
  accent_col <- unname(cb_palette["Accent"])
  
  # ----------------------------------------
  # Begin ggplot
  # ----------------------------------------
  p <- ggplot2::ggplot()
  
  # Observed points (jittered)
  p <- p + ggplot2::geom_point(
    data = obs,
    mapping = ggplot2::aes(risk_last, y, color = grp),
    position = ggplot2::position_jitter(height = 0.05, width = 0),
    alpha = POINT_ALPHA, size = 1.2
  )
  
  # Confidence bands (ribbons)
  p <- p + ggplot2::geom_ribbon(
    data = ss,
    mapping = ggplot2::aes(risk_last, ymin = l, ymax = u, fill = grp),
    alpha = 0.18, color = NA, show.legend = FALSE
  )
  
  # Fitted curves
  p <- p + ggplot2::geom_line(
    data = ss,
    mapping = ggplot2::aes(risk_last, mean, color = grp, alpha = alpha, linetype = "solid"),
    linewidth = 1.1
  )
  
  # Dummy point for ATE legend entry (invisible)
  p <- p + ggplot2::geom_point(
    data = legend_df,
    mapping = ggplot2::aes(x, y, shape = key),
    color = accent_col, alpha = 0, size = 4, inherit.aes = FALSE
  )
  
  # Reference line: referral threshold
  p <- p + ggplot2::geom_vline(
    xintercept = cutoff, linetype = 2, linewidth = 0.9, color = accent_col
  )
  
  # Ground truth curves (dashed)
  p <- p + ggplot2::geom_line(
    data = cc_truth, linetype = "dashed",
    mapping = ggplot2::aes(risk_last, mean, color = grp),
    linewidth = 1.0, alpha = 0.9, show.legend = FALSE
  )
  
  # ----------------------------------------
  # Scales & guides
  # ----------------------------------------
  p <- p + ggplot2::scale_color_manual(
    name = "", values = setNames(unname(cb_palette[c("Control", "Treated")]), lvl),
    breaks = lvl, labels = lvl
  )
  
  p <- p + ggplot2::scale_fill_manual(
    values = setNames(unname(cb_palette[c("Control", "Treated")]), lvl),
    breaks = lvl, labels = lvl, guide = "none"
  )
  
  p <- p + ggplot2::scale_shape_manual(
    name = "", values = c("ATE" = 42), labels = c("ATE" = lab)
  )
  
  p <- p + ggplot2::scale_linetype_identity(guide = "none")
  p <- p + ggplot2::scale_alpha_identity()
  
  # ----------------------------------------
  # Axis labels and limits
  # ----------------------------------------
  p <- p + ggplot2::coord_cartesian(xlim = c(zoom_min, zoom_max)) +
    ggplot2::labs(
      x = bquote(bold("Risk prediction, " * R[.(n_patients)])),
      y = y_lab
    ) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(face = "bold"))
  
  # ----------------------------------------
  # Theme and legend
  # ----------------------------------------
  p <- p + ggplot2::guides(
    color = ggplot2::guide_legend(order = 1, override.aes = list(linewidth = 1.1)),
    shape = "none"
  )
  
  p <- p + theme_modern() +
    ggplot2::theme(legend.position = "bottom")
  
  # ----------------------------------------
  # Save if requested
  # ----------------------------------------
  if (!is.null(save_path)) {
    if (grepl("\\.rds$", save_path, ignore.case = TRUE)) {
      saveRDS(p, save_path)
    } else if (grepl("\\.png$", save_path, ignore.case = TRUE)) {
      ggplot2::ggsave(save_path, plot = p, width = 8, height = 5, dpi = 300)
    } else {
      warning("save_path must end with .rds or .png to trigger saving.")
    }
  }
  
  return(p)
}


#------------------------------------------------------------------------------#
# plot_cohen_d_to_nnt
#------------------------------------------------------------------------------#
# Plots the relationship between Cohen's d and the number needed to treat (NNT)
# using the approximation: NNT = 1 / (2 * Φ(d / √2) - 1)
# Adds benchmark vertical lines for small, medium, and large effects.
#
# Inputs:
#   - target_nnt: numeric, the desired NNT to highlight on the plot
#   - save_path: optional file path to save the plot:
#         * If ends with ".rds" → saves the ggplot object
#         * If ends with ".png", ".pdf", or ".svg" → saves as image using ggsave()
#
# Output:
#   - ggplot object (also saved if save_path provided)

plot_cohen_d_to_nnt <- function(target_nnt, save_path = NULL) {
  
  # ----------------------------------------
  # Define transformation: d → NNT
  # ----------------------------------------
  d_to_nnt <- function(d) 1 / (2 * pnorm(d / sqrt(2)) - 1)
  
  # ----------------------------------------
  # Create data for curve
  # ----------------------------------------
  d_vals   <- seq(0.1, 0.9, by = 0.01)
  nnt_vals <- d_to_nnt(d_vals)
  df       <- data.frame(d = d_vals, nnt = nnt_vals)
  
  # ----------------------------------------
  # Cohen's d benchmarks
  # ----------------------------------------
  benchmarks <- data.frame(
    d     = c(0.2, 0.5, 0.8),
    label = c("Small", "Medium", "Large")
  )
  
  # ----------------------------------------
  # Build plot with swapped axes
  # ----------------------------------------
  p <- ggplot(df, aes(x = nnt, y = d)) +
    
    # Curve
    geom_line(color = cb_palette["Treated"], linewidth = 1.2) +
    
    # Target NNT vertical reference line
    geom_vline(
      xintercept = target_nnt,
      color = cb_palette["Accent"],
      linetype = "dashed", linewidth = 1.2
    ) +
    
    # Target NNT label
    annotate("text",
             x = target_nnt + .02, y = max(d_vals),
             label = paste0("Target NNT = ", target_nnt),
             vjust = 1, hjust = 0,
             color = cb_palette["Accent"],
             fontface = "bold"
    ) +
    
    # Horizontal benchmark lines for d
    geom_hline(
      data = benchmarks,
      aes(yintercept = d),
      linetype = "dotted", color = "gray50"
    ) +
    
    # Benchmark labels
    geom_text(
      data = benchmarks,
      aes(x = max(nnt_vals), y = d, label = label),
      color = "gray30",
      size = BASE_SIZE / ggplot2::.pt,
      hjust = 1, vjust = -0.5
    ) +
    
    # Axis labels and theme
    labs(
      x = "Number needed to treat (NNT)",
      y = "Cohen's d"
    ) +
    theme_modern()
  
  # ----------------------------------------
  # Save if requested
  # ----------------------------------------
  if (!is.null(save_path)) {
    if (grepl("\\.rds$", save_path, ignore.case = TRUE)) {
      saveRDS(p, save_path)
    } else if (grepl("\\.(png|pdf|svg)$", save_path, ignore.case = TRUE)) {
      ggsave(save_path, plot = p, width = 7, height = 4.5, dpi = 300)
    } else {
      warning("save_path must end with .rds, .png, .pdf, or .svg to trigger saving.")
    }
  }
  
  return(p)
}


#------------------------------------------------------------------------------#
# plot_risk_to_cohen_d
#------------------------------------------------------------------------------#
# Plots how baseline risk maps to Cohen's d, assuming a treatment effect 
# that scales linearly with baseline risk. Highlights a target d value.
#
# Inputs:
#   - target_d: Numeric. Desired Cohen’s d to highlight on the plot.
#   - treatment_effect_slope: Slope relating baseline risk to treatment effect.
#   - outcome_sd: Standard deviation of the outcome variable.
#   - save_path: Optional path to save the plot:
#         * If ends with ".rds" → saves ggplot object
#         * If ends with ".png", ".pdf", or ".svg" → saves plot image
#
# Output:
#   - ggplot object (also saved if save_path is provided)
plot_risk_to_cohen_d <- function(target_d, treatment_effect_slope, outcome_sd, save_path = NULL) {
  
  # ----------------------------------------
  # Risk grid and effect size calculation
  # ----------------------------------------
  risk_vals     <- seq(0, 1, by = 0.01)                          # Range of baseline risks
  effect_vals   <- treatment_effect_slope * risk_vals           # Treatment effect at each risk
  cohen_d_vals  <- abs(effect_vals) / outcome_sd                # Convert to Cohen's d
  df            <- data.frame(risk = risk_vals, cohen_d = cohen_d_vals)
  
  # Compute and print the risk level needed to achieve target_d
  required_risk <- target_d * outcome_sd / treatment_effect_slope
  message(sprintf("To achieve Cohen's d = %.2f, risk must be approximately %.3f", target_d, required_risk))
  
  # ----------------------------------------
  # Build plot with axes swapped
  # ----------------------------------------
  p <- ggplot(df, aes(x = cohen_d, y = risk)) +
    
    # Line for predicted risk
    geom_line(color = cb_palette["Treated"], linewidth = 1.2) +
    
    # Target d vertical threshold
    geom_vline(
      xintercept = target_d,
      color = cb_palette["Accent"],
      linetype = "dashed", linewidth = 1.2
    ) +
    
    # Annotate target d
    annotate("text",
             x = target_d, y = max(risk_vals),
             label = paste0("Target d = ", round(target_d, 2)),
             hjust = -0.1, vjust = 1.2,
             color = cb_palette["Accent"], fontface = "bold"
    ) +
    
    # Labels and theme
    labs(
      x = "Cohen's d",
      y = "Predicted risk"
    ) +
    theme_modern()
  
  # ----------------------------------------
  # Save if requested
  # ----------------------------------------
  if (!is.null(save_path)) {
    if (grepl("\\.rds$", save_path, ignore.case = TRUE)) {
      saveRDS(p, save_path)
    } else if (grepl("\\.(png|pdf|svg)$", save_path, ignore.case = TRUE)) {
      ggsave(save_path, plot = p, width = 7, height = 4.5, dpi = 300)
    } else {
      warning("save_path must end with .rds, .png, .pdf, or .svg to trigger saving.")
    }
  }
  
  return(p)
}

#------------------------------------------------------------------------------#
# plot_risk_trajectory
#------------------------------------------------------------------------------#
# Plots how patient risk scores evolve over time (across blocks).
# Shows median and interquartile range across patients, with a spaghetti
# plot for a subset of individual risk trajectories.
#
# Inputs:
#   - sim_out: A simulation output object containing either `risks_mat` or `risk_mat`
#              (each row = patient, each col = time block)
#   - n_max:   Max number of individual trajectories to show as spaghetti lines
#   - save_path: Optional path to save plot (.rds, .png, .pdf, .svg supported)
#
# Output:
#   - ggplot object (also saved if save_path is provided)
plot_risk_trajectory <- function(sim_out, n_max = 100, save_path = NULL) {
  
  # ----------------------------------------
  # Extract matrix of risks
  # ----------------------------------------
  R <- if (!is.null(sim_out$risks_mat)) sim_out$risks_mat else sim_out$risk_mat
  stopifnot(is.matrix(R))
  
  # ----------------------------------------
  # Convert to long format (patient, block, risk)
  # ----------------------------------------
  df_all <- as.data.frame(R)
  df_all$person <- seq_len(nrow(df_all))
  df_all <- tidyr::pivot_longer(df_all, -person, names_to = "block", values_to = "risk")
  df_all$block <- as.integer(gsub("\\D", "", df_all$block))  # strip non-numeric
  
  # ----------------------------------------
  # Summary statistics: median and IQR per block
  # ----------------------------------------
  df_summary <- df_all |>
    dplyr::group_by(block) |>
    dplyr::summarise(
      median = median(risk, na.rm = TRUE),
      p25    = quantile(risk, 0.25, na.rm = TRUE),
      p75    = quantile(risk, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ----------------------------------------
  # Select subset of patients for spaghetti lines
  # ----------------------------------------
  df_spaghetti <- dplyr::filter(
    df_all,
    person %in% sim_out$results$person_id[1:min(n_max, nrow(sim_out$results))]
  )
  
  # ----------------------------------------
  # Plot construction
  # ----------------------------------------
  p <- ggplot2::ggplot() +
    
    # Spaghetti lines
    ggplot2::geom_line(
      data = df_spaghetti,
      aes(x = block, y = risk, group = person),
      color = "gray70", alpha = 0.5, linewidth = 0.4
    ) +
    
    # # Interquartile ribbon
    # ggplot2::geom_ribbon(
    #   data = df_summary,
    #   aes(x = block, ymin = p25, ymax = p75),
    #   fill = cb_palette["Accent"], alpha = 0.15
    # ) +
    # 
    # # Median line
    # ggplot2::geom_line(
    #   data = df_summary,
    #   aes(x = block, y = median),
    #   color = cb_palette["Accent"], linewidth = 1.2
    # ) +
    
    # Axis labels and formatting
    ggplot2::labs(x = "Update", y = "Risk") +
    theme_modern()
  
  # ----------------------------------------
  # Save if requested
  # ----------------------------------------
  if (!is.null(save_path)) {
    if (grepl("\\.rds$", save_path, ignore.case = TRUE)) {
      saveRDS(p, save_path)
    } else if (grepl("\\.(png|pdf|svg)$", save_path, ignore.case = TRUE)) {
      ggsave(save_path, plot = p, width = 7, height = 4.5, dpi = 300)
    } else {
      warning("save_path must end with .rds, .png, .pdf, or .svg to trigger saving.")
    }
  }
  
  return(p)
}


#------------------------------------------------------------------------------#
# plot_comparison
#------------------------------------------------------------------------------#
# Plots Monte Carlo estimation errors across scenarios and methods.
# Three panels are created based on groupings of scenarios:
# - Row 1: Scenario S1
# - Row 2: Scenarios S2 and S3
# - Row 3: Scenarios S4 and S5
#
# Inputs:
#   - mc_df:           Data frame with columns: method, scenario, est, truth
#   - method_labels:   Optional named vector for relabeling methods
#   - scenario_labels: Optional named vector for relabeling scenarios
#   - ylim_list:       Optional list of y-axis limits for each panel row
#   - save_path:       Optional file path to save plot (e.g., "plot.pdf" or "plot.png")
#
# Output:
#   - Combined patchwork plot with estimation error boxplots, faceted by scenario
plot_comparison <- function(mc_df, 
                            method_labels = NULL, 
                            scenario_labels = NULL,
                            width = 16,
                            height = 10,
                            dpi = 600,
                            ylim_list = list(
                              c(-0.2, 0.21),  # For row 1 (df1)
                              c(-3, 3),      # For row 2 (df23)
                              c(-0.2, 0.3)   # For row 3 (df45)
                            ),
                            save_path = NULL, 
                            save_path2 = NULL) {
  
  # ----------------------------------------
  # Define custom color palette by method
  # ----------------------------------------
  # Accent color for spline, grays for others
  color_palette <- c(
    spline = cb_palette[["Accent"]],
    naive  = "gray15",
    outreg = "gray40",
    ipw    = "gray65",
    aipw  = "gray90"
  )
  
  # Optional patterns to differentiate comparators
  pattern_palette <- c(
    "Adaptive RD" = "none",
    "Outcome Regression" = "none",
    "Naive" = "none",
    "Inverse Probability Weighting" = "none",
    "Augmented Inverse Probability Weighting" = "none"
  )
  
  # ----------------------------------------
  # Compute estimation error
  # ----------------------------------------
  mc_df <- mc_df %>%
    mutate(est_error = est - truth)
  
  # ----------------------------------------
  # Subset into three scenario groups for layout
  # ----------------------------------------
  df1  <- filter(mc_df, scenario == "S1")
  df23 <- filter(mc_df, scenario %in% c("S2", "S3"))
  df45 <- filter(mc_df, scenario %in% c("S4", "S5"))
  
  # ----------------------------------------
  # Apply optional relabeling for methods and scenarios
  # ----------------------------------------
  if (!is.null(method_labels)) {
    mc_df$method <- factor(mc_df$method, levels = names(method_labels), labels = method_labels)
    df1$method   <- factor(df1$method,   levels = names(method_labels), labels = method_labels)
    df23$method  <- factor(df23$method,  levels = names(method_labels), labels = method_labels)
    df45$method  <- factor(df45$method,  levels = names(method_labels), labels = method_labels)
    
    # Update color palette to match new labels
    color_palette <- setNames(
      color_palette[names(method_labels)],
      method_labels
    )
  }
  
  if (!is.null(scenario_labels)) {
    mc_df$scenario <- factor(mc_df$scenario, levels = names(scenario_labels), labels = scenario_labels)
    df1$scenario   <- factor(df1$scenario,   levels = names(scenario_labels), labels = scenario_labels)
    df23$scenario  <- factor(df23$scenario,  levels = names(scenario_labels), labels = scenario_labels)
    df45$scenario  <- factor(df45$scenario,  levels = names(scenario_labels), labels = scenario_labels)
  }
  
  # ----------------------------------------
  # Helper function to make boxplots for one set of scenarios
  # ----------------------------------------
  base_plot <- function(df, ylim = NULL) {
    p <- ggplot(df, aes(x = method, y = est_error, fill = method, color = method, pattern = method)) +
      ggpattern::geom_boxplot_pattern(
        outlier.shape = NA,
        alpha = 0.9,
        pattern_density = 0.5,
        pattern_spacing = 0.04,
        pattern_key_scale_factor = 0.6,
        linewidth = 0.3
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
      scale_fill_manual(values = color_palette) +
      scale_color_manual(values = color_palette) +
      scale_pattern_manual(values = pattern_palette) +
      labs(y = "Estimation Error", x = NULL) +
      theme_modern() +
      theme(
        strip.background = element_rect(fill = "gray90", color = "gray50"),
        strip.text = element_text(face = "bold", size = 18),
        panel.border = element_rect(color = "gray50", fill = NA),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()
      )
    
    if (!is.null(ylim)) {
      p <- p + coord_cartesian(ylim = ylim)
    }
    
    return(p)
  }
  
  # ----------------------------------------
  # Dummy plot for extracting the legend
  # ----------------------------------------
  legend_plot <- ggplot(mc_df, aes(x = method, y = est_error, fill = method, color = method, pattern = method)) +
    ggpattern::geom_boxplot_pattern(
      outlier.shape = NA,
      alpha = 0.9,
      pattern_density = 0.5,
      pattern_spacing = 0.04,
      pattern_key_scale_factor = 0.6,
      linewidth = 0.3
    ) +
    scale_fill_manual(values = color_palette, name = "Method") +
    scale_color_manual(values = color_palette, name = "Method") +
    scale_pattern_manual(values = pattern_palette, name = "Method") +
    theme_modern() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text  = element_text(size = 14)
    )
  
  legend <- cowplot::get_legend(legend_plot)
  
  # ----------------------------------------
  # Generate panel rows
  # ----------------------------------------
  p1  <- base_plot(df1,  ylim = ylim_list[[1]]) +
    facet_wrap(~scenario, ncol = 1) +
    theme(legend.position = "none")
  
  p23 <- base_plot(df23, ylim = ylim_list[[2]]) +
    facet_wrap(~scenario, ncol = 2) +
    theme(legend.position = "none")
  
  p45 <- base_plot(df45, ylim = ylim_list[[3]]) +
    facet_wrap(~scenario, ncol = 2) +
    theme(legend.position = "none")
  
  # ----------------------------------------
  # Assemble final layout
  # ----------------------------------------
  #final_plot <- (p1 | patchwork::wrap_elements(legend)) / p23 / p45 +
  #  patchwork::plot_layout(heights = c(1, 1, 1))
  # ----------------------------------------
  # Assemble final layout
  # ----------------------------------------
  
  # Row 1: scenarios 1–3
  row1 <- p1 + p23 +
    plot_layout(ncol = 2, widths = c(1, 2))
  
  # Row 2: scenarios 4–5 and legend
  row2 <- patchwork::wrap_elements(legend) + p45 + 
    plot_layout(ncol = 2, widths = c(1, 2))
  
  # Combine rows vertically
  final_plot <- row1 / row2 +
    plot_layout(heights = c(1, 1))
  
  # ----------------------------------------
  # Optional: Save to file
  # ----------------------------------------
  if (!is.null(save_path)) {
    saveRDS(final_plot, save_path)
    
    
    # ---- Save to file ----
    ggplot2::ggsave(save_path2,
                    plot = final_plot,
                    device = "pdf",
                    width = width,
                    height = height,
                    dpi = dpi,
                    units = "in",
                    family = "Helvetica")
    
  }
  
  
  
  return(final_plot)
}

# ==========================================
# Plot the affine transformation of PCE risk
# ------------------------------------------
# This function illustrates how an affine transformation 
# (e.g., from a better-calibrated model) can shift PCE-predicted 
# cardiovascular risk. It compares predicted risk from the PCE model 
# to the "true" risk under a transformed scale.
# 
# Arguments:
# - slope: multiplicative adjustment to the linear predictor
# - intercept: additive shift to the linear predictor
# - save_path: optional file path to save the plot (e.g., "plot.pdf")
#
# Output: ggplot object (and optionally saves to file)
# ==========================================
plot_affine_risk_transform <- function(slope, intercept, save_path = NULL) {
  # Simulate predicted risk range
  risk_pce <- seq(0.001, 0.999, length.out = 500)
  
  # Transform to linear predictor scale (log-log link)
  lp_pce <- log(-log(1 - risk_pce))
  
  # Apply affine transformation
  lp_trans <- intercept + slope * lp_pce
  
  # Transform back to risk scale
  risk_true <- 1 - exp(-exp(lp_trans))
  
  # Prepare plot data
  df_plot <- data.frame(
    PCE_Risk  = risk_pce,
    True_Risk = risk_true,
    Group     = "Accent"
  )
  
  # Build plot
  p <- ggplot(df_plot, aes(x = PCE_Risk, y = True_Risk, color = Group)) +
    geom_line(linewidth = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = cb_palette["Accent"], guide = "none") +
    labs(x = "Predicted Risk", y = "Actual Risk") +
    theme_modern()
  
  # Optional: Save to file
  if (!is.null(save_path)) {
    saveRDS(p, save_path)
  }
  
  return(p)
}

# ---------------------------------------------------------------
# Assemble and save a panel of ggplot figures
#
# Arguments:
# - fig_matrix: matrix of filenames (each should be an .rds file containing a ggplot object)
# - filename:   name of the output EPS file (e.g., "figure_panel.eps")
# - width:      total width of output figure in inches
# - height:     total height of output figure in inches
# - dpi:        resolution in dots per inch (for raster elements; EPS is vector-based)
#
# This function:
# - Loads each plot from the specified .rds file
# - Assembles them into a panel layout (rows and columns)
# - Saves the panel as a high-resolution EPS file
# ---------------------------------------------------------------
save_panel_from_rds <- function(fig_matrix, filename,
                                width = 6, height = 6, dpi = 600) {
  # ---- Helper to fix font ----
  set_font <- function(p) {
    p + ggplot2::theme(text = ggplot2::element_text(family = "Helvetica"))
  }
  
  dims <- dim(fig_matrix)
  n_panels <- prod(dims)
  labels <- LETTERS[1:n_panels]
  
  # ---- Load plots row by row ----
  plot_list <- vector("list", n_panels)
  idx <- 1
  for (i in 1:dims[1]) {       # rows
    for (j in 1:dims[2]) {     # columns
      plot_list[[idx]] <- set_font(readRDS(fig_matrix[i, j]))
      idx <- idx + 1
    }
  }
  
  # ---- Wrap with manual tags ----
  panel_plot <- patchwork::wrap_plots(plotlist = plot_list,
                                      ncol = dims[2],
                                      byrow = TRUE) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold"))
  
  # Manually assign tags using patchwork's `tag_prefix` + `tag_suffix`
  panel_plot <- panel_plot +
    patchwork::plot_annotation(tag_levels = list(labels),
                               tag_prefix = "", tag_suffix = "")
  
  # ---- Save to file ----
  ggplot2::ggsave(filename,
                  plot = panel_plot,
                  device = "pdf",
                  width = width,
                  height = height,
                  dpi = dpi,
                  units = "in",
                  family = "Helvetica")
  
  invisible(panel_plot)
}