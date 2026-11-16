# =========================
# Styling helpers (reusable)
# =========================

# Okabe–Ito palette (color-blind friendly)
cb_palette <- c(
  "Control" = "#666666", # Dark gray, strong but neutral
  "Treated" = "#0072B2", # Blue (instead of green, closer luminance)
  "Accent"  = "#D55E00" # Vermillion/orange-red
)

# Shared point, grid, and font settings.
POINT_ALPHA = .1
N_POINTS = 200
BASE_SIZE = 14
BASE_FAMILY = "sans"

# Shared figure theme; applied to each plot without changing the global theme.
theme_modern <- function() {
  ggplot2::theme_minimal(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "gray30"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "gray20"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = "gray90",
        linewidth = 0.4
      ),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "none",
      legend.title = ggplot2::element_text(size = BASE_SIZE, face = "bold"),
      legend.text = ggplot2::element_text(size = BASE_SIZE),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

# Plot finite thresholds in patient order; circles mark threshold changes.
# An optional save_path writes the ggplot object as RDS.
plot_threshold_trajectory <- function(sim_out, save_path = NULL) {
  thr <- sim_out$results %>%
    dplyr::arrange(block) %>%
    dplyr::mutate(
      patient_num = dplyr::row_number(), # running patient index
      # flag threshold jumps
      jump = threshold_used != dplyr::lag(threshold_used),
      jump = tidyr::replace_na(jump, TRUE) # first point is always a jump
    )

  # An infinite threshold denotes the initial block with no assignments.
  thr_plot <- dplyr::filter(thr, is.finite(threshold_used))

  p <- ggplot2::ggplot(thr_plot, ggplot2::aes(
    x = patient_num,
    y = threshold_used
  )) +
    ggplot2::geom_step(linewidth = 1.1, color = cb_palette["Accent"]) +
    ggplot2::geom_point(
      data = dplyr::filter(thr_plot, jump),
      size = 2.6, shape = 21, stroke = 0.9, fill = "white",
      color = cb_palette["Accent"]
    ) +
    ggplot2::labs(
      x = "Patient number (order seen)",
      y = "Referral threshold"
    ) +
    theme_modern()

  if (!is.null(save_path)) {
    saveRDS(p, save_path)
  }

  return(p)
}

# Plot fitted means, confidence bands, true means, and observed outcomes.
# fit_out comes from estimate_spline(); truth uses original PCE risk.
# x_limits controls the probability-risk axis; labels names the two assignments.
# Optional save_path supports RDS objects or PNG images.
plot_spline_curves <- function(
  fit_out, sim_out, truth_risk_fn, truth_mean_fn, cohort_df,
  ci = 0.95,
  y_lab = "Probability of attendance",
  labels = c("Standard letter", "Telephone invitation"),
  x_limits = c(0, 0.4),
  save_path = NULL
) {
  lbl_ctrl <- labels[1]
  lbl_trt <- labels[2]
  dat <- fit_out$data
  n_patients <- nrow(dat)

  # Get the last threshold used (for reference line)
  cutoff <- tail(sim_out$results$threshold_used, 1)

  # Compute plotting window
  zoom_min <- x_limits[1]
  zoom_max <- x_limits[2]

  # Generate sequence of x-values for smooth plotting
  x_seq <- seq(zoom_min, zoom_max, length.out = N_POINTS)

  # Outcomes always use the original model, including during model revision.
  original_risk <- truth_risk_fn(cohort_df)
  truth_mean <- function(cohort_df, idx, risk_baseline, treat) {
    truth_mean_fn(cohort_df, idx, original_risk[idx], treat)
  }
  cc_truth <- conditional_curves_truth(
    cohort_df, sim_out$risks_mat, truth_mean, x_seq
  )
  cc_truth$grp <- ifelse(cc_truth$group == "Control", lbl_ctrl, lbl_trt)

  cc <- conditional_curves(fit_out, x_seq, level = ci)
  cc$grp <- ifelse(cc$grp == "Control", lbl_ctrl, lbl_trt)

  # Split curves by side of threshold for alpha transparency
  s_ctrl <- subset(cc, grp == lbl_ctrl)
  s_trt <- subset(cc, grp == lbl_trt)
  s_ctrl_left <- subset(s_ctrl, risk_last <= cutoff); s_ctrl_left$alpha <- 1
  s_ctrl_right <- subset(s_ctrl, risk_last > cutoff); s_ctrl_right$alpha <- 0.7
  s_trt_left <- subset(s_trt, risk_last < cutoff); s_trt_left$alpha <- 0.7
  s_trt_right <- subset(s_trt, risk_last >= cutoff); s_trt_right$alpha <- 1
  ss <- rbind(s_ctrl_left, s_ctrl_right, s_trt_left, s_trt_right)

  obs <- subset(dat, risk_last >= zoom_min & risk_last <= zoom_max)
  obs$grp <- ifelse(obs$t == 0, lbl_ctrl, lbl_trt)

  # Set consistent factor levels
  lvl <- c(lbl_ctrl, lbl_trt)
  obs$grp <- factor(obs$grp, levels = lvl)
  ss$grp <- factor(ss$grp, levels = lvl)
  cc_truth$grp <- factor(cc_truth$grp, levels = lvl)

  accent_col <- unname(cb_palette["Accent"])

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
    mapping = ggplot2::aes(risk_last, mean,
      color = grp, alpha = alpha,
      linetype = "solid"
    ),
    linewidth = 1.1
  )

  # Reference line: final assignment threshold.
  p <- p + ggplot2::geom_vline(
    xintercept = cutoff, linetype = 2, linewidth = 0.9, color = accent_col
  )

  # Ground truth curves (dashed)
  p <- p + ggplot2::geom_line(
    data = cc_truth, linetype = "dashed",
    mapping = ggplot2::aes(risk_last, mean, color = grp),
    linewidth = 1.0, alpha = 0.9, show.legend = FALSE
  )

  p <- p + ggplot2::scale_color_manual(
    name = "", values = setNames(
      unname(cb_palette[c("Control", "Treated")]),
      lvl
    ),
    breaks = lvl, labels = lvl
  )

  p <- p + ggplot2::scale_fill_manual(
    values = setNames(unname(cb_palette[c("Control", "Treated")]), lvl),
    breaks = lvl, labels = lvl, guide = "none"
  )

  p <- p + ggplot2::scale_linetype_identity(guide = "none")
  p <- p + ggplot2::scale_alpha_identity()

  p <- p + ggplot2::coord_cartesian(xlim = c(zoom_min, zoom_max)) +
    ggplot2::labs(
      x = bquote(bold("Risk score, " * r[.(n_patients)])),
      y = y_lab
    ) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(face = "bold"))

  p <- p + ggplot2::guides(
    color = ggplot2::guide_legend(order = 1, override.aes = list(
      linewidth = 1.1
    ))
  )

  p <- p + theme_modern() +
    ggplot2::theme(legend.position = "bottom")

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

# plot_ldl_reduction_by_risk
# Shows the assignment-related LDL reduction specified by the outcome model and
# the risk at which it reaches the target reduction.
plot_ldl_reduction_by_risk <- function(
  target_reduction,
  save_path = NULL
) {
  risk <- seq(0, 1, by = 0.01)
  reduction_at <- function(r) {
    cholest_mean(NULL, seq_along(r), r, 0) -
      cholest_mean(NULL, seq_along(r), r, 1)
  }
  reduction <- reduction_at(risk)
  target_risk <- uniroot(function(r) reduction_at(r) - target_reduction,
    interval = c(0, 1)
  )$root

  df <- data.frame(risk = risk, reduction = reduction)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = risk, y = reduction)) +
    ggplot2::geom_line(color = cb_palette["Treated"], linewidth = 1.2) +
    ggplot2::geom_hline(
      yintercept = target_reduction,
      color = cb_palette["Accent"],
      linetype = "dashed",
      linewidth = 1
    ) +
    ggplot2::geom_vline(
      xintercept = target_risk,
      color = cb_palette["Accent"],
      linetype = "dashed",
      linewidth = 1
    ) +
    ggplot2::labs(
      x = "Predicted 10-year cardiovascular risk",
      y = "Additional LDL reduction\n(mg/dL)"
    ) +
    theme_modern()

  if (!is.null(save_path)) {
    if (grepl("\\.rds$", save_path, ignore.case = TRUE)) {
      saveRDS(p, save_path)
    } else if (grepl("\\.(png|pdf|svg)$", save_path, ignore.case = TRUE)) {
      ggplot2::ggsave(save_path, plot = p, width = 7, height = 4.5, dpi = 300)
    } else {
      warning(
        "save_path must end with .rds, .png, .pdf, or .svg to trigger saving."
      )
    }
  }

  p
}

# Plot risk histories for the first n_max sampled patients.
# Rows of risks_mat index source-cohort individuals; columns index block models.
# Optional save_path supports RDS, PNG, PDF, or SVG.
plot_risk_trajectory <- function(sim_out, n_max = 100, save_path = NULL) {
  R <- sim_out$risks_mat

  df_all <- as.data.frame(R)
  df_all$person <- seq_len(nrow(df_all))
  df_all <- tidyr::pivot_longer(df_all, -person,
    names_to = "block",
    values_to = "risk"
  )
  df_all$block <- as.integer(gsub("\\D", "", df_all$block)) # strip non-numeric

  df_spaghetti <- dplyr::filter(
    df_all,
    person %in% sim_out$results$person_id[1:min(n_max, nrow(sim_out$results))]
  )

  p <- ggplot2::ggplot() +

    # Spaghetti lines
    ggplot2::geom_line(
      data = df_spaghetti,
      aes(x = block, y = risk, group = person),
      color = "gray70", alpha = 0.5, linewidth = 0.4
    ) +

    # Axis labels and formatting
    ggplot2::labs(x = "Update", y = "Risk") +
    theme_modern()

  if (!is.null(save_path)) {
    if (grepl("\\.rds$", save_path, ignore.case = TRUE)) {
      saveRDS(p, save_path)
    } else if (grepl("\\.(png|pdf|svg)$", save_path, ignore.case = TRUE)) {
      ggsave(save_path, plot = p, width = 7, height = 4.5, dpi = 300)
    } else {
      warning(
        "save_path must end with .rds, .png, .pdf, or .svg to trigger saving."
      )
    }
  }

  return(p)
}

# ============================
# Plot the change for Black participants when White PCE equations are used
# ============================
plot_initial_revision_shift <- function(
  cohort_df,
  coef_original,
  coef_initial,
  save_path = NULL
) {
  original_risk <- pce_predict(cohort_df, coefs = coef_original)
  revised_risk <- pce_predict(
    cohort_df,
    coefs = coef_original,
    coef_adjust = coef_initial
  )

  group <- as.character(cohort_df$group)
  plot_data <- data.frame(
    original_risk = original_risk,
    revised_risk = revised_risk,
    race = sub("_[FM]$", "", group),
    sex = ifelse(grepl("_F$", group), "Women", "Men")
  )

  plot_data <- plot_data[plot_data$race == "Black", ]

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = original_risk,
      y = revised_risk
    )
  ) +
    ggplot2::geom_abline(
      ggplot2::aes(slope = 1, intercept = 0, linetype = "No change"),
      color = "gray40"
    ) +
    ggplot2::geom_point(
      color = cb_palette[["Treated"]], alpha = 0.45,
      size = 1
    ) +
    ggplot2::facet_wrap(~sex) +
    ggplot2::scale_linetype_manual(
      name = NULL, values = c("No change" = "dashed")
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 0.4), ylim = c(0, 0.4)) +
    ggplot2::labs(
      x = "Original PCE risk",
      y = "Risk using White PCE equation"
    ) +
    theme_modern() +
    ggplot2::theme(legend.position = "bottom")

  if (!is.null(save_path)) {
    saveRDS(p, save_path)
  }

  p
}

# Four equal-sized error panels with replication points and a bottom legend.
# S1, S3, S4 use pp; S2 uses mg/dL. S3 and S4 share limits.
# Optional label vectors set display names/order; save_path stores RDS and
# save_path2 stores PDF. width, height, and dpi control the saved figure.
plot_comparison <- function(mc_df,
                            method_labels = NULL,
                            scenario_labels = NULL,
                            width = 16,
                            height = 10,
                            dpi = 600,
                            ylim_list = list(
                              c(-35, 30), # Percentage points, Scenario 1
                              c(-23, 15), # mg/dL, Scenario 2
                              # Shared percentage-point limits, Scenarios 3 and
                              # 4
                              c(-8, 14)
                            ),
                            save_path = NULL,
                            save_path2 = NULL) {
  # Accent color for spline, grays for others
  color_palette <- c(
    spline = cb_palette[["Accent"]],
    naive = "gray15",
    outreg = "gray40",
    outreg_separate = "gray50",
    ipw = "gray65",
    aipw = "gray90"
  )

  # Optional patterns to differentiate comparators
  pattern_palette <- c(
    "Adaptive RD" = "none",
    "Outcome regression (pooled)" = "none",
    "Outcome regression (separate arms)" = "none",
    "Naive" = "none",
    "Inverse Probability Weighting" = "none",
    "Augmented Inverse Probability Weighting" = "none"
  )

  mc_df <- mc_df %>%
    mutate(est_error = (est - truth) * if_else(scenario == "S2", 1, 100))

  # Save original scenario identifiers, then apply display labels once.
  scenario_id <- mc_df$scenario
  if (!is.null(method_labels)) {
    mc_df$method <- factor(mc_df$method,
      levels = names(method_labels), labels = method_labels
    )
    color_palette <- setNames(
      color_palette[names(method_labels)], method_labels
    )
  }
  if (!is.null(scenario_labels)) {
    mc_df$scenario <- factor(mc_df$scenario,
      levels = names(scenario_labels), labels = scenario_labels
    )
  }

  # Split the data by outcome scale, retaining identical method factor levels.
  df1 <- mc_df[scenario_id == "S1", , drop = FALSE]
  df2 <- mc_df[scenario_id == "S2", , drop = FALSE]
  df34 <- mc_df[scenario_id %in% c("S3", "S4"), , drop = FALSE]

  # Shared boxplot and replication-point styling for all four panels.
  base_plot <- function(df, ylim = NULL, y_label = "Estimation error (pp)") {
    p <- ggplot(df, aes(
      x = method, y = est_error, fill = method,
      color = method, pattern = method
    )) +
      geom_point(
        position = position_jitter(width = 0.16, height = 0, seed = 1234),
        alpha = 0.14, size = 0.5, show.legend = FALSE
      ) +
      ggpattern::geom_boxplot_pattern(
        outlier.shape = NA,
        color = "black",
        alpha = 0.55,
        pattern_density = 0.5,
        pattern_spacing = 0.04,
        pattern_key_scale_factor = 0.6,
        linewidth = 0.3
      ) +
      geom_hline(
        yintercept = 0, linetype = "dashed", color = "black",
        linewidth = 0.6
      ) +
      scale_fill_manual(values = color_palette) +
      scale_color_manual(values = color_palette) +
      scale_pattern_manual(values = pattern_palette) +
      labs(y = y_label, x = NULL) +
      theme_modern() +
      theme(
        strip.background = element_rect(fill = "gray90", color = "gray50"),
        strip.text = element_text(face = "bold", size = 18),
        panel.border = element_rect(color = "gray50", fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )

    if (!is.null(ylim)) {
      p <- p + coord_cartesian(ylim = ylim)
    }

    return(p)
  }

  legend_plot <- ggplot(mc_df, aes(
    x = method, y = est_error, fill = method,
    color = method, pattern = method
  )) +
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
      legend.text = element_text(size = 14)
    )

  # Equal-size panels, with a shared legend below the grid.
  p1 <- base_plot(df1, ylim = ylim_list[[1]]) +
    facet_wrap(~scenario, ncol = 1) + theme(legend.position = "none")
  p2 <- base_plot(df2,
    ylim = ylim_list[[2]],
    y_label = "Estimation error (mg/dL)"
  ) +
    facet_wrap(~scenario, ncol = 1) + theme(legend.position = "none")
  p3 <- base_plot(filter(df34, scenario == unique(df34$scenario)[1]),
    ylim = ylim_list[[3]]
  ) +
    facet_wrap(~scenario, ncol = 1) + theme(legend.position = "none")
  p4 <- base_plot(filter(df34, scenario == unique(df34$scenario)[2]),
    ylim = ylim_list[[3]]
  ) +
    facet_wrap(~scenario, ncol = 1) + theme(legend.position = "none")

  legend <- cowplot::get_plot_component(
    legend_plot + theme(legend.position = "bottom") +
      guides(
        fill = guide_legend(nrow = 2, byrow = TRUE),
        color = guide_legend(nrow = 2, byrow = TRUE),
        pattern = guide_legend(nrow = 2, byrow = TRUE)
      ),
    "guide-box-bottom"
  )
  panels <- (p1 + p2) / (p3 + p4)
  final_plot <- patchwork::wrap_elements(panels) / patchwork::wrap_elements(
    legend
  ) +
    plot_layout(heights = c(1, 0.17))

  if (!is.null(save_path)) {
    saveRDS(final_plot, save_path)

      ggplot2::ggsave(save_path2,
      plot = final_plot,
      device = "pdf",
      width = width,
      height = height,
      dpi = dpi,
      units = "in",
      family = "Helvetica"
    )
  }

  return(final_plot)
}

# Compare PCE probability risk with an affine adjustment on the cloglog scale.
# The dashed diagonal marks no change; optional save_path writes an RDS plot.
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
    geom_abline(
      slope = 1, intercept = 0, linetype = "dashed",
      color = "gray50"
    ) +
    scale_color_manual(values = cb_palette["Accent"], guide = "none") +
    labs(x = "Predicted Risk", y = "Actual Risk") +
    theme_modern()

  # Optional: Save to file
  if (!is.null(save_path)) {
    saveRDS(p, save_path)
  }

  return(p)
}

# Load saved plots, apply the paper font, assemble panels, and save a PDF.
# Use rows for unequal panel counts per row, design for a patchwork layout,
# or a filename matrix for a regular grid. Dimensions are in inches.
save_panel_from_rds <- function(fig_matrix, filename,
                                width = 6, height = 6, dpi = 600,
                                design = NULL, rows = NULL) {
  # Apply the same font and panel-tag styling to each saved plot.
  set_font <- function(p) {
    p + ggplot2::theme(
      text = ggplot2::element_text(family = "Helvetica"),
      plot.tag.position = "topleft",
      plot.tag = ggplot2::element_text(
        face = "bold",
        hjust = 0,
        vjust = 1
      )
    )
  }

  dims <- dim(fig_matrix)
  n_panels <- length(fig_matrix)
  labels <- LETTERS[1:n_panels]

  # Load plots in display order.
  paths <- if (is.null(dims)) fig_matrix else as.vector(t(fig_matrix))
  plot_list <- lapply(paths, function(path) set_font(readRDS(path)))

  # Build rows independently when their panel counts differ.
  if (!is.null(rows)) {
    plot_list <- Map(
      function(p, label) p + ggplot2::labs(tag = label),
      plot_list, labels
    )
    row_plots <- lapply(rows, function(i) {
      patchwork::wrap_plots(
        plotlist = plot_list[i], nrow = 1,
        widths = rep(1, length(i))
      )
    })
    panel_plot <- patchwork::wrap_plots(
      plotlist = row_plots, ncol = 1, heights = rep(1, length(rows))
    )
  } else {
    panel_plot <- if (is.null(design)) {
      patchwork::wrap_plots(plotlist = plot_list, ncol = dims[2], byrow = TRUE)
    } else {
      patchwork::wrap_plots(plotlist = plot_list, design = design)
    }
    panel_plot <- panel_plot + patchwork::plot_annotation(
      tag_levels = list(labels), tag_prefix = "", tag_suffix = ""
    )
  }

  ggplot2::ggsave(filename,
    plot = panel_plot,
    device = "pdf",
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    family = "Helvetica"
  )

  invisible(panel_plot)
}

# Rebuild the paper comparison figure from saved Monte Carlo results.
plot_comparison_paper <- function(mc) {
  plot_comparison(
    mc,
    ylim_list = list(c(-35, 30), c(-23, 15), c(-8, 14)),
    method_labels = c(
      spline = "Adaptive RD", naive = "Naive",
      outreg = "Outcome regression (pooled)",
      outreg_separate = "Outcome regression (separate arms)",
      ipw = "Inverse Probability Weighting",
      aipw = "Augmented Inverse Probability Weighting"
    ),
    scenario_labels = c(
      S1 = "Scenario 1", S2 = "Scenario 2",
      S3 = "Scenario 3", S4 = "Scenario 4"
    ),
    save_path = "figures/fig_comparison.rds",
    save_path2 = "figures/fig_comparison.pdf"
  )
}
