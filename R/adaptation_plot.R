# Build the paper's six-panel adaptation figure from saved results.
plot_adaptation_paper <- function(
    results = readRDS("results/adaptation_paper/adaptation_results.rds"),
    output_dir = "figures",
    filename = "fig_adaptation"
) {
  # Measure panel dimensions on a temporary device before aligning them.
  grDevices::pdf(NULL)
  device <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(device))

  x <- results

  # Convert setting names to ramp lengths or penalty values.
  value <- function(setting) {
    as.numeric(sub(".*_", "", ifelse(setting == "fixed", "0", setting)))
  }
  ramp_levels <- sort(unique(value(
    x$thresholds$setting[x$thresholds$scenario == "S1"]
  )))
  penalty_levels <- sort(unique(value(x$pca$setting[x$pca$scenario == "S4"])))

  # Share setting labels, line styles, and point shapes across panels.
  ramp_labels <- c("Fixed", "3-block ramp", "6-block ramp")[
    match(ramp_levels, c(0, 3, 6))
  ]
  ramp_type <- c(`0` = "solid", `3` = "dotdash", `6` = "longdash")[
    as.character(ramp_levels)
  ]
  penalty_type <- stats::setNames(c("solid", "dashed", "dotdash"),
    penalty_levels
  )
  penalty_shape <- stats::setNames(c(16, 17, 15), penalty_levels)
  ramp_scale <- ggplot2::scale_linetype_manual(
    values = ramp_type, name = "Threshold rule", labels = ramp_labels
  )
  penalty_scale <- ggplot2::scale_linetype_manual(
    values = penalty_type, name = "Revision penalty"
  )

  # Use the paper's theme, with room for panel letters and bottom legends.
  base <- theme_modern() + ggplot2::theme(
    text = ggplot2::element_text(family = "Helvetica"),
    legend.position = "bottom",
    plot.tag.position = c(0.12, 0.98),
    plot.tag = ggplot2::element_text(face = "bold", hjust = 0, vjust = 1)
  )

  # Panels C and F share symmetric error limits, with 10% headroom.
  plotted_errors <- x$performance$error[
    x$performance$scenario %in% c("S1", "S4")
  ]
  error_bound <- max(
    0.01, ceiling(110 * max(abs(plotted_errors), na.rm = TRUE)) / 100
  )
  error_limits <- c(-error_bound, error_bound)

  # Both error panels show replication points over boxplots, in pp.
  error_panel <- function(scenario, levels, label) {
    data <- x$performance |>
      dplyr::filter(.data$scenario == .env$scenario) |>
      dplyr::mutate(
        setting_value = factor(value(setting), levels = levels),
        error = 100 * error
      )
    ggplot2::ggplot(data, ggplot2::aes(setting_value, error)) +
      ggplot2::geom_hline(
        yintercept = 0, linetype = "dashed", color = "gray50"
      ) +
      ggplot2::geom_boxplot(
        outlier.shape = NA, fill = "gray90",
        color = cb_palette["Control"], linewidth = 0.6
      ) +
      ggplot2::geom_point(
        position = ggplot2::position_jitter(width = 0.07, seed = 1),
        size = 1.6, alpha = 0.5, color = cb_palette["Control"]
      ) +
      ggplot2::coord_cartesian(ylim = 100 * error_limits) +
      ggplot2::labs(x = label, y = "Estimation error (pp)") +
      base
  }

  # A: Mean threshold by block. Omit the infinite startup threshold.
  thresholds <- x$thresholds |>
    dplyr::filter(scenario == "S1", is.finite(threshold)) |>
    dplyr::mutate(ramp = factor(value(setting), levels = ramp_levels)) |>
    dplyr::group_by(ramp, block) |>
    dplyr::summarise(threshold = mean(threshold), .groups = "drop") |>
    dplyr::mutate(patient_num = (block - 1) * x$config$block_size + 1)
  a <- ggplot2::ggplot(
    thresholds, ggplot2::aes(patient_num, threshold, linetype = ramp)
  ) +
    ggplot2::geom_step(linewidth = 1.1, color = cb_palette["Accent"]) +
    ramp_scale +
    ggplot2::labs(
      x = "Patient number (order seen)", y = "Mean risk threshold"
    ) +
    base

  # B: Center risk bins on each run's final cutoff. Average assignment rates
  # within runs first, then across runs, giving each run equal weight.
  assignments <- x$assignments |>
    dplyr::filter(scenario == "S1") |>
    dplyr::mutate(
      ramp = factor(value(setting), levels = ramp_levels),
      centered_risk = floor((mean_risk - final_threshold) / 0.01) * 0.01 + 0.005
    ) |>
    dplyr::group_by(ramp, replication, centered_risk) |>
    dplyr::summarise(
      assignment_proportion = sum(n_assigned) / sum(n_patients),
      .groups = "drop"
    ) |>
    dplyr::group_by(ramp, centered_risk) |>
    dplyr::summarise(proportion = mean(assignment_proportion), .groups = "drop")
  b <- ggplot2::ggplot(
    assignments,
    ggplot2::aes(centered_risk, 100 * proportion,
      linetype = ramp, group = ramp
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed", color = "gray50"
    ) +
    ggplot2::geom_line(linewidth = 1.1, color = cb_palette["Accent"]) +
    ramp_scale +
    ggplot2::scale_y_continuous(limits = c(0, 100)) +
    ggplot2::coord_cartesian(xlim = c(-0.2, 0.2)) +
    ggplot2::labs(
      x = "Final risk score minus final threshold",
      y = "Mean assignment rate (%)"
    ) +
    base

  # C: Estimation errors under the fixed and ramped threshold rules.
  c <- error_panel("S1", ramp_levels, "Threshold rule") +
    ggplot2::scale_x_discrete(labels = c(
      `0` = "Fixed", `3` = "3-block ramp", `6` = "6-block ramp"
    ))

  # D: Mean absolute change from the starting model, across replications.
  trajectories <- dplyr::bind_rows(x$risk_changes) |>
    dplyr::mutate(penalty = factor(value(setting), levels = penalty_levels)) |>
    dplyr::group_by(penalty, block) |>
    dplyr::summarise(change = mean(mean_absolute_change), .groups = "drop") |>
    dplyr::mutate(patient_num = (block - 1) * x$config$block_size + 1)
  d <- ggplot2::ggplot(
    trajectories, ggplot2::aes(patient_num, 100 * change, linetype = penalty)
  ) +
    ggplot2::geom_step(linewidth = 1.1, color = cb_palette["Treated"]) +
    penalty_scale +
    ggplot2::labs(
      x = "Patient number (order seen)",
      y = "Mean absolute risk change\nfrom starting model (pp)"
    ) +
    base

  # E: Mean cumulative variance of the first five residual-history PCs.
  pca <- x$pca |>
    dplyr::filter(scenario == "S4") |>
    dplyr::mutate(penalty = factor(value(setting), levels = penalty_levels)) |>
    dplyr::group_by(penalty, component) |>
    dplyr::summarise(
      cumulative = mean(cumulative_variance, na.rm = TRUE), .groups = "drop"
    )
  e <- ggplot2::ggplot(
    pca, ggplot2::aes(component, 100 * cumulative,
      linetype = penalty, shape = penalty
    )
  ) +
    ggplot2::geom_line(linewidth = 1.1, color = cb_palette["Treated"]) +
    ggplot2::geom_point(size = 2.6, color = cb_palette["Treated"]) +
    ggplot2::scale_x_continuous(breaks = 1:5) +
    penalty_scale +
    ggplot2::scale_shape_manual(
      values = penalty_shape, name = "Revision penalty"
    ) +
    ggplot2::labs(
      x = "Number of principal components",
      y = "Mean cumulative variance explained (%)"
    ) +
    base

  # F: Estimation errors under the three revision penalties.
  f <- error_panel("S4", penalty_levels, "Revision penalty")

  # Extract one bottom legend for each scenario before hiding panel legends.
  bottom_legend <- function(plot) {
    grob <- ggplot2::ggplotGrob(plot)
    index <- which(grob$layout$name == "guide-box-bottom")
    grob$grobs[[index[1L]]]
  }
  ramp_legend <- bottom_legend(a)
  revision_legend <- bottom_legend(e)

  # Align all six panels, then arrange two rows of three equally sized panels.
  panels <- lapply(list(a, b, c, d, e, f), function(plot) {
    plot + ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 24, r = 8, b = 8, l = 8)
    )
  })
  panels <- cowplot::align_plots(plotlist = panels, align = "hv", axis = "tblr")
  top <- cowplot::plot_grid(
    plotlist = panels[1:3], nrow = 1, labels = c("A", "B", "C"),
    label_size = BASE_SIZE, label_fontfamily = "Helvetica",
    label_x = 0.04, label_y = 0.99, align = "hv", axis = "tblr"
  )
  bottom <- cowplot::plot_grid(
    plotlist = panels[4:6], nrow = 1, labels = c("D", "E", "F"),
    label_size = BASE_SIZE, label_fontfamily = "Helvetica",
    label_x = 0.04, label_y = 0.99, align = "hv", axis = "tblr"
  )
  legends <- cowplot::plot_grid(ramp_legend, revision_legend, nrow = 1)
  figure <- cowplot::plot_grid(
    top, bottom, legends, ncol = 1, rel_heights = c(1, 1, 0.12)
  )

  # Save the plot object and the PDF/PNG versions used by the paper and Quarto.
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(figure, file.path(output_dir, paste0(filename, ".rds")))
  ggplot2::ggsave(
    file.path(output_dir, paste0(filename, ".pdf")), figure,
    width = 16, height = 10, family = "Helvetica"
  )
  ggplot2::ggsave(
    file.path(output_dir, paste0(filename, ".png")), figure,
    width = 16, height = 10, dpi = 150, bg = "white"
  )
  invisible(figure)
}
