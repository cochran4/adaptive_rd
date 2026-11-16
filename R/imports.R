# imports.R

# Toggle: install missing packages automatically
INSTALL_MISSING <- TRUE

pkgs <- c(
  "dplyr",
  "ggplot2",
  "ggpattern",
  "purrr",
  "foreign",
  "mice",
  "tibble",
  "knitr",
  "future.apply",
  "progressr",
  "patchwork",
  "cowplot"
)

# Pick a CRAN mirror if none set (avoids interactive prompt)
if (isTRUE(INSTALL_MISSING) &&
    (is.null(getOption("repos")) ||
     getOption("repos")["CRAN"] == "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

# Install any missing packages (quietly)
if (INSTALL_MISSING) {
  missing_pkgs <- pkgs[
    !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_pkgs)) {
    install.packages(missing_pkgs)
  }
}

# Load packages without startup chatter
suppressPackageStartupMessages({
  invisible(lapply(pkgs, library, character.only = TRUE))
})

invisible(TRUE)

# Data processing and risk calculation
source("R/data_io.R")           # read_nhanes_2017_2018(), build_cohort()
source("R/pce.R")               # pce_predict(), pce_coefs

# Threshold adaptation methods
# adapt_threshold_quantile(), adapt_threshold_benefit(), adapt_threshold()
source("R/adapt_threshold.R")

# Model adaptation methods
# adapt_model_revise(), adapt_model_recalibrate(), adapt_model()
source("R/adapt_model.R")

# Estimation methods
source("R/estimate.R")   # estimate_spline()

# Main simulation
source("R/simulate.R")   # simulate_design()

# Plotting functions
source("R/plot.R")   # plot_spline_curves()

# Tables
source("R/tables.R")   # compare_ate_at_threshold()

# Functions needed for simulation scenarios
source("R/scenario_functions.R")   # shared adaptation settings and outcomes
