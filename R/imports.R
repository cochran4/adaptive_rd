# imports.R

# Toggle: install missing packages automatically
INSTALL_MISSING <- TRUE

pkgs <- c(
  "dplyr",
  "ggplot2",
  "ggpattern",
  "purrr",
  "readr",
  "foreign",  
  "mice",
  "tibble",
  "knitr",
  "skimr",
  "scales",
  "future.apply",
  "progressr",
  "patchwork",
  "cowplot"
)


# Load ggpp separately to suppress S3 overwrite warning
suppressWarnings(suppressPackageStartupMessages(library(ggpp))) 

# Pick a CRAN mirror if none set (avoids interactive prompt)
if (isTRUE(INSTALL_MISSING) && (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

# Install any missing packages (quietly)
if (INSTALL_MISSING) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    install.packages(missing_pkgs, dependencies = TRUE)
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
source("R/adapt_threshold.R")   # adapt_threshold_quantile(), adapt_threshold_nnt(), adapt_threshold()

# Model adaptation methods
source("R/adapt_model.R")   # adapt_model_revise(), adapt_model_recalibrate(), adapt_model() 

# Estimation methods
source("R/estimate.R")   # estimate_spline()

# Main simulation
source("R/simulate.R")   # adapt_threshold_quantile(), adapt_threshold_nnt(), adapt_threshold()

# Plotting functions
source("R/plot.R")   # plot_spline_curves()

# Tables
source("R/tables.R")   # plot_spline_curves()

# Functions needed for simulation scenarios
source("R/scenario_functions.R")   # plot_spline_curves()