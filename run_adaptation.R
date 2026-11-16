# Run the adaptation simulation study from the project root.
#
# From RStudio:
#   source("run_adaptation.R")
#
# From a terminal, optionally supplying the number of replications:
#   Rscript run_adaptation.R 500

source("R/imports.R")
source("R/adaptation_runs.R")
source("R/adaptation_plot.R")

arguments <- commandArgs(trailingOnly = TRUE)
n_replications <- if (length(arguments) >= 1L) {
  as.integer(arguments[[1]])
} else {
  500L
}

if (is.na(n_replications) || n_replications < 1L) {
  stop("The number of replications must be a positive integer.")
}

# S4 is the model-revision scenario.
adaptation_results <- run_adaptation(
  n_replications = n_replications,
  scenarios = c("S1", "S4")
)
plot_adaptation_paper(adaptation_results)
