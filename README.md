# Adaptive regression discontinuity

Simulation code accompanying *Designing a quasi-experiment to study the clinical
impact of risk prediction models that adapt after deployment*.

[View the demonstration](https://cochran4.github.io/adaptive_rd/).

## Setup

Download the NHANES 2017–March 2020 pre-pandemic files from
https://wwwn.cdc.gov/nchs/nhanes/ into `data_raw/`:

`P_DEMO.XPT`, `P_BPXO.XPT`, `P_DIQ.XPT`, `P_BPQ.XPT`, `P_HDL.XPT`,
`P_MCQ.XPT`, `P_SMQ.XPT`, `P_TCHOL.XPT`.

Run commands from the repository root. `R/imports.R` loads the shared functions
and installs missing R packages. Quarto is required to render the demonstration.

## Run the simulations

```r
source("run_comparison.R")  # 2,000 replications per scenario
source("run_adaptation.R")  # 500 per setting; reuses completed settings
```

Results are saved in `results/comparison/` and `results/adaptation_paper/`.
Shared scenario definitions are in `R/scenario_functions.R`.

## Rebuild figures from saved results

```r
source("R/imports.R")
source("R/adaptation_plot.R")
plot_comparison_paper(readRDS("results/comparison/comparison_results.rds"))
plot_adaptation_paper(readRDS("results/adaptation_paper/adaptation_results.rds"))
```

Render the four-scenario demonstration and its two illustrative paper figures:

```sh
quarto render adaptive_rd_demo.qmd
```

Rendering runs one example per scenario and reads the Monte Carlo figures.
It does not rerun the Monte Carlo studies. The four paper figures are in
`figures/`; the rendered website is in `docs/`.

MIT License. See `LICENSE`.
