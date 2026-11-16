# Adaptive regression discontinuity

**View the Adaptive RD Simulation Demo:**
[cochran4.github.io/adaptive_rd](https://cochran4.github.io/adaptive_rd/)

## Overview

This repository accompanies the study by **Odeh-Couvertier V., Zayas-Cabán G.,
Patterson B., and Cochran A.**:

*Designing a quasi-experiment to study the clinical impact of risk prediction
models that adapt after deployment.*

The project develops and demonstrates a **regression discontinuity (RD)**
framework for healthcare settings where both the **risk prediction model** and
**decision threshold** can adapt over time. The simulations evaluate local
assignment effects while thresholds target capacity or estimated benefit, and
prediction models are recalibrated or revised as patient data accumulate.

A Quarto-generated HTML demonstration explains the scenarios, threshold and
model updates, estimation, and results.

## Repository structure

- `adaptive_rd_demo.qmd` — Quarto demonstration and narrative.
- `docs/` — GitHub Pages landing page, rendered demonstration, and assets.
- `R/` — Data preparation, prediction, adaptation, estimation, and plotting.
- `run_comparison.R` — Estimator-comparison simulation runner.
- `run_adaptation.R` — Adaptation-study simulation runner.
- `results/` — Final saved simulation results, organized by analysis.
- `figures/` — Publication figures and saved figure objects.
- `references.bib` — References cited in the demonstration.
- `data_raw/` *(ignored)* — Local NHANES `.XPT` files, not tracked by Git.
- `LICENSE` — MIT License.
- `README.md` — Overview and instructions.

> Raw NHANES data is intentionally excluded. See setup below.

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
