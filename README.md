# adaptive_rd  
Designing a quasi-experiment to study the clinical impact of adaptive risk prediction models

---

**View the Adaptive RD Simulation Demo:**  
https://cochran4.github.io/adaptive_rd/

---

## Overview

This repository accompanies the study:

**Odeh-Couvertier V., Zayas-Cabán G., Patterson B., Cochran A.**  
*Designing a quasi-experiment to study the clinical impact of adaptive risk prediction models* (Preprint, 2025).

The project develops and demonstrates a **regression discontinuity (RD)** framework for modern healthcare settings where both the **risk prediction model** and **decision threshold** can adapt over time. Traditional quasi-experimental designs assume fixed policies; here, we simulate adaptive policies and estimate local treatment effects when models are updated or thresholds are tuned (e.g., to target capacity, NNT, or calibration).

A Quarto-generated HTML demo illustrates the method, threshold dynamics, and estimation results.

---

## Repository Structure

- `adaptive_rd_demo.qmd` — Main Quarto document with the simulation demonstration and narrative.  
- `docs/` — GitHub Pages output (landing page `index.html`, demo HTML, assets).  
- `R/` — Core R source code (data I/O, prediction, threshold/model adaptation, estimation, plotting).  
- `figures/` — Saved figure objects and panels (`.rds`, `.pdf`) used in the demo/manuscript.  
- `data_raw/` *(ignored)* — Local-only NHANES `.XPT` files (not tracked).  
- `LICENSE` — MIT License.  
- `README.md` — This file.

> **Note:** Raw NHANES data is intentionally excluded. See below for setup.

---

## NHANES Data (Local Only — Not in Repo)

The simulation uses NHANES 2017–2018 public data to construct a cohort.  
Download the following `.XPT` files and place them in a local folder named `data_raw/` at the project root (already in `.gitignore`):

Required files:  
- `P_DEMO.XPT`  
- `P_BPXO.XPT`  
- `P_DIQ.XPT`  
- `P_BPQ.XPT`  
- `P_HDL.XPT`  
- `P_TCHOL.XPT`  
- `P_SMQ.XPT`  
- `P_MCQ.XPT` *(optional)*

Download from: https://wwwn.cdc.gov/nchs/nhanes/

**Example (R):**
```r
source("R/imports.R")
nh_list   <- read_nhanes_2017_2018("data_raw")
cohort_df <- build_cohort(nh_list, impute = TRUE)
```

---

## Core Components

- `simulate_design()` — Streams patients through blocks; assigns treatment; adapts threshold/model.  
- `adapt_threshold.R` — Quantile- and NNT-based adaptive rules.  
- `adapt_model.R` — Optional model recalibration/revision across blocks.  
- `estimate_spline()` — Spline-based estimator for the local ATE at the threshold.  
- `plot.R` — Plotting utilities (threshold trajectory, conditional means, comparisons).  
- `adaptive_rd_demo.qmd` — End-to-end demonstration (simulation → estimation → visualization).

---

## Reproduce the Demo Locally

```r
# 1) Install/load dependencies and project functions
source("R/imports.R")

# 2) Prepare cohort (requires NHANES files in data_raw/)
nh_list   <- read_nhanes_2017_2018("data_raw")
cohort_df <- build_cohort(nh_list)

# 3) Run a sample adaptive simulation
sim_out <- simulate_design(
  cohort_df              = cohort_df,
  risk_fn                = pce_predict,
  initial_block_size     = 400,
  block_size             = 100,
  n_blocks               = 26,
  threshold_adapt_method = list(type = "quantile", initial_threshold = 0.10, desired_treat_rate = 0.30),
  model_adapt_method     = list(type = "none"),
  outcome_fn             = outcome_attend
)

# 4) Estimate the local treatment effect at the (adaptive) threshold
fit_out <- estimate_spline(sim_out, family = binomial())
compare_ate_at_threshold(fit_out, sim_out, cohort_df, pce_predict, attendance_prob)
```

---

## Citing This Work

If you use this repository or build upon its methods, please cite:

**Odeh-Couvertier V., Zayas-Cabán G., Patterson B., Cochran A.**  
*Designing a quasi-experiment to study the clinical impact of adaptive risk prediction models.*  
Preprint, 2025.

---

## License

This project is released under the **MIT License**. See `LICENSE` for details.