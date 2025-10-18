# adaptive_rd
Designing a quasi-experiment to study the clinical impact of adaptive risk prediction models

---

## Overview
This repository contains code, analyses, and documentation used to design and evaluate quasi‑experimental studies that measure the clinical impact of adaptive risk prediction models (models that are periodically retrained/updated in production). The repository includes data‑preparation scripts, simulation code, model training/evaluation code, and a Quarto document with rendered HTML that walk through the simulations reported in the associated paper.

Project site (Quarto): https://cochran4.github.io/adaptive_rd

## Repository layout (expected)
- README.md — this file
- data_processed/ — outputs produced by preprocessing scripts (do not commit large derived files)
- scripts/ or R/ or src/ — preprocessing, model training, and simulation scripts
- analysis/quarto/ — Quarto source (.qmd) files used to generate the site and HTML
- docs/ or analysis/quarto/_site/ — rendered HTML pages for the Quarto site
- results/ — generated figures, tables, and model artifacts
- configs/ — YAML/JSON configuration files for experiments
- environment.yml or requirements.txt — environment specification for reproducibility
- LICENSE — MIT license (present in this repository)

## Handling NHANES data (do not add raw data to the repo)
NHANES data files must not be committed to the public repository. Instead, place the required NHANES files in a local directory outside the repository and point preprocessing scripts to that directory via a config value or environment variable.

Required NHANES files (exact filenames expected by preprocessing):
- P_DEMO.XPT
- P_BPXO.XPT
- P_DIQ.XPT
- P_BPQ.XPT
- P_HDL.XPT
- P_MCQ.XPT    (optional: not used by PCE but safe to include)
- P_SMQ.XPT
- P_TCHOL.XPT

Where to download:
- NHANES data portal: https://wwwn.cdc.gov/nchs/nhanes/
  - Select the relevant cycle(s) (e.g., 2017‑2018, 2019‑2020, or combined releases) and download the component files above (SAS transport .XPT files).

Recommended pattern: configure a local NHANES directory and add it to your run-time configuration rather than adding files to the repo.

Environment variable approach (example)
- Set the path to your local NHANES directory:
  - macOS / Linux / WSL:
    export NHANES_DIR="/path/to/local/nhanes_files"
  - Windows (PowerShell):
    $env:NHANES_DIR = "C:\path\to\local\nhanes_files"

R example (reading files from NHANES_DIR)

```r
nhanes_dir <- Sys.getenv("NHANES_DIR")
demo   <- file.path(nhanes_dir, "P_DEMO.XPT")
bpx    <- file.path(nhanes_dir, "P_BPXO.XPT")
diq    <- file.path(nhanes_dir, "P_DIQ.XPT")
bpq    <- file.path(nhanes_dir, "P_BPQ.XPT")
hdl    <- file.path(nhanes_dir, "P_HDL.XPT")
mcq    <- file.path(nhanes_dir, "P_MCQ.XPT")
smq    <- file.path(nhanes_dir, "P_SMQ.XPT")
tchol  <- file.path(nhanes_dir, "P_TCHOL.XPT")

# e.g. read one file
library(haven)
demo_df <- read_xpt(demo)
```

Python example (reading files from NHANES_DIR)

```python
import os
import pyreadstat

nhanes_dir = os.environ.get("NHANES_DIR")
demo_path = os.path.join(nhanes_dir, "P_DEMO.XPT")
demo_df, meta = pyreadstat.read_xport(demo_path)
```

If you convert XPT to CSV for local use, keep consistent filenames and update your config to point to the CSV versions.

## Preprocessing and pipeline (example commands)
- Typical R preprocessing (example - replace with your actual script name and flags):
  ```
  Rscript scripts/prep_nhanes.R --input-dir "$NHANES_DIR" --output-dir data_processed
  ```
- Typical Python preprocessing (example):
  ```
  python scripts/prep_nhanes.py --input-dir "$NHANES_DIR" --output data_processed
  ```

## Rendering the Quarto site and reproducing simulations
- Install Quarto: https://quarto.org/docs/get-started/
- From repository root, render the site or specific pages:
  ```
  quarto render analysis/quarto --to html
  ```
  or
  ```
  quarto render analysis/quarto/simulations.qmd --to html
  ```
- Open the generated HTML in `analysis/quarto/_site/` or `docs/`.

The published site is available at: https://cochran4.github.io/adaptive_rd

## Paper and citation
This repository documents the simulations and analyses for the manuscript "Designing a quasi-experiment to study the clinical impact of adaptive risk prediction models." Please cite the manuscript and the project site when reusing methods or figures. Add a BibTeX entry or CITATION file if available.

## License — MIT
This repository is distributed under the MIT License (a copy of the LICENSE file is included in the repo). Summary:
- You may use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the software.
- You must include the original copyright and license notice in any substantial portions of the software.
- The software is provided "as is", without warranty; authors are not liable for damages.

## Reproducibility & privacy
- Pin package versions in `environment.yml` or `requirements.txt`.
- Store deterministic seeds and experiment configs in `configs/`.
- Do not commit PHI or other sensitive data. Keep NHANES files outside the repo and reference them via an environment variable or configuration file.

## Contact / maintainer
Maintainer: cochran4  
For questions or contributions, open an issue or submit a pull request on GitHub.
```
