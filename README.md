# Abdominal Perfusion Subphenotypes in Acute Pancreatitis — Analysis Code

R code for the main analyses of a study identifying longitudinal
mean-arterial-pressure (MAP) and intra-abdominal-pressure (IAP) trajectory
subphenotypes in acute pancreatitis, and evaluating their prognostic value and
differential response to fluid resuscitation.

## Overview

The pipeline has two parts:

**Phase 1 — Subphenotype derivation and prognosis**
- Group-based multivariate trajectory modelling (GBMT) of joint MAP/IAP
  trajectories, yielding four subphenotypes.
- Kaplan–Meier survival and Cox proportional-hazards models for 28-day mortality.
- External validation of the trajectory subphenotypes in an independent cohort.

**Phase 2 — Fluid resuscitation response**
- Descriptive analysis of daily/cumulative fluid volume by subphenotype.
- Linear mixed-effects (LME) models of the fluid × subphenotype interaction on
  next-day change in MAP and IAP.
- Bayesian multivariate mixed-effects model (brms/Stan) as a sensitivity check.

## Repository layout

```
.
├── R/                                  # analysis scripts (run in numeric order)
│   ├── 01_trajectory_gbmt.R            # GBMT, 2–5 groups, fit indices, plots
│   ├── 02_kaplan_meier.R               # KM curves + log-rank (28-day mortality)
│   ├── 03_cox_nested_models.R          # Cox nested models + VIF/PH diagnostics
│   ├── 04_cox_forest_plot.R            # forest plots (reads output of 03)
│   ├── 05_external_validation_gbmt.R   # MICE + GBMT in external cohort
│   ├── 06_phase2_vfr_descriptive.R     # fluid volume descriptives
│   ├── 07_phase2_lme.R                 # LME interaction models + phase space
│   ├── 08_phase2_bayesian_fit.R        # brms model fit (long-running)
│   ├── 09_phase2_bayesian_summary.R    # brms diagnostics & posterior effects
│   ├── 10_sensitivity_gbmt_complete_case.R
│   └── 11_sensitivity_lme_complete_case.R
├── data/                               # input data (not provided — see data/README.md)
├── models/                             # cached fitted models (created at runtime)
├── output/                             # results, tables, figures (created at runtime)
├── install_packages.R                  # one-off dependency installer
├── sessionInfo.txt                     # exact R/OS/package versions used
├── packages.tsv                        # tab-separated package version table
├── CITATION.cff                        # citation metadata (fill in on publication)
├── LICENSE                             # MIT
└── README.md
```

## Requirements

- R (≥ 4.3 recommended)
- Script 08 additionally requires a working C++ toolchain and Stan
  (via the `brms` / `rstan` or `cmdstanr` backend).

Install all CRAN dependencies:

```r
source("install_packages.R")
```

## Running

Scripts read from `data/` and write to `output/` (and `models/`) using paths
relative to the repository root, so run R with the working directory set to the
repository root:

```r
setwd("/path/to/this/repo")
source("R/01_trajectory_gbmt.R")
```

Order matters where scripts consume each other's output:
- `04_cox_forest_plot.R` reads the Excel file produced by `03_cox_nested_models.R`.
- `09_phase2_bayesian_summary.R` reads the model cached by `08_phase2_bayesian_fit.R`.

All other scripts are independent given the input data.

## Data availability

The individual-level patient data are not included in this repository. The
expected input files and their columns are documented in
[`data/README.md`](data/README.md). The dataset used and/or analysed during
the current study is available from the corresponding author on reasonable
request.

## Ethics approval

In accordance with local laws and regulations, the Institutional Ethics
Committee of The First Affiliated Hospital of Nanchang University and Jinling
Hospital, Affiliated Hospital of Medical School, Nanjing University, after
reviewing the design of the study and the anonymised research dataset,
authorised its implementation. The study was conducted in line with the
Declaration of Helsinki and adhered to the STROBE reporting guidelines.

## Funding

This work was supported by the National Natural Science Foundation of China
(No. 82370661, No. 81960128); the Double-Thousand Plan of Jiangxi Province
(No. jxsc2019201028); the Jiangxi Medicine Academy of Nutrition and Health
Management (No. 2022-PYXM-01); the Science and Technology Innovation Team
Cultivation Project of the First Affiliated Hospital of Nanchang University
(YFYKCTDPY202202); the Project for Academic and Technical Leaders of Major
Disciplines in Jiangxi Province (20243BCE51144); and the Jiangxi Provincial
Natural Science Foundation (20242BAB25438).

## Reproducibility

- Random seeds are set for all stochastic procedures (MICE imputation, GBMT,
  Stan sampling).
- The exact R version, OS, locale, and full package versions used to develop
  and validate this code are recorded in [`sessionInfo.txt`](sessionInfo.txt)
  and [`packages.tsv`](packages.tsv). Reproducing published results requires
  matching these versions closely (in particular `brms`, which sanitises
  factor level names differently across releases).
- All input paths are relative placeholders; substitute your own data files
  following the data dictionary in [`data/README.md`](data/README.md).

## Citation

If you use this code, please cite the associated publication and this
software release; the citation metadata is provided in
[`CITATION.cff`](CITATION.cff) (the author list, journal, and DOI are placed
after acceptance).

## License

Released under the MIT License — see [`LICENSE`](LICENSE).
