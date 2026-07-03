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
├── models/                            # cached fitted models (created at runtime)
├── output/                            # results, tables, figures (created at runtime)
├── install_packages.R                 # one-off dependency installer
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
[`data/README.md`](data/README.md). Data may be available from the corresponding
author on reasonable request, subject to institutional and ethics approvals.

## Notes on reproducibility

- Random seeds are set where results depend on stochastic procedures (MICE
  imputation, GBMT, Stan sampling).
- File paths are relative placeholders; substitute your own data files following
  the data dictionary.
