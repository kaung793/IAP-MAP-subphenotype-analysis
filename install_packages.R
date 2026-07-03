# One-off installer for all CRAN packages used across the analysis scripts.
# Run once:  source("install_packages.R")

pkgs <- c(
  # data handling
  "readr", "readxl", "writexl", "dplyr", "tidyr", "tibble",
  # trajectory modelling & imputation
  "gbmt", "mice",
  # survival
  "survival", "survminer", "car",
  # mixed models
  "lme4", "lmerTest",
  # Bayesian
  "brms", "tidybayes",
  # figures & stats
  "ggplot2", "patchwork", "cowplot", "rstatix", "ggpubr"
)

to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  message("All required packages already installed.")
}

# Note: brms (script 08) needs a working C++ toolchain and a Stan backend
# (rstan or cmdstanr). See https://mc-stan.org for platform-specific setup.
