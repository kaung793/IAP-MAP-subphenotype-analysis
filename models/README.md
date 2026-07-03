# Fitted models

This directory holds cached model objects created at runtime.

- `brms_mv_model.rds` — the Bayesian multivariate mixed-effects model, written by
  `R/08_phase2_bayesian_fit.R` and read back by `R/09_phase2_bayesian_summary.R`.

The `.rds` files are large and are not committed to version control (see
`.gitignore`). Run `R/08_phase2_bayesian_fit.R` to (re)generate them.
