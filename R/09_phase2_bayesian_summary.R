# ============================================================
# 09 | Phase 2: Bayesian model summary & diagnostics
# Reads the cached brms fit from script 08, reports convergence
# (R-hat, ESS), subtype-specific posterior marginal effects,
# and posterior predictive checks.
#
# Input : models/brms_mv_model.rds  (produced by script 08)
# Output: output/09_phase2_bayesian_summary/
# ============================================================

library(brms); library(dplyr); library(tidybayes); library(ggplot2); library(writexl); library(patchwork)

out_dir <- "output/09_phase2_bayesian_summary"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

PHENO_LEVELS <- c("Compensated", "Hypodynamic", "Decompensated", "Pressure-Compensated")

cat("Loading brms model...\n")
model_mv <- readRDS("models/brms_mv_model.rds")

## ---- convergence diagnostics ----
rhat_max <- max(rhat(model_mv), na.rm = TRUE)
ess_min  <- min(neff_ratio(model_mv), na.rm = TRUE)
n_rhat   <- sum(rhat(model_mv) > 1.01, na.rm = TRUE)
n_ess    <- sum(neff_ratio(model_mv) < 0.1, na.rm = TRUE)
cat(sprintf("R-hat max = %.3f | parameters > 1.01: %d\n", rhat_max, n_rhat))
cat(sprintf("ESS ratio min = %.3f | parameters < 0.1: %d\n", ess_min, n_ess))

## ---- posterior marginal effects (VFR x Subtype) ----
post <- as_draws_df(model_mv)

# brms sanitises factor level names with make.names() (spaces -> ".", "-" -> "M"
# in current versions), but the rule has varied. Resolve interaction columns by
# position: they appear in the same order as levels(Subtype)[-1].
interaction_col <- function(post, treat, resp_suffix, level_index) {
  pat <- paste0("^b_", resp_suffix, "_", treat, ":Subtype[^:]*$")
  hits <- grep(pat, names(post), value = TRUE)
  if (length(hits) == 0) return(NA_character_)
  if (level_index > length(hits)) return(NA_character_)
  hits[level_index]
}

extract_marginal <- function(post, treat, resp_suffix, subtypes) {
  ref_col <- paste0("b_", resp_suffix, "_", treat)
  if (!ref_col %in% names(post)) stop("Reference column not found: ", ref_col)
  b_ref <- post[[ref_col]]
  res <- list(data.frame(
    Subtype = subtypes[1], Treatment = treat, Outcome = resp_suffix,
    Mean = mean(b_ref), SD = sd(b_ref),
    Lower95 = quantile(b_ref, .025), Upper95 = quantile(b_ref, .975),
    P_gt0 = mean(b_ref > 0), P_lt0 = mean(b_ref < 0)
  ))
  for (i in seq_along(subtypes)[-1]) {
    st      <- subtypes[i]
    int_col <- interaction_col(post, treat, resp_suffix, i - 1)
    if (is.na(int_col)) { cat("WARNING: column not found for", st, "\n"); next }
    b_tot <- b_ref + post[[int_col]]
    res[[st]] <- data.frame(
      Subtype = st, Treatment = treat, Outcome = resp_suffix,
      Mean = mean(b_tot), SD = sd(b_tot),
      Lower95 = quantile(b_tot, .025), Upper95 = quantile(b_tot, .975),
      P_gt0 = mean(b_tot > 0), P_lt0 = mean(b_tot < 0)
    )
  }
  bind_rows(res) |> mutate(effect_ci = sprintf("%.2f (%.2f, %.2f)", Mean, Lower95, Upper95))
}

bayes_VFR_dMAP <- extract_marginal(post, "VFR_L", "dMAP", PHENO_LEVELS)
bayes_VFR_dIAP <- extract_marginal(post, "VFR_L", "dIAP", PHENO_LEVELS)

cat("\n=== Bayesian VFR +1L -> dMAP ===\n"); print(bayes_VFR_dMAP[, c("Subtype", "effect_ci", "P_gt0", "P_lt0")])
cat("\n=== Bayesian VFR +1L -> dIAP ===\n"); print(bayes_VFR_dIAP[, c("Subtype", "effect_ci", "P_gt0", "P_lt0")])

## ---- comparison table ----
comparison <- bayes_VFR_dMAP |>
  transmute(Subtype, Bayes_dMAP = effect_ci, P_gt0_dMAP = round(P_gt0, 3)) |>
  left_join(
    bayes_VFR_dIAP |> transmute(Subtype, Bayes_dIAP = effect_ci, P_gt0_dIAP = round(P_gt0, 3)),
    by = "Subtype")
cat("\n=== Posterior effect summary ===\n"); print(comparison)

## ---- posterior predictive checks ----
pdf(file.path(out_dir, "Figure_Phase2_brms_diagnostics.pdf"), width = 10, height = 5)
p1 <- pp_check(model_mv, resp = "dMAP", ndraws = 50) +
  labs(title = "Posterior predictive check: dMAP") + theme_classic(base_size = 10)
p2 <- pp_check(model_mv, resp = "dIAP", ndraws = 50) +
  labs(title = "Posterior predictive check: dIAP") + theme_classic(base_size = 10)
print(p1 + p2)
dev.off()

## ---- save ----
write_xlsx(list(
  Bayes_VFR_dMAP = bayes_VFR_dMAP,
  Bayes_VFR_dIAP = bayes_VFR_dIAP,
  Comparison     = comparison,
  Convergence    = data.frame(rhat_max, ess_min, n_rhat_gt1.01 = n_rhat, n_ess_lt0.1 = n_ess)
), file.path(out_dir, "Phase2_Bayesian_results.xlsx"))

cat("\nSaved:", out_dir, "\n")
