# ============================================================
# 11 | Sensitivity: complete-case vs imputed LME (Phase 2)
# Compares the VFR x Subphenotype interaction and marginal
# effects between the imputed (primary) fluid dataset and the
# observed-only (complete-case) dataset.
#
# Input : data/phase2_lagged_primary.csv
#         data/phase2_lagged_complete_case.csv
# Output: output/11_sensitivity_lme_complete_case/
# ============================================================

library(readr); library(dplyr); library(lme4); library(lmerTest); library(writexl)

out_dir <- "output/11_sensitivity_lme_complete_case"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

PHENO_LEVELS <- c("Compensated", "Hypodynamic", "Decompensated", "Pressure-Compensated")
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

fmt <- function(b, lo, hi, p)
  sprintf("%.2f (%.2f, %.2f)%s", b, lo, hi, ifelse(p < 0.001, "*", ""))

marginal <- function(model, treat, subtypes) {
  fe <- fixef(model); vc <- vcov(model)
  ref_b <- fe[treat]; ref_v <- vc[treat, treat]
  res <- list(data.frame(Subtype = subtypes[1], Effect = ref_b, SE = sqrt(ref_v),
    Lower95 = ref_b - 1.96 * sqrt(ref_v), Upper95 = ref_b + 1.96 * sqrt(ref_v)))
  for (st in subtypes[-1]) {
    int <- paste0(treat, ":Subtype", st)
    if (!int %in% names(fe)) int <- paste0("Subtype", st, ":", treat)
    if (!int %in% names(fe)) next
    b <- ref_b + fe[int]; v <- ref_v + vc[int, int] + 2 * vc[treat, int]
    res[[st]] <- data.frame(Subtype = st, Effect = b, SE = sqrt(v),
      Lower95 = b - 1.96 * sqrt(v), Upper95 = b + 1.96 * sqrt(v))
  }
  bind_rows(res) |> mutate(p = 2 * pnorm(-abs(Effect / SE)),
    effect_ci = fmt(Effect, Lower95, Upper95, p))
}

run_model <- function(df, label) {
  cat("\n===", label, "===\n")
  cat("n_obs:", nrow(df), "| n_patients:", n_distinct(df$ID), "\n")
  f  <- dMAP ~ VFR_L * Subtype + PCD * Subtype + NE + CRRT + MAP + IAP + Age + SEX + APACHEII1 + Etiology + (1 | ID)
  f0 <- dMAP ~ VFR_L + PCD + Subtype + NE + CRRT + MAP + IAP + Age + SEX + APACHEII1 + Etiology + (1 | ID)
  g  <- dIAP ~ VFR_L * Subtype + PCD * Subtype + NE + CRRT + MAP + IAP + Age + SEX + APACHEII1 + Etiology + (1 | ID)
  g0 <- dIAP ~ VFR_L + PCD + Subtype + NE + CRRT + MAP + IAP + Age + SEX + APACHEII1 + Etiology + (1 | ID)

  mM  <- lmer(f,  df, REML = FALSE, control = ctrl)
  mI  <- lmer(g,  df, REML = FALSE, control = ctrl)
  m0M <- lmer(f0, df, REML = FALSE, control = ctrl)
  m0I <- lmer(g0, df, REML = FALSE, control = ctrl)

  lrt_M <- anova(m0M, mM); lrt_I <- anova(m0I, mI)
  cat(sprintf("LRT dMAP: Chi2=%.2f df=%d p=%.4f\n", lrt_M$Chisq[2], lrt_M$Df[2], lrt_M$`Pr(>Chisq)`[2]))
  cat(sprintf("LRT dIAP: Chi2=%.2f df=%d p=%.4f\n", lrt_I$Chisq[2], lrt_I$Df[2], lrt_I$`Pr(>Chisq)`[2]))

  eff_M <- marginal(mM, "VFR_L", PHENO_LEVELS) |> mutate(Outcome = "dMAP", Dataset = label)
  eff_I <- marginal(mI, "VFR_L", PHENO_LEVELS) |> mutate(Outcome = "dIAP", Dataset = label)
  lrt_df <- data.frame(Dataset = label,
    LRT_dMAP_chi2 = round(lrt_M$Chisq[2], 2), LRT_dMAP_p = round(lrt_M$`Pr(>Chisq)`[2], 4),
    LRT_dIAP_chi2 = round(lrt_I$Chisq[2], 2), LRT_dIAP_p = round(lrt_I$`Pr(>Chisq)`[2], 4))
  list(eff = bind_rows(eff_M, eff_I), lrt = lrt_df)
}

base_vars <- c("dMAP", "dIAP", "VFR_L", "Subtype", "PCD", "NE", "CRRT",
               "MAP", "IAP", "Age", "SEX", "APACHEII1", "Etiology")

## ---- imputed (primary) ----
d_imp <- read_csv("data/phase2_lagged_primary.csv", show_col_types = FALSE) |>
  mutate(VFR_L = coalesce(VFR_L_imputed, VFR_L),
         Subtype = factor(Subtype, levels = PHENO_LEVELS),
         SEX = factor(SEX), Etiology = factor(Etiology),
         PCD = as.numeric(PCD), NE = as.numeric(NE), CRRT = as.numeric(CRRT)) |>
  filter(complete.cases(across(all_of(base_vars))))

## ---- complete-case (observed FRV only) ----
d_cc <- read_csv("data/phase2_lagged_complete_case.csv", show_col_types = FALSE) |>
  mutate(VFR_L = VFR_L,
         Subtype = factor(Subtype, levels = PHENO_LEVELS),
         SEX = factor(SEX), Etiology = factor(Etiology),
         PCD = as.numeric(PCD), NE = as.numeric(NE), CRRT = as.numeric(CRRT)) |>
  filter(complete.cases(across(all_of(base_vars))))

res_imp <- run_model(d_imp, "Imputed (primary)")
res_cc  <- run_model(d_cc,  "Complete-case")

## ---- comparison ----
compare_eff <- bind_rows(res_imp$eff, res_cc$eff) |>
  select(Dataset, Outcome, Subtype, effect_ci, p) |>
  arrange(Outcome, Subtype, Dataset)
compare_lrt <- bind_rows(res_imp$lrt, res_cc$lrt)

cat("\n=== LRT comparison ===\n"); print(compare_lrt)
cat("\n=== Marginal effects comparison (VFR +1L) ===\n"); print(compare_eff |> select(-p))

write_xlsx(list(
  LRT_comparison       = compare_lrt,
  Marginal_effects     = compare_eff,
  Imputed_details      = res_imp$eff,
  CompleteCase_details = res_cc$eff
), file.path(out_dir, "Sensitivity_LME_complete_case.xlsx"))

cat("\nSaved:", out_dir, "\n")
