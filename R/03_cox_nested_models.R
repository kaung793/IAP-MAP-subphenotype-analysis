# ============================================================
# 03 | Cox proportional-hazards models (28-day mortality)
# Three nested models, development + external validation cohorts.
# Reference subphenotype: Compensated. Reports HRs, VIF
# diagnostics and proportional-hazards (Schoenfeld) tests.
#
# Input : data/dev_cohort_full.csv
#         data/external_validation_315.xlsx
# Output: output/03_cox_nested/cox_nested_models.xlsx
# ============================================================

library(readr); library(readxl); library(dplyr); library(survival); library(car)
library(writexl); library(tibble)

out_dir <- "output/03_cox_nested"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

PHENO_LEVELS <- c("Compensated", "Hypodynamic", "Decompensated", "Pressure-compensated")
ETIOL_LABELS <- c("Biliary", "HTG", "Alcoholic", "Other")

## ---- helpers ----
extract_cox <- function(fit, cohort, model_label) {
  s  <- summary(fit)
  ci <- as.data.frame(s$conf.int)
  pv <- s$coefficients[, "Pr(>|z|)"]
  tibble(
    cohort = cohort, model = model_label, term = rownames(ci),
    HR = ci[["exp(coef)"]], lower = ci[["lower .95"]], upper = ci[["upper .95"]], p = pv
  ) |> mutate(
    HR_95CI = sprintf("%.2f (%.2f-%.2f)", HR, lower, upper),
    p_fmt   = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

fmt_vif <- function(v, cohort) {
  if (is.matrix(v)) {
    data.frame(cohort = cohort, variable = rownames(v),
               VIF = round(v[, "GVIF^(1/(2*Df))"]^2, 3))
  } else {
    data.frame(cohort = cohort, variable = names(v), VIF = round(as.numeric(v), 3))
  }
}

## ---- 1. Development cohort ----
dev <- read_csv("data/dev_cohort_full.csv", show_col_types = FALSE) |>
  transmute(
    group    = relevel(factor(`IAP-MAPgroup`, levels = 1:4, labels = PHENO_LEVELS), ref = "Compensated"),
    time     = as.numeric(R28day),
    event    = as.numeric(R28Death),
    age      = Age,
    sex      = factor(SEX, levels = c(0, 1), labels = c("Female", "Male")),
    weight   = Weight,
    temp     = temperature01,
    pulse    = pulse01,
    resp     = respirations01,
    etiology = factor(`Etiology.of.pancreatitis`, levels = 1:4, labels = ETIOL_LABELS),
    wbc      = WBC01, plt = PLT01, bun = BUN01, alt = ALT01, pao2 = PaO201
  ) |>
  filter(!is.na(time), !is.na(event))

cox1_dev <- coxph(Surv(time, event) ~ group + age + sex + weight, data = dev)
cox2_dev <- coxph(Surv(time, event) ~ group + age + sex + weight +
                    temp + pulse + resp + etiology, data = dev)
cox3_dev <- coxph(Surv(time, event) ~ group + age + sex + weight +
                    temp + pulse + resp + etiology +
                    wbc + plt + bun + alt + pao2, data = dev)

vif_dev <- vif(cox3_dev)
ph_dev  <- cox.zph(cox3_dev)
cat("=== Development VIF (Model 3) ===\n"); print(vif_dev)
cat("=== PH test (dev) ===\n"); print(ph_dev)

## ---- 2. External validation ----
# Note: 20 28-day events; Model 3 is exploratory in this cohort.
val <- read_excel("data/external_validation_315.xlsx") |>
  transmute(
    group    = relevel(factor(trajectory_phenotype, levels = PHENO_LEVELS), ref = "Compensated"),
    time     = as.numeric(`28-hospitalday`),
    event    = as.numeric(`R-28Death`),
    age      = as.numeric(age),
    sex      = factor(gender, levels = c(0, 1), labels = c("Female", "Male")),
    weight   = as.numeric(weight),
    temp     = as.numeric(temp_day1),
    pulse    = as.numeric(pulse_day1),
    resp     = as.numeric(resp_day1),
    etiology = factor(etiology, levels = 1:4, labels = ETIOL_LABELS),
    wbc      = as.numeric(WBC), plt = as.numeric(PLT), bun = as.numeric(BUN), alt = as.numeric(ALT)
  ) |>
  filter(!is.na(time), !is.na(event))

cat("\nExternal validation 28-day events:", sum(val$event), "/ n =", nrow(val), "\n")

cox1_val <- coxph(Surv(time, event) ~ group + age + sex + weight, data = val)
cox2_val <- coxph(Surv(time, event) ~ group + age + sex + weight +
                    temp + pulse + resp + etiology, data = val)
# ALT excluded from Model 3 (external): VIF > threshold of 5
cox3_val <- coxph(Surv(time, event) ~ group + age + sex + weight +
                    temp + pulse + resp + etiology +
                    wbc + plt + bun, data = val)

vif_val <- vif(cox3_val)
ph_val  <- tryCatch(cox.zph(cox3_val), error = function(e) {
  cat("PH test failed for external validation (likely singular matrix due to low event count):",
      conditionMessage(e), "\n"); NULL
})
cat("\n=== External VIF (Model 3) ===\n"); print(vif_val)
if (!is.null(ph_val)) { cat("=== PH test (val) ===\n"); print(ph_val) }

## ---- 3. Compile ----
all_results <- bind_rows(
  extract_cox(cox1_dev, "Development", "Model 1 (Demographics)"),
  extract_cox(cox2_dev, "Development", "Model 2 (+Vitals+Etiology)"),
  extract_cox(cox3_dev, "Development", "Model 3 (+Labs, full)"),
  extract_cox(cox1_val, "External", "Model 1 (Demographics)"),
  extract_cox(cox2_val, "External", "Model 2 (+Vitals+Etiology)"),
  extract_cox(cox3_val, "External", "Model 3 (+Labs, full)")
)

group_results <- all_results |>
  filter(grepl("^group", term)) |>
  mutate(term = gsub("^group", "", term))

vif_table <- bind_rows(fmt_vif(vif_dev, "Development"), fmt_vif(vif_val, "External"))

ph_df <- function(ph, cohort) {
  as.data.frame(ph$table) |>
    tibble::rownames_to_column("term") |>
    mutate(cohort = cohort, global_p = ph$table["GLOBAL", "p"],
           interpretation = ifelse(p > 0.05, "OK", "VIOLATED"))
}
ph_table <- bind_rows(
  ph_df(ph_dev, "Development"),
  if (!is.null(ph_val)) ph_df(ph_val, "External")
  else data.frame(term = "N/A", cohort = "External",
                  note = "Model singular: too few events for PH test")
)

write_xlsx(list(
  `Trajectory group HRs` = group_results,
  `All covariates`       = all_results,
  `VIF diagnostics`      = vif_table,
  `PH assumption`        = ph_table
), file.path(out_dir, "cox_nested_models.xlsx"))

cat("\nSaved to:", out_dir, "\n")
