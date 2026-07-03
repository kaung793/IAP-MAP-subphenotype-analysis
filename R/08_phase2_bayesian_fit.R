# ============================================================
# 08 | Phase 2: Bayesian multivariate mixed-effects model (brms)
# Bivariate (dMAP, dIAP) model with correlated random effects,
# fitted with Stan. This is a long-running fit (hours). The
# fitted object is cached to models/brms_mv_model.rds and read
# back by script 09.
#
# Input : data/phase2_lagged_primary.csv
# Output: models/brms_mv_model.rds
#         output/08_phase2_bayesian/
# ============================================================

library(readr); library(dplyr); library(brms); library(tidybayes); library(writexl)

cat("brms version:", as.character(packageVersion("brms")), "\n")

dir.create("models", recursive = TRUE, showWarnings = FALSE)
out_dir <- "output/08_phase2_bayesian"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data <- read_csv("data/phase2_lagged_primary.csv", show_col_types = FALSE) %>%
  mutate(
    VFR_L   = coalesce(VFR_L_imputed, VFR_L),
    Subtype = factor(Subtype, levels = c("Compensated", "Hypodynamic",
                                         "Decompensated", "Pressure-Compensated")),
    SEX      = factor(SEX),
    Etiology = factor(Etiology),
    PCD      = as.numeric(PCD),
    NE       = as.numeric(NE),
    CRRT     = as.numeric(CRRT)
  ) %>%
  filter(!is.na(dMAP), !is.na(dIAP), !is.na(VFR_L), !is.na(Subtype))

cat("N obs:", nrow(data), "| N patients:", n_distinct(data$ID), "\n")

## ---- priors ----
priors <- c(
  prior(normal(0, 10),     class = "b",     resp = "dMAP"),
  prior(normal(0, 10),     class = "b",     resp = "dIAP"),
  prior(lkj(2),            class = "cor"),
  prior(student_t(3, 0, 5), class = "sigma", resp = "dMAP"),
  prior(student_t(3, 0, 5), class = "sigma", resp = "dIAP"),
  prior(student_t(3, 0, 5), class = "sd",    resp = "dMAP"),
  prior(student_t(3, 0, 5), class = "sd",    resp = "dIAP")
)

## ---- multivariate formula ----
formula_mv <- bf(dMAP ~ VFR_L * Subtype + PCD * Subtype +
                   NE + CRRT + MAP + IAP + Age + SEX + APACHEII1 + Etiology + (1 | p | ID)) +
              bf(dIAP ~ VFR_L * Subtype + PCD * Subtype +
                   NE + CRRT + MAP + IAP + Age + SEX + APACHEII1 + Etiology + (1 | p | ID)) +
              set_rescor(TRUE)

cat("Fitting multivariate brms model...\n")

set.seed(2026)
model_mv <- brm(
  formula_mv, data = data, family = gaussian(), prior = priors,
  chains = 4, iter = 6000, warmup = 2000, cores = 4, seed = 2026,
  file = "models/brms_mv_model"   # auto-saves .rds
)

cat("=== Convergence summary ===\n")
cat("R-hat > 1.01:", sum(rhat(model_mv) > 1.01, na.rm = TRUE), "parameters\n")
cat("ESS ratio < 0.1:", sum(neff_ratio(model_mv) < 0.1, na.rm = TRUE), "parameters\n")

## ---- marginal effects by subtype ----
post <- as_draws_df(model_mv)
subtypes <- c("Compensated", "Hypodynamic", "Decompensated", "Pressure-Compensated")

# brms sanitises factor level names with make.names() (dashes -> "M",
# spaces -> "."), and the exact rule has varied across versions. Resolve the
# interaction columns by their position: they appear in the same order as
# levels(Subtype)[-1] (i.e. Hypodynamic, Decompensated, Pressure-Compensated).
interaction_col <- function(post, treat, resp_suffix, level_index) {
  pat  <- paste0("^b_", resp_suffix, "_", treat, ":Subtype[^:]*$")
  hits <- grep(pat, names(post), value = TRUE)
  if (length(hits) == 0 || level_index > length(hits)) return(NA_character_)
  hits[level_index]
}

extract_brms_marginal <- function(post, treat, outcome_suffix) {
  ref_col <- paste0("b_", outcome_suffix, "_", treat)
  if (!ref_col %in% names(post)) stop("Reference column not found: ", ref_col)
  results <- list()
  b_ref   <- post[[ref_col]]
  results[["Compensated"]] <- data.frame(
    Subtype = "Compensated", Treatment = treat,
    Mean = mean(b_ref), SD = sd(b_ref),
    Lower95 = quantile(b_ref, 0.025), Upper95 = quantile(b_ref, 0.975),
    P_positive = mean(b_ref > 0)
  )
  for (i in seq_along(subtypes)[-1]) {
    st      <- subtypes[i]
    int_col <- interaction_col(post, treat, outcome_suffix, i - 1)
    if (is.na(int_col)) { cat("Missing interaction column for:", st, "\n"); next }
    b_tot <- b_ref + post[[int_col]]
    results[[st]] <- data.frame(
      Subtype = st, Treatment = treat,
      Mean = mean(b_tot), SD = sd(b_tot),
      Lower95 = quantile(b_tot, 0.025), Upper95 = quantile(b_tot, 0.975),
      P_positive = mean(b_tot > 0)
    )
  }
  bind_rows(results)
}

eff_VFR_dMAP <- extract_brms_marginal(post, "VFR_L", "dMAP")
eff_VFR_dIAP <- extract_brms_marginal(post, "VFR_L", "dIAP")
eff_PCD_dMAP <- extract_brms_marginal(post, "PCD",   "dMAP")
eff_PCD_dIAP <- extract_brms_marginal(post, "PCD",   "dIAP")

## ---- posterior net-benefit probability ----
net_benefit <- lapply(seq_along(subtypes), function(i) {
  st    <- subtypes[i]
  b_map <- post[["b_dMAP_VFR_L"]]
  b_iap <- post[["b_dIAP_VFR_L"]]
  if (i > 1) {
    int_map <- interaction_col(post, "VFR_L", "dMAP", i - 1)
    int_iap <- interaction_col(post, "VFR_L", "dIAP", i - 1)
    if (!is.na(int_map)) b_map <- b_map + post[[int_map]]
    if (!is.na(int_iap)) b_iap <- b_iap + post[[int_iap]]
  }
  data.frame(
    Subtype = st,
    P_MAP_improve  = mean(b_map > 0),
    P_IAP_increase = mean(b_iap > 0),
    P_net_benefit  = mean(b_map > 5 & b_iap < 3)
  )
}) |> bind_rows()

cat("\n=== Posterior net benefit probabilities ===\n"); print(net_benefit)

write_xlsx(
  list(VFR_dMAP = eff_VFR_dMAP, VFR_dIAP = eff_VFR_dIAP,
       PCD_dMAP = eff_PCD_dMAP, PCD_dIAP = eff_PCD_dIAP,
       Net_Benefit = net_benefit),
  path = file.path(out_dir, "brms_treatment_effects.xlsx")
)

cat("\nbrms analysis complete. Model: models/brms_mv_model.rds | Output:", out_dir, "\n")
