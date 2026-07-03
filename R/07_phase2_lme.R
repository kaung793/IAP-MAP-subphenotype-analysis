# ============================================================
# 07 | Phase 2: linear mixed-effects models (LME)
# VFR x Subphenotype interaction on next-day change in MAP/IAP
# (dMAP, dIAP). Likelihood-ratio test for the interaction,
# subtype-specific marginal effects, and a MAP-IAP phase-space
# vector figure.
#
# Input : data/phase2_lagged_complete_case.csv
# Output: output/07_phase2_lme/
# ============================================================

library(readr); library(dplyr); library(tidyr); library(ggplot2)
library(lme4); library(lmerTest); library(writexl)

out_dir <- "output/07_phase2_lme"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

COLORS <- c("Compensated" = "#2E7D32", "Hypodynamic" = "#1565C0",
            "Decompensated" = "#D32F2F", "Pressure-Compensated" = "#F57C00")
PHENO_LEVELS <- c("Compensated", "Hypodynamic", "Decompensated", "Pressure-Compensated")

## ---- load & fit ----
d <- read_csv("data/phase2_lagged_complete_case.csv", show_col_types = FALSE) |>
  mutate(Subtype = factor(Subtype, levels = PHENO_LEVELS)) |>
  filter(!is.na(dMAP), !is.na(dIAP), !is.na(VFR_L), !is.na(Subtype),
         !is.na(NE), !is.na(CRRT), !is.na(MAP), !is.na(IAP),
         !is.na(Age), !is.na(SEX), !is.na(APACHEII1), !is.na(Etiology))

cat("n_obs =", nrow(d), "| n_patients =", n_distinct(d$ID), "\n\n")
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

f_dMAP  <- dMAP ~ VFR_L * Subtype + NE + CRRT + MAP + IAP + Age + SEX + APACHEII1 + Etiology + (1 | ID)
f_dIAP  <- dIAP ~ VFR_L * Subtype + NE + CRRT + MAP + IAP + Age + SEX + APACHEII1 + Etiology + (1 | ID)
f0_dMAP <- dMAP ~ VFR_L + Subtype + NE + CRRT + MAP + IAP + Age + SEX + APACHEII1 + Etiology + (1 | ID)
f0_dIAP <- dIAP ~ VFR_L + Subtype + NE + CRRT + MAP + IAP + Age + SEX + APACHEII1 + Etiology + (1 | ID)

m_dMAP      <- lmer(f_dMAP,  data = d, REML = FALSE, control = ctrl)
m_dIAP      <- lmer(f_dIAP,  data = d, REML = FALSE, control = ctrl)
m_dMAP_null <- lmer(f0_dMAP, data = d, REML = FALSE, control = ctrl)
m_dIAP_null <- lmer(f0_dIAP, data = d, REML = FALSE, control = ctrl)
cat("Models fitted.\n\n")

## ---- LRT: interaction significance ----
lrt_MAP <- anova(m_dMAP_null, m_dMAP)
lrt_IAP <- anova(m_dIAP_null, m_dIAP)
cat("LRT dMAP  Chi2 =", round(lrt_MAP$Chisq[2], 2), "df =", lrt_MAP$Df[2],
    "p =", round(lrt_MAP$`Pr(>Chisq)`[2], 4), "\n")
cat("LRT dIAP  Chi2 =", round(lrt_IAP$Chisq[2], 2), "df =", lrt_IAP$Df[2],
    "p =", round(lrt_IAP$`Pr(>Chisq)`[2], 4), "\n\n")

## ---- marginal effects (delta method) ----
marginal <- function(model, treat, subtypes) {
  fe <- fixef(model); vc <- vcov(model)
  ref_b <- fe[treat]; ref_v <- vc[treat, treat]
  res <- list(data.frame(Subtype = subtypes[1], Treatment = treat,
    Effect = ref_b, SE = sqrt(ref_v),
    Lower95 = ref_b - 1.96 * sqrt(ref_v), Upper95 = ref_b + 1.96 * sqrt(ref_v)))
  for (st in subtypes[-1]) {
    int <- paste0(treat, ":Subtype", st)
    if (!int %in% names(fe)) int <- paste0("Subtype", st, ":", treat)
    if (!int %in% names(fe)) next
    b <- ref_b + fe[int]
    v <- ref_v + vc[int, int] + 2 * vc[treat, int]
    res[[st]] <- data.frame(Subtype = st, Treatment = treat, Effect = b, SE = sqrt(v),
      Lower95 = b - 1.96 * sqrt(v), Upper95 = b + 1.96 * sqrt(v))
  }
  bind_rows(res) |> mutate(p_value = 2 * pnorm(-abs(Effect / SE)),
    p_fmt = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)),
    effect_ci = sprintf("%.2f (%.2f, %.2f)", Effect, Lower95, Upper95))
}

subtypes <- PHENO_LEVELS
eff_VFR_dMAP <- marginal(m_dMAP, "VFR_L", subtypes)
eff_VFR_dIAP <- marginal(m_dIAP, "VFR_L", subtypes)

cat("=== VFR +1L -> dMAP ===\n"); print(eff_VFR_dMAP[, c("Subtype", "effect_ci", "p_fmt")])
cat("=== VFR +1L -> dIAP ===\n"); print(eff_VFR_dIAP[, c("Subtype", "effect_ci", "p_fmt")])

## ---- phase-space vectors ----
starts <- d |> filter(Day == 3) |>
  group_by(Subtype) |>
  summarise(MAP0 = median(MAP, na.rm = TRUE), IAP0 = median(IAP, na.rm = TRUE), .groups = "drop")

SCALE <- 3
plot_phase <- function(pdata, title) {
  ggplot(pdata, aes(color = Subtype)) +
    geom_abline(slope = 1, intercept = 70, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    annotate("text", x = 95, y = 26, label = "APP = 70 mmHg", size = 3, color = "gray50", angle = 45) +
    geom_point(aes(x = MAP0, y = IAP0), size = 3, shape = 16) +
    geom_segment(aes(x = MAP0, y = IAP0, xend = MAP1, yend = IAP1),
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed"), linewidth = 0.9) +
    geom_text(aes(x = MAP0 - 1, y = IAP0 + 0.8, label = Subtype), size = 3, hjust = 1,
              fontface = "italic", check_overlap = TRUE) +
    scale_color_manual(values = COLORS) +
    scale_x_continuous(limits = c(60, 110), breaks = seq(60, 110, 10)) +
    scale_y_continuous(limits = c(6, 22), breaks = seq(6, 22, 4)) +
    labs(x = "MAP (mmHg)", y = "IAP (mmHg)", title = title,
         caption = paste0("Arrow = effect of +1L fluid (x", SCALE, " scale); start = median Day-3")) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 12),
          plot.caption = element_text(size = 8, color = "gray50"))
}

phase_VFR <- eff_VFR_dMAP |>
  rename(dMAP = Effect, dMAP_lo = Lower95, dMAP_hi = Upper95) |>
  left_join(eff_VFR_dIAP |> rename(dIAP = Effect, dIAP_lo = Lower95, dIAP_hi = Upper95),
            by = c("Subtype", "Treatment")) |>
  left_join(starts, by = "Subtype") |>
  mutate(MAP1 = MAP0 + dMAP * SCALE, IAP1 = IAP0 + dIAP * SCALE)

p_VFR <- plot_phase(phase_VFR, "Effect of +1 L fluid resuscitation on MAP-IAP trajectory")

pdf(file.path(out_dir, "Figure_Phase2_PhaseSpace.pdf"), width = 6, height = 5)
print(p_VFR); dev.off()
tiff(file.path(out_dir, "Figure_Phase2_PhaseSpace.tiff"),
     width = 150, height = 125, units = "mm", res = 600, compression = "lzw", type = "cairo")
print(p_VFR); dev.off()

write_xlsx(list(
  VFR_dMAP  = eff_VFR_dMAP,
  VFR_dIAP  = eff_VFR_dIAP,
  Phase_VFR = phase_VFR,
  LRT_DMAP  = as.data.frame(lrt_MAP),
  LRT_DIAP  = as.data.frame(lrt_IAP)
), file.path(out_dir, "Phase2_LME_results.xlsx"))

cat("\nSaved:", out_dir, "\n")
