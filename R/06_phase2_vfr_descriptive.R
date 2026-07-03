# ============================================================
# 06 | Phase 2: fluid resuscitation (VFR) descriptive analysis
# Daily and cumulative fluid volume by subphenotype;
# Kruskal-Wallis + Dunn post-hoc; trend and boxplot figures.
#
# Input : data/phase2_lagged_primary.csv
# Output: output/06_phase2_vfr_descriptive/
# ============================================================

library(readr); library(dplyr); library(tidyr); library(ggplot2)
library(writexl); library(rstatix); library(ggpubr); library(patchwork)

out_dir <- "output/06_phase2_vfr_descriptive"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

COLORS <- c("Compensated" = "#2E7D32", "Hypodynamic" = "#1565C0",
            "Decompensated" = "#D32F2F", "Pressure-Compensated" = "#F57C00")
PHENO_LEVELS <- c("Compensated", "Hypodynamic", "Decompensated", "Pressure-Compensated")

d <- read_csv("data/phase2_lagged_primary.csv", show_col_types = FALSE) |>
  mutate(
    VFR_L   = coalesce(VFR_L_imputed, VFR_L),
    Subtype = factor(Subtype, levels = PHENO_LEVELS)
  )

## ---- 1. daily VFR by subtype ----
daily <- d |>
  group_by(Subtype, Day) |>
  summarise(mean = mean(VFR_L, na.rm = TRUE), sd = sd(VFR_L, na.rm = TRUE),
            n = n(), se = sd / sqrt(n), .groups = "drop")

## ---- 2. cumulative VFR per patient ----
cum_vfr <- d |>
  group_by(ID, Subtype) |>
  summarise(total_VFR = sum(VFR_L, na.rm = TRUE), .groups = "drop")

kw   <- kruskal.test(total_VFR ~ Subtype, data = cum_vfr)
dunn <- cum_vfr |> dunn_test(total_VFR ~ Subtype, p.adjust.method = "bonferroni")
cat("Kruskal-Wallis: H =", round(kw$statistic, 2), "p =", round(kw$p.value, 4), "\n")
print(dunn)

## ---- 3. summary tables ----
summary_tbl <- d |>
  group_by(Subtype, Day) |>
  summarise(
    N = n(),
    VFR_mean = round(mean(VFR_L, na.rm = TRUE), 2),
    VFR_sd   = round(sd(VFR_L, na.rm = TRUE), 2),
    VFR_med  = round(median(VFR_L, na.rm = TRUE), 2),
    VFR_q1   = round(quantile(VFR_L, .25, na.rm = TRUE), 2),
    VFR_q3   = round(quantile(VFR_L, .75, na.rm = TRUE), 2),
    .groups = "drop"
  )

cum_summary <- cum_vfr |>
  group_by(Subtype) |>
  summarise(N = n(), mean = round(mean(total_VFR), 2), sd = round(sd(total_VFR), 2),
            med = round(median(total_VFR), 2),
            q1 = round(quantile(total_VFR, .25), 2), q3 = round(quantile(total_VFR, .75), 2),
            .groups = "drop")

## ---- 4. figures ----
p_trend <- ggplot(daily, aes(x = Day, y = mean, color = Subtype, group = Subtype)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = Subtype), alpha = 0.12, color = NA) +
  geom_point(size = 2) +
  scale_color_manual(values = COLORS) + scale_fill_manual(values = COLORS) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5), limits = c(1, 5)) +
  labs(x = "Day", y = "Fluid intake (L/day)", title = "Daily fluid resuscitation by subphenotype") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.direction = "horizontal")

p_box <- ggplot(cum_vfr, aes(x = Subtype, y = total_VFR, fill = Subtype)) +
  geom_boxplot(outlier.size = 0.8, linewidth = 0.5, width = 0.6) +
  scale_fill_manual(values = COLORS) +
  stat_compare_means(method = "kruskal.test", label.y = max(cum_vfr$total_VFR) * 1.05, size = 3.5) +
  labs(x = NULL, y = "Cumulative fluid intake (L, Days 1-7)", title = "Total fluid intake by subphenotype") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))

p_combined <- p_trend / p_box + plot_layout(heights = c(1.2, 1))

pdf(file.path(out_dir, "Figure_Phase2_VFR_descriptive.pdf"), width = 7, height = 8)
print(p_combined); dev.off()
tiff(file.path(out_dir, "Figure_Phase2_VFR_descriptive.tiff"),
     width = 175, height = 200, units = "mm", res = 600, compression = "lzw", type = "cairo")
print(p_combined); dev.off()

write_xlsx(list(Daily_VFR = summary_tbl, Cumulative_VFR = cum_summary,
                Dunn_test = as.data.frame(dunn)),
           file.path(out_dir, "Phase2_VFR_descriptive.xlsx"))

cat("\nSaved:", out_dir, "\n")
