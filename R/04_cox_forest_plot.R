# ============================================================
# 04 | Cox forest plots (custom ggplot2 table + forest layout)
# Development and external validation cohorts.
#
# Input : output/03_cox_nested/cox_nested_models.xlsx  (from script 03)
# Output: output/04_cox_forest/
# ============================================================

library(readxl); library(dplyr); library(tidyr); library(ggplot2)
library(patchwork); library(cowplot)

out_dir <- "output/04_cox_forest"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

M_COLORS <- c("Model 1" = "#0D47A1", "Model 2" = "#1B5E20", "Model 3" = "#BF360C")
M_SHAPES <- c("Model 1" = 16L, "Model 2" = 15L, "Model 3" = 18L)
PHENO_ORDER <- c("Hypertensive compensated", "Hypertensive decompensated", "Hypoperfusion", "Stable")
Y_POS <- setNames(c(1, 2, 3, 4), PHENO_ORDER)

fmt_hr <- function(hr, lo, hi, p) {
  d <- function(x) if (x >= 10) 1L else 2L
  suffix <- if (p < 0.001) "*" else ""
  sprintf("%.*f (%.*f-%.*f)%s", d(hr), hr, d(lo), lo, d(hi), hi, suffix)
}

## ---- load ----
raw <- read_excel("output/03_cox_nested/cox_nested_models.xlsx",
                  sheet = "Trajectory group HRs") |>
  rowwise() |>
  mutate(
    phenotype = case_when(
      grepl("Hypod",    term) ~ "Hypoperfusion",
      grepl("Decomp",   term) ~ "Hypertensive decompensated",
      grepl("Pressure", term) ~ "Hypertensive compensated"
    ),
    model_f    = case_when(
      grepl("Model 1", model) ~ "Model 1",
      grepl("Model 2", model) ~ "Model 2",
      TRUE ~ "Model 3"
    ),
    hr_ci      = fmt_hr(HR, lower, upper, p),
    upper_clip = pmin(upper, 150)
  ) |> ungroup() |>
  mutate(
    y_base = Y_POS[phenotype],
    model_f = factor(model_f, levels = c("Model 1", "Model 2", "Model 3"))
  )

## ---- build panels ----
make_panels <- function(coh, title, xlim_max, xbreaks) {

  d <- raw |> filter(cohort == coh) |>
    mutate(y_dodge = y_base + case_when(model_f == "Model 1" ~ 0.25,
                                        model_f == "Model 2" ~ 0,
                                        TRUE                 ~ -0.25))

  txt <- d |> select(phenotype, model_f, hr_ci, y_base) |>
    pivot_wider(names_from = model_f, values_from = hr_ci)

  HDR_Y <- 4.75
  FSLAB <- 3.7
  FSHDR <- 3.9

  left <- ggplot() +
    geom_hline(yintercept = 4.5, linewidth = 0.5, color = "black") +
    geom_hline(yintercept = 0.4, linewidth = 0.5, color = "black") +
    annotate("text", x = 0.05, y = HDR_Y, label = "Phenotype", fontface = "bold", size = FSHDR, hjust = 0) +
    annotate("text", x = 0.50, y = HDR_Y, label = "Model 1", fontface = "bold", size = FSHDR, hjust = 0.5, color = M_COLORS["Model 1"]) +
    annotate("text", x = 0.70, y = HDR_Y, label = "Model 2", fontface = "bold", size = FSHDR, hjust = 0.5, color = M_COLORS["Model 2"]) +
    annotate("text", x = 0.90, y = HDR_Y, label = "Model 3", fontface = "bold", size = FSHDR, hjust = 0.5, color = M_COLORS["Model 3"]) +
    annotate("text", x = 0.05, y = 4, label = "Stable (ref)", size = FSLAB, hjust = 0, color = "gray45", fontface = "italic") +
    annotate("text", x = 0.50, y = 4, label = "Reference", size = FSLAB, hjust = 0.5, color = "gray45") +
    annotate("text", x = 0.70, y = 4, label = "Reference", size = FSLAB, hjust = 0.5, color = "gray45") +
    annotate("text", x = 0.90, y = 4, label = "Reference", size = FSLAB, hjust = 0.5, color = "gray45") +
    geom_text(data = txt, aes(x = 0.05, y = y_base, label = phenotype), size = FSLAB, hjust = 0, fontface = "bold") +
    geom_text(data = txt, aes(x = 0.50, y = y_base, label = `Model 1`), size = FSLAB - 0.3, hjust = 0.5, color = M_COLORS["Model 1"]) +
    geom_text(data = txt, aes(x = 0.70, y = y_base, label = `Model 2`), size = FSLAB - 0.3, hjust = 0.5, color = M_COLORS["Model 2"]) +
    geom_text(data = txt, aes(x = 0.90, y = y_base, label = `Model 3`), size = FSLAB - 0.3, hjust = 0.5, color = M_COLORS["Model 3"]) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.4, 5.1)) +
    labs(x = "", y = NULL) +
    theme_classic(base_size = 11, base_family = "sans") +
    theme(
      axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
      axis.title.x = element_text(size = 10.5, color = "white"), axis.title.y = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 5, r = 0, b = 5, l = 5)
    )

  right <- ggplot(d, aes(x = HR, xmin = lower, xmax = upper_clip,
                         y = y_dodge, color = model_f, shape = model_f)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4, color = "gray40") +
    geom_errorbarh(aes(xmin = lower, xmax = upper_clip), height = 0, linewidth = 0.55) +
    geom_point(size = 2.2) +
    scale_x_log10(breaks = xbreaks,
                  labels = function(x) ifelse(x == as.integer(x), as.integer(x), x),
                  limits = c(0.25, xlim_max)) +
    scale_color_manual(values = M_COLORS) +
    scale_shape_manual(values = M_SHAPES) +
    scale_y_continuous(limits = c(0.4, 5.1), breaks = 1:4, labels = NULL) +
    annotate("text", x = 1, y = 4.75, label = "HR (95% CI)", size = FSHDR, fontface = "bold", hjust = 0.5) +
    coord_cartesian(clip = "off") +
    labs(x = "Hazard ratio (log scale)", y = NULL) +
    theme_classic(base_size = 11, base_family = "sans") +
    theme(
      legend.position = "none",
      axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 10.5, margin = margin(t = 4)),
      axis.line.x = element_line(linewidth = 0.5), axis.ticks.x = element_line(linewidth = 0.4),
      panel.grid.major.x = element_line(linewidth = 0.2, color = "gray88", linetype = "dashed"),
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 0)
    ) +
    ggtitle(title)

  list(left = left, right = right)
}

panels_dev <- make_panels("Development", "Development cohort (n = 1008)", 22, c(0.5, 1, 2, 5, 10, 20))
panels_val <- make_panels("External", "External validation cohort (n = 315)", 150, c(0.5, 1, 5, 20, 100))

## ---- legend ----
leg_data <- data.frame(
  m = factor(c("Model 1", "Model 2", "Model 3"), levels = c("Model 1", "Model 2", "Model 3")),
  x = 1:3, y = 1
)
leg_plt <- ggplot(leg_data, aes(x, y, color = m, shape = m)) +
  geom_point(size = 3) +
  scale_color_manual(values = M_COLORS, labels = c("Model 1", "Model 2", "Model 3")) +
  scale_shape_manual(values = M_SHAPES, labels = c("Model 1", "Model 2", "Model 3")) +
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_blank(), legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines"), legend.spacing.x = unit(0.4, "lines"))
leg_grob <- get_legend(leg_plt)

## ---- assemble & save ----
save_fig <- function(panels, fname, w_mm = 220, h_mm = 115) {
  aligned  <- cowplot::align_plots(panels$left, panels$right, align = "hv", axis = "tblr")
  combined <- cowplot::plot_grid(aligned[[1]], aligned[[2]], ncol = 2, rel_widths = c(2.2, 1))
  g <- cowplot::plot_grid(combined, leg_grob, ncol = 1, rel_heights = c(9, 1))

  pdf(file.path(out_dir, paste0(fname, ".pdf")), width = w_mm / 25.4, height = h_mm / 25.4)
  print(g); dev.off()
  tiff(file.path(out_dir, paste0(fname, ".tiff")), width = w_mm, height = h_mm, units = "mm",
       res = 600, compression = "lzw", type = "cairo")
  print(g); dev.off()
}

save_fig(panels_dev, "Figure_Cox_Development")
save_fig(panels_val, "Figure_Cox_Validation", h_mm = 115)

cat("Saved to:", out_dir, "\n")
