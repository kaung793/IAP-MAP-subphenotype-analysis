# ============================================================
# 05 | External validation: MICE imputation + GBMT (ng = 4)
# Re-derives the four MAP/IAP trajectory subphenotypes in the
# external cohort and cross-tabulates against original labels.
#
# Input : data/external_validation_merged.xlsx  (wide; columns
#         id, source, trajectory_group, MAP_day1..7, IAP_day1..7)
# Output: output/05_external_gbmt/
# ============================================================

library(readxl); library(dplyr); library(tidyr); library(ggplot2)
library(patchwork); library(gbmt); library(mice); library(writexl)

out_dir <- "output/05_external_gbmt"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 1. load ----
merged <- read_excel("data/external_validation_merged.xlsx")

imp_vars <- c("MAP_day1", "MAP_day2", "MAP_day3", "MAP_day5", "MAP_day7",
              "IAP_day1", "IAP_day2", "IAP_day3", "IAP_day5", "IAP_day7")

## ---- 2. MICE imputation (PMM, m = 5, first completed set) ----
set.seed(42)
imp <- mice(merged[, imp_vars], m = 5, method = "pmm", printFlag = FALSE)
complete_data <- complete(imp, action = 1)
merged_imp <- merged
merged_imp[, imp_vars] <- complete_data
cat("Missing after imputation:", sum(is.na(merged_imp[, imp_vars])), "\n")

## ---- 3. long format ----
long <- merged_imp |>
  select(id, source, trajectory_group, all_of(imp_vars)) |>
  pivot_longer(cols = all_of(imp_vars),
               names_to = c(".value", "Time"),
               names_pattern = "^(MAP|IAP)_day(\\d+)$") |>
  mutate(Time = as.integer(Time)) |>
  rename(ID = id) |>
  mutate(ID = as.character(ID))
cat("Unique IDs in long format:", n_distinct(long$ID), "\n")

## ---- 4. GBMT d = 3, ng = 4 ----
set.seed(42)
m4 <- gbmt(x.names = c("MAP", "IAP"), unit = "ID", time = "Time",
           d = 3, ng = 4, data = as.data.frame(long), scaling = 0)
cat("=== ICs ===\n");   print(m4$ic)
cat("=== APPA ===\n");  print(m4$appa)
cat("=== Prior ===\n"); print(m4$prior)

## ---- 5. phenotype labelling ----
assign_df <- data.frame(ID = names(m4$assign), new_group = as.integer(m4$assign),
                        stringsAsFactors = FALSE)

char <- long |>
  inner_join(assign_df, by = "ID") |>
  group_by(new_group) |>
  summarise(mean_MAP = mean(MAP), mean_IAP = mean(IAP),
            n = n_distinct(ID), .groups = "drop")

# Map groups to phenotypes by MAP/IAP profile (matching the development cohort)
char <- char |>
  arrange(desc(mean_IAP)) |>
  mutate(
    phenotype = case_when(
      new_group == 1 ~ "Decompensated",
      new_group == 2 ~ "Hypodynamic",
      new_group == 3 ~ "Pressure-compensated",
      new_group == 4 ~ "Compensated"
    )
  )
cat("=== Group characterization ===\n"); print(char)

compare <- assign_df |>
  left_join(select(merged, id, trajectory_group, source) |>
              mutate(id = as.character(id)),
            by = c("ID" = "id"))
cat("=== New vs original cross-tab ===\n")
print(table(New = compare$new_group, Original = compare$trajectory_group, useNA = "ifany"))

## ---- 6. summary for plot ----
plot_data <- long |>
  inner_join(assign_df, by = "ID") |>
  inner_join(select(char, new_group, phenotype, n), by = "new_group") |>
  group_by(phenotype, n, Time) |>
  summarise(
    MAP_m = mean(MAP), MAP_se = sd(MAP) / sqrt(n()),
    IAP_m = mean(IAP), IAP_se = sd(IAP) / sqrt(n()),
    .groups = "drop"
  ) |>
  mutate(
    MAP_lo = MAP_m - 1.96 * MAP_se, MAP_hi = MAP_m + 1.96 * MAP_se,
    IAP_lo = IAP_m - 1.96 * IAP_se, IAP_hi = IAP_m + 1.96 * IAP_se
  )

## ---- 7. figure ----
COLORS <- c("Compensated" = "#2E7D32", "Hypodynamic" = "#1565C0",
            "Decompensated" = "#D32F2F", "Pressure-compensated" = "#F57C00")

plot_data <- plot_data |>
  mutate(label = paste0(phenotype, " (n=", n, ")"),
         phenotype = factor(phenotype,
           levels = c("Compensated", "Hypodynamic", "Decompensated", "Pressure-compensated")))

label_map   <- plot_data |> distinct(phenotype, label) |> arrange(phenotype) |> pull(label, name = phenotype)
color_named <- setNames(COLORS[names(label_map)], label_map)

cc_theme <- function() {
  theme_classic(base_size = 9, base_family = "sans") +
    theme(
      axis.line  = element_line(linewidth = 0.5, color = "gray30"),
      axis.ticks = element_line(linewidth = 0.4, color = "gray30"),
      axis.text  = element_text(size = 8, color = "gray10"),
      axis.title = element_text(size = 9, face = "bold", color = "gray10"),
      legend.title = element_blank(), legend.text = element_text(size = 8),
      legend.position = "bottom", legend.key.width = unit(1.2, "cm"),
      aspect.ratio = 1,
      panel.grid.major.y = element_line(linewidth = 0.25, color = "gray85"),
      panel.grid.minor = element_blank(),
      plot.tag = element_text(face = "bold", size = 11, color = "gray10")
    )
}

fig_MAP <- ggplot(plot_data, aes(Time, MAP_m, color = label, fill = label)) +
  geom_ribbon(aes(ymin = MAP_lo, ymax = MAP_hi), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.1) + geom_point(size = 2.5, shape = 16) +
  geom_hline(yintercept = 65, linetype = "dotted", color = "gray50", linewidth = 0.35) +
  scale_color_manual(values = color_named) + scale_fill_manual(values = color_named) +
  scale_x_continuous(breaks = 1:7, limits = c(0.8, 7.2)) +
  labs(x = "Day", y = "MAP (mmHg)", tag = "A") +
  guides(fill = "none", color = guide_legend(nrow = 2)) + cc_theme()

fig_IAP <- ggplot(plot_data, aes(Time, IAP_m, color = label, fill = label)) +
  geom_ribbon(aes(ymin = IAP_lo, ymax = IAP_hi), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.1) + geom_point(size = 2.5, shape = 16) +
  geom_hline(yintercept = 12, linetype = "dotted", color = "gray50", linewidth = 0.35) +
  scale_color_manual(values = color_named) + scale_fill_manual(values = color_named) +
  scale_x_continuous(breaks = 1:7, limits = c(0.8, 7.2)) +
  labs(x = "Day", y = "IAP (mmHg)", tag = "B") +
  guides(fill = "none", color = guide_legend(nrow = 2)) + cc_theme()

fig <- fig_MAP | fig_IAP + plot_layout(guides = "collect") & theme(legend.position = "bottom")
fig <- fig + plot_annotation(
  caption = "Shaded areas: 95% CI. Dotted lines: MAP 65 mmHg (A) and IAP 12 mmHg (B). MICE-imputed (m=1).",
  theme = theme(plot.caption = element_text(size = 7.5, color = "gray40", hjust = 0))
)

ggsave(file.path(out_dir, "Figure_external_validation_gbmt.pdf"), fig,
       width = 180 / 25.4, height = 125 / 25.4, device = cairo_pdf)
ggsave(file.path(out_dir, "Figure_external_validation_gbmt.tiff"), fig,
       width = 180 / 25.4, height = 125 / 25.4, dpi = 600, compression = "lzw", device = "tiff")

write_xlsx(
  compare |> left_join(select(char, new_group, phenotype), by = "new_group"),
  file.path(out_dir, "external_gbmt_assignments.xlsx")
)

cat("Saved to:", out_dir, "\n")
