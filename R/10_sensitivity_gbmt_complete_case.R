# ============================================================
# 10 | Sensitivity: complete-case GBMT
# Re-fits the ng = 4 trajectory model on patients with all five
# MAP + IAP time points observed, and checks agreement with the
# main (imputed) subphenotype assignment.
#
# Input : data/dev_cohort_full.csv
# Output: output/10_sensitivity_gbmt_complete_case/
# ============================================================

library(readr); library(dplyr); library(tidyr); library(gbmt); library(writexl); library(ggplot2)

out_dir <- "output/10_sensitivity_gbmt_complete_case"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

COLORS <- c("Compensated" = "#2E7D32", "Hypodynamic" = "#1565C0",
            "Decompensated" = "#D32F2F", "Pressure-Compensated" = "#F57C00")

d <- read_csv("data/dev_cohort_full.csv", show_col_types = FALSE)

iap_cols <- c("IAP 01", "IAP 02", "IAP 03", "IAP 05", "IAP 07")
map_cols <- c("MeanArterialPressure01", "meanarterialpressure02",
              "MeanArterialPressure03", "MeanArterialPressure05", "MeanArterialPressure07")

cc <- d |> filter(if_all(all_of(c(iap_cols, map_cols)), ~ !is.na(.x)))
cat("Complete cases:", nrow(cc), "/", nrow(d), sprintf("(%.1f%%)\n", 100 * nrow(cc) / nrow(d)))

## ---- complete vs incomplete baseline ----
missing_flag <- d |>
  mutate(complete = if_all(all_of(c(iap_cols, map_cols)), ~ !is.na(.x))) |>
  group_by(complete) |>
  summarise(n = n(), death_pct = round(mean(R28Death, na.rm = TRUE) * 100, 1),
            age = round(mean(Age, na.rm = TRUE), 1),
            apache = round(mean(APACHEII1, na.rm = TRUE), 1), .groups = "drop")
cat("\nComplete vs incomplete baseline:\n"); print(missing_flag)

## ---- reshape to long ----
long <- cc |>
  select(ID = `ID.study.group.`, `IAP-MAPgroup`, all_of(iap_cols), all_of(map_cols)) |>
  pivot_longer(cols = all_of(iap_cols), names_to = "iap_var", values_to = "IAP") |>
  mutate(Time = c("IAP 01" = 1, "IAP 02" = 2, "IAP 03" = 3, "IAP 05" = 5, "IAP 07" = 7)[iap_var]) |>
  left_join(
    cc |> select(ID = `ID.study.group.`, all_of(map_cols)) |>
      pivot_longer(all_of(map_cols), names_to = "map_var", values_to = "MAP") |>
      mutate(Time = c("MeanArterialPressure01" = 1, "meanarterialpressure02" = 2,
                      "MeanArterialPressure03" = 3, "MeanArterialPressure05" = 5,
                      "MeanArterialPressure07" = 7)[map_var]) |>
      select(ID, Time, MAP),
    by = c("ID", "Time")
  ) |>
  select(ID, Time, IAP, MAP) |>
  mutate(ID = as.character(ID)) |>
  filter(!is.na(IAP), !is.na(MAP))

cat("\nLong format:", nrow(long), "obs,", n_distinct(long$ID), "patients\n")

## ---- GBMT ng = 4, d = 3 ----
set.seed(2026)
fit_cc <- gbmt(x.names = c("MAP", "IAP"), unit = "ID", time = "Time",
               d = 3, ng = 4, scaling = 0, data = as.data.frame(long))

cat("\n=== Complete-case GBMT (ng = 4) ===\n")
cat("BIC:", fit_cc$ic["bic"], "\n")

assign_cc <- data.frame(ID = as.character(names(fit_cc$assign)), group_cc = as.integer(fit_cc$assign))

appa <- fit_cc$appa
cat("APPA by group:", round(appa, 3), "\n")
cat("Group sizes:", table(assign_cc$group_cc), "\n")

## ---- auto-label (same logic as main analysis) ----
group_chars <- long |>
  left_join(assign_cc, by = "ID") |>
  group_by(group_cc) |>
  summarise(mean_IAP = mean(IAP, na.rm = TRUE), mean_MAP = mean(MAP, na.rm = TRUE)) |>
  arrange(desc(mean_IAP)) |>
  mutate(iap_tier = if_else(row_number() <= 2, "high", "low")) |>
  group_by(iap_tier) |>
  mutate(phenotype = case_when(
    iap_tier == "high" & mean_MAP == max(mean_MAP) ~ "Pressure-Compensated",
    iap_tier == "high"                             ~ "Decompensated",
    iap_tier == "low"  & mean_MAP == max(mean_MAP) ~ "Compensated",
    TRUE                                           ~ "Hypodynamic"
  )) |> ungroup()

cat("\nComplete-case phenotype mapping:\n")
print(group_chars[, c("group_cc", "phenotype", "mean_MAP", "mean_IAP")])

## ---- agreement with main analysis ----
orig <- cc |>
  transmute(ID = as.character(`ID.study.group.`), group_orig = `IAP-MAPgroup`) |>
  left_join(assign_cc, by = "ID") |>
  left_join(group_chars |> select(group_cc, phenotype_cc = phenotype), by = "group_cc")

agree_tbl <- table(orig = orig$group_orig, cc = orig$group_cc)
cat("\nCross-tabulation (original group vs complete-case group):\n"); print(agree_tbl)
pct_agree <- sum(diag(agree_tbl)) / sum(agree_tbl)
cat(sprintf("Agreement: %.1f%%\n", pct_agree * 100))

## ---- save ----
write_xlsx(list(
  APPA        = data.frame(group = 1:4, APPA = round(appa, 3), n = as.integer(table(assign_cc$group_cc))),
  Assignments = assign_cc |> left_join(group_chars |> select(group_cc, phenotype), by = "group_cc"),
  CrossTab    = as.data.frame.matrix(agree_tbl),
  Baseline_compare = missing_flag,
  Group_chars = as.data.frame(group_chars)
), file.path(out_dir, "Sensitivity_GBMT_complete_case.xlsx"))

cat("\nSaved:", out_dir, "\n")
