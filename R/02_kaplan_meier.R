# ============================================================
# 02 | Kaplan-Meier survival by trajectory subphenotype
# 28-day mortality, development cohort, log-rank test.
#
# Input : data/dev_cohort_wide.xlsx   (must contain columns
#         'IAP-MAPgroup' (1-4), 'R28day', 'R28Death')
# Output: output/02_kaplan_meier/
# ============================================================

options(stringsAsFactors = FALSE)

## ---- packages ----
pkgs <- c("readxl", "writexl", "dplyr", "survival", "survminer", "ggplot2")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

## ---- paths ----
input_file <- "data/dev_cohort_wide.xlsx"
output_dir <- "output/02_kaplan_meier"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ---- read data ----
df <- as.data.frame(readxl::read_excel(input_file, .name_repair = "minimal"))

time_col  <- "R28day"
event_candidates <- c("R28Death", "R28Death...456", "R28Death...457")
event_col <- event_candidates[event_candidates %in% names(df)][1]
group_col <- "IAP-MAPgroup"

if (!time_col  %in% names(df)) stop("R28day column not found in input data.")
if (!group_col %in% names(df)) stop("IAP-MAPgroup column not found in input data.")
if (is.na(event_col))          stop("R28Death column not found in input data.")

# Labels A-D follow the manuscript's trajectory numbering (1-4)
km_df <- df |>
  dplyr::transmute(
    IAP.MAPgroup = factor(.data[[group_col]], levels = c(1, 2, 3, 4),
                          labels = c("D", "A", "C", "B")),
    R28day   = suppressWarnings(as.numeric(.data[[time_col]])),
    R28Death = suppressWarnings(as.numeric(.data[[event_col]]))
  ) |>
  dplyr::filter(!is.na(IAP.MAPgroup), !is.na(R28day), !is.na(R28Death))

## ---- log-rank test ----
fit       <- survival::survfit(survival::Surv(R28day, R28Death) ~ IAP.MAPgroup, data = km_df)
surv_diff <- survival::survdiff(survival::Surv(R28day, R28Death) ~ IAP.MAPgroup, data = km_df)
p_value   <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

summary_fit <- summary(fit)
km_points <- data.frame(
  strata   = summary_fit$strata,
  time     = summary_fit$time,
  n.risk   = summary_fit$n.risk,
  n.event  = summary_fit$n.event,
  survival = summary_fit$surv,
  std.err  = summary_fit$std.err,
  lower    = summary_fit$lower,
  upper    = summary_fit$upper
)
group_counts <- km_df |> dplyr::count(IAP.MAPgroup, name = "n")

## ---- plot ----
plot_obj <- survminer::ggsurvplot(
  fit, data = km_df,
  pval = TRUE, pval.coord = c(0, 0.15),
  risk.table = TRUE, risk.table.col = "strata",
  xlab = "Time (days)", ylab = "Cumulative hazard",
  legend.title = "IAP-MAPgroup", legend.labs = c("D", "A", "C", "B"),
  surv.median.line = "hv", ggtheme = ggplot2::theme_minimal(),
  palette = c("#e7298a", "#1b9e77", "#7570b3", "#d95f02"),
  risk.table.height = 0.3, risk.table.y.text.col = TRUE,
  risk.table.title = "Number at risk: n (%)",
  font.risk.table = c(12, "plain", "black"),
  font.legend = c(12, "plain", "black"),
  break.time.by = 7, xlim = c(0, 28), ylim = c(0, 1),
  conf.int = TRUE, risk.table.y.text = FALSE, legend = c(0.85, 0.25)
)

pdf_file <- file.path(output_dir, "KM_curve_R28Death.pdf")
png_file <- file.path(output_dir, "KM_curve_R28Death.png")
risk_pdf <- file.path(output_dir, "KM_curve_R28Death_with_risktable.pdf")
risk_png <- file.path(output_dir, "KM_curve_R28Death_with_risktable.png")

grDevices::pdf(risk_pdf, width = 14, height = 8); print(plot_obj); grDevices::dev.off()
ggplot2::ggsave(risk_png, plot = survminer::arrange_ggsurvplots(list(plot_obj), print = FALSE),
                width = 14, height = 8, dpi = 300)
grDevices::pdf(pdf_file, width = 14, height = 6); print(plot_obj$plot); grDevices::dev.off()
ggplot2::ggsave(png_file, plot = plot_obj$plot, width = 14, height = 6, dpi = 300)

## ---- export ----
stats_df <- data.frame(chisq = surv_diff$chisq, df = length(surv_diff$n) - 1,
                       p_value = p_value, n_used = nrow(km_df))

write.csv(km_points,    file.path(output_dir, "KM_curve_points.csv"),  row.names = FALSE, fileEncoding = "UTF-8")
write.csv(group_counts, file.path(output_dir, "KM_group_counts.csv"),  row.names = FALSE, fileEncoding = "UTF-8")
write.csv(stats_df,     file.path(output_dir, "KM_logrank_stats.csv"), row.names = FALSE, fileEncoding = "UTF-8")
writexl::write_xlsx(list(km_points = km_points, group_counts = group_counts, logrank_stats = stats_df),
                    path = file.path(output_dir, "KM_results.xlsx"))
writeLines(capture.output(print(summary(fit))), file.path(output_dir, "KM_summary.txt"))

cat("Log-rank p =", format(p_value, digits = 6), "\n")
cat("Saved to:", output_dir, "\n")
