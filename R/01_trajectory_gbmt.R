# ============================================================
# 01 | Group-based multivariate trajectory modelling (GBMT)
# Joint MAP + IAP trajectories, development cohort
# Fits 2-5 group solutions and exports fit indices (AIC/BIC/CIC),
# average posterior probabilities (APPA), OCC, group proportions,
# trajectory plots, and per-subject group assignments.
#
# Input : data/dev_longitudinal_map_iap.csv   (long format; see data/README.md)
# Output: output/01_trajectory_gbmt/
# ============================================================

options(stringsAsFactors = FALSE)

## ---- packages ----
pkgs <- c("readr", "dplyr", "tidyr", "gbmt", "writexl", "ggplot2", "tibble")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

## ---- paths ----
input_file <- "data/dev_longitudinal_map_iap.csv"
output_dir <- "output/01_trajectory_gbmt"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ---- helpers ----
calculate_OCC <- function(appa, prior) {
  (appa / (1 - appa)) / (prior / (1 - prior))
}

pad_named_vector <- function(x, n, prefix) {
  out <- rep(NA_real_, n)
  out[seq_along(x)] <- as.numeric(x)
  names(out) <- paste0(prefix, seq_len(n))
  out
}

extract_ic <- function(model) {
  # Manual AIC / BIC / CIC for gbmt models
  logLik <- if (!is.null(model$logLik)) model$logLik else NA_real_
  n  <- length(model$assign)                 # number of subjects
  ng <- length(unique(model$assign))
  n_vars   <- length(model$x.names)
  n_params <- ng * n_vars * (model$d + 1) + (ng - 1)

  if (is.na(logLik)) return(data.frame(AIC = NA_real_, BIC = NA_real_, CIC = NA_real_))

  AIC <- -2 * logLik + 2 * n_params
  BIC <- -2 * logLik + n_params * log(n)
  CIC <- -2 * logLik + n_params * (log(n) + 1)
  data.frame(AIC = AIC, BIC = BIC, CIC = CIC)
}

build_metrics_row <- function(model, ng) {
  model$d <- 3
  ic <- extract_ic(model)
  max_groups <- 5
  appa  <- model$appa
  prior <- model$prior
  occ   <- rep(NA_real_, length(appa))
  valid <- !is.na(appa) & !is.na(prior) & appa < 1 & prior < 1 & appa > 0 & prior > 0
  occ[valid] <- calculate_OCC(appa[valid], prior[valid])
  proportions <- as.numeric(prop.table(table(model$assign)))
  cbind(
    data.frame(
      ng = ng,
      degree = 3,
      n_subjects = length(model$assign),
      logLik = if (!is.null(model$logLik)) model$logLik else NA_real_
    ),
    ic,
    as.data.frame(as.list(pad_named_vector(appa, max_groups, "AvePP_group"))),
    as.data.frame(as.list(pad_named_vector(occ, max_groups, "OCC_group"))),
    as.data.frame(as.list(pad_named_vector(prior, max_groups, "Prior_group"))),
    as.data.frame(as.list(pad_named_vector(proportions, max_groups, "Prop_group")))
  )
}

plot_model <- function(model, ng, output_dir) {
  png_file <- file.path(output_dir, paste0("joint_gbmt_", ng, "_groups.png"))
  pdf_file <- file.path(output_dir, paste0("joint_gbmt_", ng, "_groups.pdf"))
  cols <- c("#0099B4B2", "#D7191C", "#7FBC41", "#FFD92F", "#984EA3")[seq_len(ng)]
  png(filename = png_file, width = 1800, height = 800, res = 150)
  par(lwd = 1)
  plot(model, mar = c(4, 7, 2, 5), bands = FALSE, xlab = "Days", ylab = "mmHg",
       titles = c("MeanArterialPressure", "IAP"),
       add.grid = TRUE, add.legend = TRUE, col = cols)
  dev.off()
  pdf(file = pdf_file, width = 12, height = 5.5)
  par(lwd = 1)
  plot(model, mar = c(4, 7, 2, 5), bands = FALSE, xlab = "Days", ylab = "mmHg",
       titles = c("MeanArterialPressure", "IAP"),
       add.grid = TRUE, add.legend = TRUE, col = cols)
  dev.off()
  data.frame(ng = ng, plot_png = png_file, plot_pdf = pdf_file)
}

## ---- read data ----
df <- readr::read_csv(input_file, show_col_types = FALSE)
names(df) <- c("ID", "Time", "MeanArterialPressure", "IAP")
df <- as.data.frame(df) |>
  mutate(
    ID = as.character(ID),
    Time = as.numeric(Time),
    MeanArterialPressure = as.numeric(MeanArterialPressure),
    IAP = as.numeric(IAP)
  ) |>
  arrange(ID, Time)

dup_df <- df |> count(ID, Time, name = "n") |> filter(n > 1)
if (nrow(dup_df) > 0) stop("Duplicate ID-Time combinations found.")

## ---- fit 2-5 group models ----
var_names       <- c("MeanArterialPressure", "IAP")
fitted_models   <- list()
metrics_list    <- list()
plot_files      <- list()
summary_text    <- character()
assignment_list <- list()
posterior_list  <- list()

for (g in 2:5) {
  model <- gbmt(x.names = var_names, unit = "ID", time = "Time",
                d = 3, ng = g, data = df, scaling = 0)
  fitted_models[[as.character(g)]] <- model
  metrics_list[[as.character(g)]]  <- build_metrics_row(model, g)
  plot_files[[as.character(g)]]    <- plot_model(model, g, output_dir)
  summary_text <- c(summary_text, paste0("===== ", g, " groups ====="),
                    paste(capture.output(summary(model)), collapse = "\n"), "")

  assign_df <- tibble::tibble(ID = names(model$assign), group = as.integer(model$assign))
  names(assign_df)[2] <- paste0("group_", g)
  assignment_list[[paste0("assign_", g, "_groups")]] <- assign_df

  post_df <- tibble::as_tibble(model$postprob, rownames = "ID")
  posterior_list[[paste0("posterior_", g, "_groups")]] <- post_df
}

metrics_all <- dplyr::bind_rows(metrics_list)
plots_df    <- dplyr::bind_rows(plot_files)

ic_table        <- metrics_all |> select(ng, degree, n_subjects, logLik, AIC, BIC, CIC)
avepp_table     <- metrics_all |> select(ng, starts_with("AvePP_group"))
occ_table       <- metrics_all |> select(ng, starts_with("OCC_group"))
prior_table     <- metrics_all |> select(ng, starts_with("Prior_group"))
proportion_table <- metrics_all |> select(ng, starts_with("Prop_group"))

## ---- export ----
writeLines(summary_text, file.path(output_dir, "joint_gbmt_model_summaries.txt"))
readr::write_csv(metrics_all,      file.path(output_dir, "joint_gbmt_metrics_all.csv"), na = "")
readr::write_csv(ic_table,         file.path(output_dir, "joint_gbmt_ic_table.csv"), na = "")
readr::write_csv(avepp_table,      file.path(output_dir, "joint_gbmt_avepp_table.csv"), na = "")
readr::write_csv(occ_table,        file.path(output_dir, "joint_gbmt_occ_table.csv"), na = "")
readr::write_csv(prior_table,      file.path(output_dir, "joint_gbmt_prior_table.csv"), na = "")
readr::write_csv(proportion_table, file.path(output_dir, "joint_gbmt_group_proportions.csv"), na = "")

excel_sheets <- c(
  list(
    metrics_all      = metrics_all,
    ic_table         = ic_table,
    avepp_table      = avepp_table,
    occ_table        = occ_table,
    prior_table      = prior_table,
    proportion_table = proportion_table,
    plot_files       = plots_df
  ),
  assignment_list,
  posterior_list
)
writexl::write_xlsx(excel_sheets, path = file.path(output_dir, "joint_gbmt_results.xlsx"))
saveRDS(fitted_models, file.path(output_dir, "joint_gbmt_models_2_to_5_groups.rds"))

cat("Saved to:", output_dir, "\n")
