## Demo: cond.correct package on Okataina and Rotorua profiler data
## Demonstrates the full workflow: analyse → detect eras → optimise → correct → plot

library(dplyr)
library(ggplot2)
library(patchwork)

# Source package functions
pkg_dir <- "C:/Users/chris/Dropbox/limnotrack/software/cond_correct/R"
for (f in list.files(pkg_dir, full.names = TRUE)) source(f)

data_dir <- "C:/Users/chris/Dropbox/limnotrack/dropbox/wqprofiler/data/raw"
fig_dir  <- "C:/Users/chris/Dropbox/limnotrack/software/cond_correct/figures"
dir.create(fig_dir, showWarnings = FALSE)


# =============================================================================
# OKATAINA
# =============================================================================
cat("\n========== OKATAINA ==========\n")

raw_ok <- readRDS(file.path(data_dir, "Okataina.rds"))
pro_ok <- raw_ok[["pro"]] |>
  filter(DateTime >= as.POSIXct("2024-12-09", tz = "Etc/GMT-12"),
         DptSns > 0, DptSns < 60,
         !is.na(Cond), Cond > 0, Cond < 9000)

# Step 1: Analyse all profiles
cat("Analysing profiles...\n")
t0 <- Sys.time()
analysis_ok <- cond_analyse(pro_ok, col_cond = "Cond", col_temp = "TmpWtr", delta_t = 0.2)
cat("  Done in", round(difftime(Sys.time(), t0, units = "secs")), "seconds\n")
cat("  Profiles:", nrow(analysis_ok), "| Fit OK:", sum(analysis_ok$fit_ok), "\n")
cat("  MLD available:", sum(!is.na(analysis_ok$mld)), "\n")
cat("  A median:", round(median(analysis_ok$A[analysis_ok$fit_ok], na.rm = TRUE), 1), "\n")
cat("  tau median:", round(median(analysis_ok$tau[analysis_ok$fit_ok], na.rm = TRUE), 1), "\n")

# Step 2: Detect eras (parameter change periods)
analysis_ok <- cond_detect_eras(analysis_ok)
cat("  Eras:", length(unique(analysis_ok$era)), "\n")

# Step 3: Optimise correction parameters per era
cat("Optimising...\n")
pid_ok <- detect_profiles(pro_ok$DateTime, pro_ok$DptSns)
t0 <- Sys.time()
opt_ok <- cond_optimise(pro_ok, analysis_ok, col_cond = "Cond", col_temp = "TmpWtr",
                        profile_id = pid_ok, model = "double")
cat("  Done in", round(difftime(Sys.time(), t0, units = "secs")), "seconds\n")
cat("  Global tau:", opt_ok$tau, "m (from", opt_ok$tau_n, "deep-mixed profiles)\n")
cat("  Optimised A per era:\n")
params_ok <- opt_ok$params
print(params_ok)

# Step 4: Apply correction (era-specific A, global tau)
era_map_ok <- analysis_ok$era[match(pid_ok, analysis_ok$profile_id)]
pro_ok$SpCond_corrected <- cond_correct(pro_ok$Cond, pro_ok$DptSns,
                                         params = params_ok, era = era_map_ok)

# Assess improvement
cat("\nDepth uniformity comparison:\n")
cat("  Raw Cond:\n")
pro_ok |> mutate(dband = cut(DptSns, breaks = c(0,5,10,20,30,40,60))) |>
  group_by(dband) |>
  summarise(med = round(median(Cond), 1), iqr = round(IQR(Cond), 1), .groups = "drop") |>
  print()
cat("  Corrected SpCond:\n")
pro_ok |> filter(!is.na(SpCond_corrected)) |>
  mutate(dband = cut(DptSns, breaks = c(0,5,10,20,30,40,60))) |>
  group_by(dband) |>
  summarise(med = round(median(SpCond_corrected), 1), iqr = round(IQR(SpCond_corrected), 1), .groups = "drop") |>
  print()

# Step 5: Diagnostic plots
cat("\nGenerating plots...\n")

# Diagnostics panel
p_diag <- cond_plot_diagnostics(analysis_ok, opt_ok)
print(p_diag)
ggsave(file.path(fig_dir, "okataina_diagnostics.png"), p_diag, width = 12, height = 10, dpi = 150)

# Example profiles
p_prof <- cond_plot_profiles(pro_ok, analysis_ok, opt_ok,
                              col_cond = "Cond", profile_id = pid_ok, n_profiles = 16)
print(p_prof)
ggsave(file.path(fig_dir, "okataina_profiles.png"), p_prof, width = 14, height = 10, dpi = 150)

# Before/after dot plot
p_raw <- ggplot(pro_ok, aes(x = DateTime, y = Cond, colour = DptSns)) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_colour_viridis_c(option = "C", direction = -1, name = "Depth (m)") +
  labs(x = NULL, y = "Raw Cond (µS/cm)", title = "Okataina — Raw conductivity") +
  theme_bw(base_size = 9)

p_corr <- ggplot(pro_ok |> filter(!is.na(SpCond_corrected)),
                  aes(x = DateTime, y = SpCond_corrected, colour = DptSns)) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_colour_viridis_c(option = "C", direction = -1, name = "Depth (m)") +
  labs(x = NULL, y = "Corrected SpCond (µS/cm)", title = "Okataina — Corrected") +
  theme_bw(base_size = 9)

p_compare <- p_raw / p_corr
print(p_compare)
ggsave(file.path(fig_dir, "okataina_before_after.png"), p_compare, width = 14, height = 8, dpi = 150)
cat("  Okataina plots saved.\n")


# =============================================================================
# OKATAINA — verify QC script uses the same correction
# =============================================================================
cat("\n--- Verifying qc1_Okataina.R correction matches ---\n")
Sys.setenv(tz = "Etc/GMT-12")
qc_funs <- list.files("C:/Users/chris/Dropbox/limnotrack/dropbox/wqprofiler/scripts_qc/funs_qc",
                       pattern = "*.R", full.names = TRUE)
invisible(lapply(qc_funs, source))
source("C:/Users/chris/Dropbox/limnotrack/software/wqprofiler-data/fun/derived_vars.R")
source("C:/Users/chris/Dropbox/limnotrack/dropbox/wqprofiler/scripts_qc/qc1_Okataina.R")
full_ok <- readRDS("C:/Users/chris/Dropbox/limnotrack/dropbox/wqprofiler/data/full/full_Okataina.rds")
qc_ok <- qc1_site(raw_ok, full_ok[["meta"]])
pro_qc_ok <- qc_ok[["pro"]] |> filter(DptSns > 0, DptSns < 60, !is.na(SpCond))
cat("  QC script SpCond depth uniformity:\n")
pro_qc_ok |> mutate(dband = cut(DptSns, breaks = c(0,5,10,20,30,40,60))) |>
  group_by(dband) |>
  summarise(med = round(median(SpCond), 1), iqr = round(IQR(SpCond), 1), .groups = "drop") |>
  print()


# =============================================================================
# ROTORUA
# =============================================================================
cat("\n========== ROTORUA ==========\n")

raw_ro <- readRDS(file.path(data_dir, "Rotorua.rds"))
pro_ro <- raw_ro[["pro"]]

# Check what conductivity column exists
cond_col <- if ("SpCond" %in% names(pro_ro) && any(!is.na(pro_ro$SpCond))) {
  "SpCond"
} else if ("Cond" %in% names(pro_ro) && any(!is.na(pro_ro$Cond))) {
  "Cond"
} else {
  NULL
}

if (is.null(cond_col)) {
  cat("  No conductivity data found for Rotorua. Trying other column names...\n")
  cond_cols <- grep("cond|Cond|SpCond|cnd", names(pro_ro), value = TRUE, ignore.case = TRUE)
  cat("  Available:", paste(cond_cols, collapse = ", "), "\n")
  if (length(cond_cols) > 0) {
    for (cc in cond_cols) {
      n <- sum(!is.na(pro_ro[[cc]]))
      cat("  ", cc, ": n_valid =", n, "\n")
      if (n > 1000) { cond_col <- cc; break }
    }
  }
}

if (!is.null(cond_col)) {
  pro_ro <- pro_ro |>
    filter(DptSns > 0, DptSns < 50, !is.na(.data[[cond_col]]),
           .data[[cond_col]] > 0, .data[[cond_col]] < 9000)

  cat("  Using column:", cond_col, "\n")
  cat("  Rows:", nrow(pro_ro), "\n")
  cat("  Date range:", format(range(pro_ro$DateTime, na.rm = TRUE)), "\n")

  # Analyse
  cat("  Analysing profiles...\n")
  t0 <- Sys.time()
  analysis_ro <- cond_analyse(pro_ro, col_cond = cond_col, col_temp = "TmpWtr", delta_t = 0.2)
  cat("  Done in", round(difftime(Sys.time(), t0, units = "secs")), "seconds\n")
  cat("  Profiles:", nrow(analysis_ro), "| Fit OK:", sum(analysis_ro$fit_ok), "\n")
  cat("  A median:", round(median(analysis_ro$A[analysis_ro$fit_ok], na.rm = TRUE), 1), "\n")
  cat("  tau median:", round(median(analysis_ro$tau[analysis_ro$fit_ok], na.rm = TRUE), 1), "\n")

  # Detect eras
  analysis_ro <- cond_detect_eras(analysis_ro)
  cat("  Eras:", length(unique(analysis_ro$era)), "\n")

  # Optimise
  cat("  Optimising...\n")
  pid_ro <- detect_profiles(pro_ro$DateTime, pro_ro$DptSns)
  t0 <- Sys.time()
  opt_ro <- cond_optimise(pro_ro, analysis_ro, col_cond = cond_col, col_temp = "TmpWtr",
                          profile_id = pid_ro, model = "double")
  cat("  Done in", round(difftime(Sys.time(), t0, units = "secs")), "seconds\n")
  cat("  Global tau:", opt_ro$tau, "\n")
  params_ro <- opt_ro$params
  print(params_ro)

  # Correct
  era_map_ro <- analysis_ro$era[match(pid_ro, analysis_ro$profile_id)]
  pro_ro$SpCond_corrected <- cond_correct(pro_ro[[cond_col]], pro_ro$DptSns,
                                           params = params_ro, era = era_map_ro)

  # Depth uniformity
  cat("\n  Depth uniformity:\n  Raw:\n")
  pro_ro |> mutate(dband = cut(DptSns, breaks = c(0,5,10,20,30,40,50))) |>
    group_by(dband) |>
    summarise(med = round(median(.data[[cond_col]]), 1), .groups = "drop") |> print()
  cat("  Corrected:\n")
  pro_ro |> filter(!is.na(SpCond_corrected)) |>
    mutate(dband = cut(DptSns, breaks = c(0,5,10,20,30,40,50))) |>
    group_by(dband) |>
    summarise(med = round(median(SpCond_corrected), 1), .groups = "drop") |> print()

  # Plots
  p_diag_ro <- cond_plot_diagnostics(analysis_ro, opt_ro)
  print(p_diag_ro)
  ggsave(file.path(fig_dir, "rotorua_diagnostics.png"), p_diag_ro, width = 12, height = 8, dpi = 150)

  p_prof_ro <- cond_plot_profiles(pro_ro, analysis_ro, opt_ro,
                                   col_cond = cond_col, profile_id = pid_ro, n_profiles = 16)
  print(p_prof_ro)
  ggsave(file.path(fig_dir, "rotorua_profiles.png"), p_prof_ro, width = 14, height = 10, dpi = 150)

  p_raw_ro <- ggplot(pro_ro, aes(x = DateTime, y = .data[[cond_col]], colour = DptSns)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_colour_viridis_c(option = "C", direction = -1, name = "Depth (m)") +
    labs(x = NULL, y = paste("Raw", cond_col), title = "Rotorua — Raw") +
    theme_bw(base_size = 9)

  p_corr_ro <- ggplot(pro_ro |> filter(!is.na(SpCond_corrected)),
                       aes(x = DateTime, y = SpCond_corrected, colour = DptSns)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_colour_viridis_c(option = "C", direction = -1, name = "Depth (m)") +
    labs(x = NULL, y = "Corrected SpCond", title = "Rotorua — Corrected") +
    theme_bw(base_size = 9)

  p_compare_ro <- p_raw_ro / p_corr_ro
  print(p_compare_ro)
  ggsave(file.path(fig_dir, "rotorua_before_after.png"), p_compare_ro, width = 14, height = 8, dpi = 150)
  cat("  Rotorua plots saved.\n")
} else {
  cat("  No usable conductivity data found for Rotorua.\n")
}

cat("\n========== DONE ==========\n")
cat("Figures saved to:", fig_dir, "\n")
