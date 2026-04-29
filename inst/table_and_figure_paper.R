# =============================================================================
# table_and_figure_paper.R
#
# Generate tables and figures for the article "When Does Unequal
# Randomization Shorten Survival Trials? Closed-Form Bounds on the Required
# Accrual Rate".
#
# The script produces:
#   * Figures 1-2: relative efficiency curves with three vertical bounds for
#                  Phase 2 (HR = 0.6) and Phase 3 (HR = 0.8) settings;
#   * Figure 3:    signed relative gaps of all three analytical bounds from
#                  rho* across (HR, r, d);
#   * Figure 4:    the four quantities rho_L^univ, rho_L^tight, rho*, rho_U
#                  as functions of r, faceted by HR x d;
#   * Table:       application-section table for three hypothetical
#                  trial examples (CV outcome, Phase 3 oncology,
#                  Phase 2 oncology).
#
# A 54-scenario grid (Phase 2 / Phase 3 with HR in {0.6, 0.8}) is used to
# verify that rho_L^tight <= rho* <= rho_U holds across realistic design
# configurations. The verification result is reported in the console and
# referenced in the article body without a dedicated table.
#
# Outputs (EPS for compatibility with the existing LaTeX pipeline):
#   fig1_re_phase2.eps            Phase 2 RE curves with three vertical bounds
#   fig2_re_phase3.eps            Phase 3 RE curves with three vertical bounds
#   fig3_univ_accuracy.eps        Signed relative gaps of all three bounds
#   fig4_bounds_vs_HR.eps         Four quantities vs r, faceted by HR x d
#   table_hypothetical.tex        Application section table (3 hypothetical
#                                 trial examples)
#   table_hypothetical.csv        Same table in CSV form
# =============================================================================

# --- Packages ---
library(dplyr)
library(tidyr)
library(ggplot2)

# --- Load functions: prefer the installed survRho package, ---
# --- fall back to source() during package development.       ---
# Recommended: library(survRho).
# Otherwise (e.g., when developing the package), set the working
# directory to the folder containing the R function files and the
# script will source them directly. Set the working directory in
# RStudio via "Session > Set Working Directory > To Source File
# Location", or via setwd() (e.g., setwd("/path/to/survRho/R")).
if (requireNamespace("survRho", quietly = TRUE)) {
  library(survRho)
} else {
  source_files <- c(
    "calc_trial_duration.R",
    "calc_relative_efficiency.R",
    "calc_relative_efficiency_grid.R",
    "find_rho_for_target_re.R",
    "find_rho_for_target_re_grid.R",
    "calc_rho_bounds.R",
    "calc_rho_bounds_grid.R",
    "print.rho_bounds.R",
    "print.rho_bounds_grid.R",
    "print.re_grid.R",
    "build_accrual_table.R"
  )
  missing_files <- source_files[!file.exists(source_files)]
  if (length(missing_files) > 0) {
    stop(
      "The survRho package is not installed, and the following source ",
      "files were not found in the working directory ('", getwd(), "'):\n  ",
      paste(missing_files, collapse = "\n  "),
      "\nEither install the package (devtools::install_github(...)) or ",
      "set the working directory to the folder containing the R function ",
      "files (e.g., via 'setwd()' or RStudio's Session menu)."
    )
  }
  invisible(lapply(source_files, source))
}

# --- Output directory ---
# All figures and tables produced by this script are written to OUT_DIR.
# By default, outputs are written to a subdirectory "outputs" of the current
# working directory; the directory is created if it does not exist. To send
# outputs elsewhere, set OUT_DIR to any existing or new path before running
# the script (e.g., OUT_DIR <- "/path/to/figures").

# --- Common design parameters shared across phases ---
ALPHA <- 0.025
POWER <- 0.9
TIME_UNIT <- "months"

# Phase 2 setting (HR = 0.6, short MSTC, short follow-up, smaller accrual)
PHASE2_HR <- 0.6
PHASE2_MST <- c(6, 9, 12)
PHASE2_F <- 8
PHASE2_OMEGA1 <- 30

# Phase 3 setting (HR = 0.8, long MSTC, longer follow-up, larger accrual)
PHASE3_HR <- 0.8
PHASE3_MST <- c(12, 18, 24)
PHASE3_F <- 12
PHASE3_OMEGA1 <- 60

# Dropout rates evaluated
DROPOUT_VALS <- c(0, 0.15, 0.3)

# Allocation ratios for the figures and the main scenario grid
R_VALS <- c(2, 3, 4)

# Output directory
if (!exists("OUT_DIR")) {
  OUT_DIR <- file.path(getwd(), "outputs")
}
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# Paper-friendly figure size (in mm)
FIG_WIDTH_MM <- 160
FIG_HEIGHT_MM <- 110

# Common base font size
BASE_SIZE <- 12

# =============================================================================
# Helper: build parse-ready facet labels
# =============================================================================
make_HR_label <- function(x) {
  paste0("HR == ", format(x, nsmall = 1))
}
make_d_label <- function(x) {
  paste0("d == ", format(x, nsmall = 2))
}
make_mst_label <- function(x) {
  paste0("MST[C] == ", format(round(x, 0), nsmall = 0))
}

# =============================================================================
# Helper: solve rho* numerically for a single scenario
# =============================================================================
solve_rho_numerical <- function(HR, hazard_c, f, d, omega1, r) {
  a_k_vec <- 0
  omega_k_vec <- omega1
  res <- tryCatch(
    find_rho_for_target_re(
      HR = HR,
      hazard_c = hazard_c,
      a_k = a_k_vec,
      omega_k = omega_k_vec,
      f = f,
      d_t = d,
      d_c = d,
      r = r,
      target_re = 1,
      metric = "tau",
      alpha = ALPHA,
      power = POWER,
      time_unit = TIME_UNIT,
      formula = "Schoenfeld"
    ),
    error = function(e) NULL
  )
  if (is.null(res)) return(NA_real_)
  if (!is.null(res$rho)) return(as.numeric(res$rho))
  if (!is.null(res$achieved_rho)) return(as.numeric(res$achieved_rho))
  if (is.numeric(res)) return(as.numeric(res))
  NA_real_
}

# =============================================================================
# Build the master scenario data frame across (phase, HR, MST_C, d, r)
# =============================================================================
# Build the 54-scenario grid used for the consistency check spanning
# Phase 2 (HR = 0.6) and Phase 3 (HR = 0.8) settings.
phase2_grid <- expand.grid(phase = "Phase 2", HR = PHASE2_HR,
                           mst_c = PHASE2_MST, f = PHASE2_F,
                           omega1 = PHASE2_OMEGA1,
                           d = DROPOUT_VALS, r = R_VALS,
                           stringsAsFactors = FALSE)
phase3_grid <- expand.grid(phase = "Phase 3", HR = PHASE3_HR,
                           mst_c = PHASE3_MST, f = PHASE3_F,
                           omega1 = PHASE3_OMEGA1,
                           d = DROPOUT_VALS, r = R_VALS,
                           stringsAsFactors = FALSE)

scenario_df <- dplyr::bind_rows(phase2_grid, phase3_grid)
scenario_df$hazard_c <- log(2) / scenario_df$mst_c

scenario_df$rho_numerical <- NA_real_
scenario_df$rho_lower_loose <- NA_real_
scenario_df$rho_lower_tight <- NA_real_
scenario_df$rho_upper <- NA_real_
scenario_df$tau_1 <- NA_real_

for (i in seq_len(nrow(scenario_df))) {
  sc <- scenario_df[i, ]
  scenario_df$rho_numerical[i] <- solve_rho_numerical(
    HR = sc$HR, hazard_c = sc$hazard_c, f = sc$f,
    d = sc$d, omega1 = sc$omega1, r = sc$r
  )
  bnd <- calc_rho_bounds(
    HR = sc$HR, hazard_c = sc$hazard_c,
    d = sc$d, r = sc$r, time_unit = TIME_UNIT
  )
  scenario_df$rho_lower_loose[i] <- bnd$bounds_df$rho_lower_loose
  scenario_df$rho_lower_tight[i] <- bnd$bounds_df$rho_lower_tight
  scenario_df$rho_upper[i] <- bnd$bounds_df$rho_upper

  # Baseline tau_1 is informative for the RE plot annotation
  td_baseline <- tryCatch(
    calc_trial_duration(
      HR = sc$HR, hazard_c = sc$hazard_c,
      a_k = 0, omega_k = sc$omega1,
      f = sc$f, d_t = sc$d, d_c = sc$d,
      r = 1, alpha = ALPHA, power = POWER,
      time_unit = TIME_UNIT, formula = "Schoenfeld"
    ),
    error = function(e) NULL
  )
  if (!is.null(td_baseline) && !is.null(td_baseline$tau)) {
    scenario_df$tau_1[i] <- td_baseline$tau
  }
}

TOL <- 1e-6
scenario_df$within_bounds <- with(
  scenario_df,
  rho_numerical >= rho_lower_tight - TOL &
    rho_numerical <= rho_upper + TOL
)
n_total <- nrow(scenario_df)
n_pass <- sum(scenario_df$within_bounds, na.rm = TRUE)
message(sprintf(
  "Consistency check: %d / %d scenarios satisfy rho_L^tight <= rho* <= rho_U (tolerance = %.0e)",
  n_pass, n_total, TOL
))
if (n_pass < n_total) {
  failed <- scenario_df[!scenario_df$within_bounds, ]
  message("Scenarios failing the bounds check:")
  print(failed)
}

# =============================================================================
# Figure 1 / Figure 2: Relative efficiency curves with three vertical bounds
# =============================================================================
build_phase_curves <- function(phase_HR, phase_MST, phase_F, phase_OMEGA1) {
  rho_seq <- seq(1, 2, length.out = 101)
  r_and_rho <- expand.grid(r = R_VALS, rho = rho_seq,
                           stringsAsFactors = FALSE)
  results <- calc_relative_efficiency_grid(
    HR = phase_HR,
    hazard_c = log(2) / phase_MST,
    f = phase_F,
    d = DROPOUT_VALS,
    a_k = 0,
    omega_k = phase_OMEGA1,
    r = r_and_rho$r,
    rho = r_and_rho$rho,
    alpha = ALPHA,
    power = POWER,
    time_unit = TIME_UNIT,
    formula = "Schoenfeld"
  )
  results
}

build_phase_bounds <- function(phase_HR, phase_MST) {
  grid_df <- expand.grid(mst_c = phase_MST, d = DROPOUT_VALS, r = R_VALS,
                         stringsAsFactors = FALSE)
  grid_df$rho_lower_loose <- NA_real_
  grid_df$rho_lower_tight <- NA_real_
  grid_df$rho_upper <- NA_real_
  for (i in seq_len(nrow(grid_df))) {
    bnd <- calc_rho_bounds(
      HR = phase_HR,
      hazard_c = log(2) / grid_df$mst_c[i],
      d = grid_df$d[i],
      r = grid_df$r[i],
      time_unit = TIME_UNIT
    )
    grid_df$rho_lower_loose[i] <- bnd$bounds_df$rho_lower_loose
    grid_df$rho_lower_tight[i] <- bnd$bounds_df$rho_lower_tight
    grid_df$rho_upper[i] <- bnd$bounds_df$rho_upper
  }
  long_df <- dplyr::bind_rows(
    data.frame(grid_df[, c("mst_c", "d", "r")],
               bound_type = "Universal lower",
               rho_value = grid_df$rho_lower_loose,
               stringsAsFactors = FALSE),
    data.frame(grid_df[, c("mst_c", "d", "r")],
               bound_type = "Tight lower",
               rho_value = grid_df$rho_lower_tight,
               stringsAsFactors = FALSE),
    data.frame(grid_df[, c("mst_c", "d", "r")],
               bound_type = "Upper",
               rho_value = grid_df$rho_upper,
               stringsAsFactors = FALSE)
  )
  long_df$bound_type <- factor(
    long_df$bound_type,
    levels = c("Universal lower", "Tight lower", "Upper")
  )
  long_df
}

plot_phase_with_three_bounds <- function(phase_result, bounds_df,
                                         out_filename) {
  re_df <- phase_result$re_df
  re_df$mst_label <- factor(
    make_mst_label(re_df$mst_c),
    levels = make_mst_label(sort(unique(re_df$mst_c)))
  )
  re_df$d_label <- factor(
    make_d_label(re_df$d),
    levels = make_d_label(sort(unique(re_df$d)))
  )
  re_df$r_factor <- factor(re_df$r)

  bounds_df$mst_label <- factor(
    make_mst_label(bounds_df$mst_c),
    levels = make_mst_label(sort(unique(bounds_df$mst_c)))
  )
  bounds_df$d_label <- factor(
    make_d_label(bounds_df$d),
    levels = make_d_label(sort(unique(bounds_df$d)))
  )
  bounds_df$r_factor <- factor(bounds_df$r)

  baseline_df <- phase_result$original_df |>
    dplyr::filter(.data$r == 1 & .data$rho == 1) |>
    dplyr::mutate(
      mst_c = log(2) / .data$hazard_c,
      mst_label = make_mst_label(.data$mst_c),
      d_label = make_d_label(.data$d),
      tau_label = paste0("tau[1] == ",
                         format(round(.data$tau, 1), nsmall = 1))
    )

  # Compute rho* (intersection of RE curve with RE = 1) for each panel
  # by linear interpolation along the simulated rho grid.
  rho_star_df <- re_df |>
    dplyr::group_by(.data$mst_c, .data$d, .data$r,
                    .data$mst_label, .data$d_label, .data$r_factor) |>
    dplyr::summarise(
      rho_star = stats::approx(x = .data$re_tau, y = .data$rho,
                               xout = 1, ties = "ordered")$y,
      .groups = "drop"
    ) |>
    dplyr::filter(!is.na(.data$rho_star))

  r_levels <- sort(unique(re_df$r))
  r_expr <- lapply(r_levels, function(v) bquote(r == .(v)))

  p <- ggplot(re_df,
              aes(x = rho, y = re_tau,
                  color = r_factor, linetype = r_factor)) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 1,
               linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_vline(data = subset(bounds_df, bound_type == "Universal lower"),
               aes(xintercept = rho_value, color = r_factor),
               linetype = "dashed", linewidth = 0.35, alpha = 0.6,
               show.legend = FALSE) +
    geom_vline(data = subset(bounds_df, bound_type == "Tight lower"),
               aes(xintercept = rho_value, color = r_factor),
               linetype = "solid", linewidth = 0.35, alpha = 0.7,
               show.legend = FALSE) +
    geom_vline(data = subset(bounds_df, bound_type == "Upper"),
               aes(xintercept = rho_value, color = r_factor),
               linetype = "dotted", linewidth = 0.45, alpha = 0.9,
               show.legend = FALSE) +
    geom_point(data = rho_star_df,
               aes(x = rho_star, y = 1, color = r_factor),
               inherit.aes = FALSE,
               size = 2.2, shape = 16, show.legend = FALSE) +
    geom_text(data = baseline_df,
              aes(x = 1.02, y = Inf, label = tau_label),
              inherit.aes = FALSE, parse = TRUE,
              hjust = 0, vjust = 1.3, size = 3.6,
              color = "grey30") +
    facet_grid(rows = vars(d_label), cols = vars(mst_label),
               labeller = label_parsed) +
    scale_color_manual(
      name = NULL,
      values = c("2" = "#c77cff", "3" = "#1f77b4", "4" = "#2ca02c"),
      labels = r_expr
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c("2" = "solid", "3" = "dashed", "4" = "dotted"),
      labels = r_expr
    ) +
    scale_x_continuous(limits = c(1, 2)) +
    labs(x = expression(rho),
         y = "Relative efficiency") +
    guides(
      color    = guide_legend(keywidth = unit(1.0, "cm")),
      linetype = guide_legend(keywidth = unit(1.0, "cm"))
    ) +
    theme_bw(base_size = BASE_SIZE) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.5, "cm"),
          strip.background = element_rect(fill = "grey92"))

  ggsave(file.path(OUT_DIR, out_filename),
         plot = p,
         width = FIG_WIDTH_MM, height = FIG_HEIGHT_MM * 1.5,
         units = "mm", device = cairo_ps, fallback_resolution = 600)
  invisible(p)
}

phase2_curves <- build_phase_curves(
  phase_HR = PHASE2_HR,
  phase_MST = PHASE2_MST,
  phase_F = PHASE2_F,
  phase_OMEGA1 = PHASE2_OMEGA1
)
phase2_bounds <- build_phase_bounds(
  phase_HR = PHASE2_HR,
  phase_MST = PHASE2_MST
)
plot_phase_with_three_bounds(
  phase_result = phase2_curves,
  bounds_df = phase2_bounds,
  out_filename = "fig1_re_phase2.eps"
)

phase3_curves <- build_phase_curves(
  phase_HR = PHASE3_HR,
  phase_MST = PHASE3_MST,
  phase_F = PHASE3_F,
  phase_OMEGA1 = PHASE3_OMEGA1
)
phase3_bounds <- build_phase_bounds(
  phase_HR = PHASE3_HR,
  phase_MST = PHASE3_MST
)
plot_phase_with_three_bounds(
  phase_result = phase3_curves,
  bounds_df = phase3_bounds,
  out_filename = "fig2_re_phase3.eps"
)

# =============================================================================
# Figure 3: Signed relative gaps of all three bounds from rho*
#
# Quantifies, for each (HR, r, d, MST_C), the signed relative gap of
# each of the three analytical bounds from rho*, all using the common
# convention (rho_quantity - rho*) / (rho* - 1) * 100:
#     gap_univ  < 0  because rho_L^univ  <= rho*
#     gap_tight < 0  because rho_L^tight <= rho*
#     gap_upper > 0  because rho_U       >= rho*
# A horizontal line at y = 0 marks rho* itself. The closer the line is to
# zero, the tighter the analytical bound; the upper bound rho_U is
# interpreted as a sufficient condition (any accrual increase beyond rho_U
# guarantees tau_r < tau_1).
# =============================================================================
hr_seq <- seq(0.5, 0.95, by = 0.025)
r_seq_fig3 <- c(2, 3, 4, 5)
mst_for_fig3 <- 12
f_for_fig3 <- 12
omega1_for_fig3 <- 50

fig3_grid <- expand.grid(HR = hr_seq, r = r_seq_fig3,
                         d = DROPOUT_VALS,
                         stringsAsFactors = FALSE)
fig3_grid$rho_numerical <- NA_real_
fig3_grid$rho_lower_loose <- NA_real_
fig3_grid$rho_lower_tight <- NA_real_
fig3_grid$rho_upper <- NA_real_
for (i in seq_len(nrow(fig3_grid))) {
  hazard_c_i <- log(2) / mst_for_fig3
  fig3_grid$rho_numerical[i] <- solve_rho_numerical(
    HR = fig3_grid$HR[i], hazard_c = hazard_c_i,
    f = f_for_fig3, d = fig3_grid$d[i],
    omega1 = omega1_for_fig3, r = fig3_grid$r[i]
  )
  bnd <- calc_rho_bounds(
    HR = fig3_grid$HR[i], hazard_c = hazard_c_i,
    d = fig3_grid$d[i], r = fig3_grid$r[i],
    time_unit = TIME_UNIT
  )
  fig3_grid$rho_lower_loose[i] <- bnd$bounds_df$rho_lower_loose
  fig3_grid$rho_lower_tight[i] <- bnd$bounds_df$rho_lower_tight
  fig3_grid$rho_upper[i] <- bnd$bounds_df$rho_upper
}
fig3_grid$gap_univ <- with(
  fig3_grid,
  (rho_lower_loose - rho_numerical) / (rho_numerical - 1) * 100
)
fig3_grid$gap_tight <- with(
  fig3_grid,
  (rho_lower_tight - rho_numerical) / (rho_numerical - 1) * 100
)
fig3_grid$gap_upper <- with(
  fig3_grid,
  (rho_upper - rho_numerical) / (rho_numerical - 1) * 100
)
fig3_grid$d_label <- factor(
  make_d_label(fig3_grid$d),
  levels = make_d_label(sort(unique(fig3_grid$d)))
)

# Long format for faceting by quantity
fig3_long <- dplyr::bind_rows(
  data.frame(fig3_grid[, c("HR", "r", "d", "d_label")],
             quantity = "univ",
             gap = fig3_grid$gap_univ,
             stringsAsFactors = FALSE),
  data.frame(fig3_grid[, c("HR", "r", "d", "d_label")],
             quantity = "tight",
             gap = fig3_grid$gap_tight,
             stringsAsFactors = FALSE),
  data.frame(fig3_grid[, c("HR", "r", "d", "d_label")],
             quantity = "upper",
             gap = fig3_grid$gap_upper,
             stringsAsFactors = FALSE)
)
fig3_long$quantity_label <- factor(
  fig3_long$quantity,
  levels = c("univ", "tight", "upper"),
  labels = c("rho[L]^'univ'", "rho[L]^'tight'", "rho[U]")
)

p_fig3 <- ggplot(fig3_long,
                 aes(x = HR, y = gap,
                     color = factor(r), shape = factor(r))) +
  geom_hline(yintercept = 0, color = "grey50",
             linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_grid(rows = vars(quantity_label), cols = vars(d_label),
             labeller = label_parsed, scales = "free_y") +
  scale_color_manual(
    name = NULL,
    values = c("2" = "#c77cff", "3" = "#1f77b4",
               "4" = "#2ca02c", "5" = "#d62728"),
    labels = lapply(r_seq_fig3, function(v) bquote(r == .(v)))
  ) +
  scale_shape_manual(
    name = NULL,
    values = c("2" = 16, "3" = 17, "4" = 15, "5" = 18),
    labels = lapply(r_seq_fig3, function(v) bquote(r == .(v)))
  ) +
  scale_x_continuous(breaks = seq(0.5, 0.95, by = 0.1)) +
  labs(x = expression("Hazard ratio " * HR),
       y = expression(atop("Signed relative gap",
                           "from " * rho * "* (%)"))) +
  guides(
    color = guide_legend(keywidth = unit(1.0, "cm")),
    shape = guide_legend(keywidth = unit(1.0, "cm"))
  ) +
  theme_bw(base_size = BASE_SIZE) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.5, "cm"),
        strip.background = element_rect(fill = "grey92"))

ggsave(file.path(OUT_DIR, "fig3_univ_accuracy.eps"),
       plot = p_fig3,
       width = FIG_WIDTH_MM, height = FIG_HEIGHT_MM * 1.5,
       units = "mm", device = cairo_ps, fallback_resolution = 600)

# =============================================================================
# Figure 4: rho_L^univ, rho_L^tight, rho*, rho_U across (r, HR, d)
#
# At fixed (MST_C, follow-up, omega_1), shows how the four quantities behave
# as functions of r. Faceted by HR (rows) and d (cols) so that each
# panel contains exactly four lines and the convergence of the three bounds
# to rho* as HR -> 1 is visible across rows.
# =============================================================================
r_seq_fig4 <- seq(2, 5, by = 0.05)
hr_set_fig4 <- c(0.5, 0.6, 0.7, 0.8, 0.9)
mst_for_fig4 <- 12
f_for_fig4 <- 12
omega1_for_fig4 <- 50

fig4_grid <- expand.grid(HR = hr_set_fig4, r = r_seq_fig4,
                         d = DROPOUT_VALS,
                         stringsAsFactors = FALSE)
fig4_grid$rho_numerical <- NA_real_
fig4_grid$rho_lower_loose <- NA_real_
fig4_grid$rho_lower_tight <- NA_real_
fig4_grid$rho_upper <- NA_real_
for (i in seq_len(nrow(fig4_grid))) {
  hazard_c_i <- log(2) / mst_for_fig4
  fig4_grid$rho_numerical[i] <- solve_rho_numerical(
    HR = fig4_grid$HR[i], hazard_c = hazard_c_i,
    f = f_for_fig4, d = fig4_grid$d[i],
    omega1 = omega1_for_fig4, r = fig4_grid$r[i]
  )
  bnd <- calc_rho_bounds(
    HR = fig4_grid$HR[i], hazard_c = hazard_c_i,
    d = fig4_grid$d[i], r = fig4_grid$r[i],
    time_unit = TIME_UNIT
  )
  fig4_grid$rho_lower_loose[i] <- bnd$bounds_df$rho_lower_loose
  fig4_grid$rho_lower_tight[i] <- bnd$bounds_df$rho_lower_tight
  fig4_grid$rho_upper[i] <- bnd$bounds_df$rho_upper
}

fig4_long <- dplyr::bind_rows(
  data.frame(fig4_grid[, c("HR", "r", "d")],
             quantity = "Universal lower",
             value = fig4_grid$rho_lower_loose,
             stringsAsFactors = FALSE),
  data.frame(fig4_grid[, c("HR", "r", "d")],
             quantity = "Tight lower",
             value = fig4_grid$rho_lower_tight,
             stringsAsFactors = FALSE),
  data.frame(fig4_grid[, c("HR", "r", "d")],
             quantity = "Numerical rho*",
             value = fig4_grid$rho_numerical,
             stringsAsFactors = FALSE),
  data.frame(fig4_grid[, c("HR", "r", "d")],
             quantity = "Upper",
             value = fig4_grid$rho_upper,
             stringsAsFactors = FALSE)
)
fig4_long$quantity <- factor(
  fig4_long$quantity,
  levels = c("Universal lower", "Tight lower", "Numerical rho*", "Upper")
)
fig4_long$d_label <- factor(
  make_d_label(fig4_long$d),
  levels = make_d_label(sort(unique(fig4_long$d)))
)
fig4_long$HR_label <- factor(
  paste0("HR == ", format(fig4_long$HR, nsmall = 1)),
  levels = paste0("HR == ", format(sort(unique(fig4_long$HR)), nsmall = 1))
)

quantity_labels_fig4 <- c(
  "Universal lower" = expression(rho[L]^"univ"),
  "Tight lower"     = expression(rho[L]^"tight"),
  "Numerical rho*"  = expression(rho * "*"),
  "Upper"           = expression(rho[U])
)

quantity_colors_fig4 <- c(
  "Universal lower" = "#7f7f7f",
  "Tight lower"     = "#1f77b4",
  "Numerical rho*"  = "#d62728",
  "Upper"           = "#2ca02c"
)
quantity_linetypes_fig4 <- c(
  "Universal lower" = "dashed",
  "Tight lower"     = "solid",
  "Numerical rho*"  = "solid",
  "Upper"           = "dotdash"
)
quantity_widths_fig4 <- c(
  "Universal lower" = 0.6,
  "Tight lower"     = 0.6,
  "Numerical rho*"  = 0.9,
  "Upper"           = 0.6
)

p_fig4 <- ggplot(fig4_long,
                 aes(x = r, y = value,
                     color = quantity,
                     linetype = quantity,
                     linewidth = quantity,
                     group = quantity)) +
  geom_line() +
  facet_grid(rows = vars(HR_label), cols = vars(d_label),
             labeller = label_parsed) +
  scale_color_manual(
    name = NULL,
    values = quantity_colors_fig4,
    labels = quantity_labels_fig4
  ) +
  scale_linetype_manual(
    name = NULL,
    values = quantity_linetypes_fig4,
    labels = quantity_labels_fig4
  ) +
  scale_linewidth_manual(
    name = NULL,
    values = quantity_widths_fig4,
    labels = quantity_labels_fig4
  ) +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) +
  labs(x = expression("Allocation ratio " * r),
       y = expression(rho)) +
  guides(
    color    = guide_legend(keywidth = unit(1.0, "cm")),
    linetype = guide_legend(keywidth = unit(1.0, "cm")),
    linewidth = guide_legend(keywidth = unit(1.0, "cm"))
  ) +
  theme_bw(base_size = BASE_SIZE) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.5, "cm"),
        strip.background = element_rect(fill = "grey92"))

ggsave(file.path(OUT_DIR, "fig4_bounds_vs_HR.eps"),
       plot = p_fig4,
       width = FIG_WIDTH_MM, height = FIG_HEIGHT_MM * 1.7,
       units = "mm", device = cairo_ps, fallback_resolution = 600)

# =============================================================================
# Table 1: Hypothetical trial examples
#
# Three illustrative time-to-event trial scenarios used in the application
# section: a cardiovascular outcome trial (HR = 0.8), a phase 3 oncology
# trial (HR = 0.7), and a randomized phase 2 oncology trial (HR = 0.6).
# Each example is reported under r = 2, r = 3, and r = 4 so that the
# implications of allocation ratio choice are visible.
# =============================================================================
hypothetical_examples <- data.frame(
  example = c("CV outcome trial",
              "Phase 3 oncology",
              "Phase 2 oncology"),
  abbrev = c("CV", "Onc-III", "Onc-II"),
  HR = c(0.8, 0.7, 0.6),
  mst_c = c(60, 18, 9),
  f = c(20, 12, 6),
  d = c(0.15, 0.10, 0.05),
  omega1 = c(100, 50, 25),
  alpha = c(0.025, 0.025, 0.10),
  power = c(0.9, 0.9, 0.8),
  stringsAsFactors = FALSE
)
hypothetical_r <- c(2, 3, 4)

solve_rho_for_example <- function(HR, hazard_c, f, d, omega1,
                                  alpha_one_sided, power_lvl, r) {
  res <- tryCatch(
    find_rho_for_target_re(
      HR = HR, hazard_c = hazard_c,
      a_k = 0, omega_k = omega1,
      f = f, d_t = d, d_c = d,
      r = r,
      target_re = 1, metric = "tau",
      alpha = alpha_one_sided, power = power_lvl,
      time_unit = TIME_UNIT, formula = "Schoenfeld"
    ),
    error = function(e) NULL
  )
  if (is.null(res)) return(NA_real_)
  if (!is.null(res$rho)) return(as.numeric(res$rho))
  if (!is.null(res$achieved_rho)) return(as.numeric(res$achieved_rho))
  if (is.numeric(res)) return(as.numeric(res))
  NA_real_
}

calc_tau_for_example <- function(HR, hazard_c, f, d, omega1,
                                 alpha_one_sided, power_lvl, r, rho) {
  td <- tryCatch(
    calc_trial_duration(
      HR = HR, hazard_c = hazard_c,
      a_k = 0, omega_k = omega1 * rho,
      f = f, d_t = d, d_c = d, r = r,
      alpha = alpha_one_sided, power = power_lvl,
      time_unit = TIME_UNIT, formula = "Schoenfeld"
    ),
    error = function(e) NULL
  )
  if (is.null(td)) {
    return(list(tau = NA_real_, E = NA_real_,
                n_t = NA_real_, n_c = NA_real_, N = NA_real_))
  }
  # Required events: ceiling so that the integer count meets the
  # continuous Lachin-Foulkes target.
  E_int <- if (!is.null(td$E)) ceiling(as.numeric(td$E)) else NA_real_
  # Sample sizes: round up the control-arm raw count to an integer, then
  # set experimental-arm = r * n_c so that the total is divisible by
  # (1 + r). This is the usual practical rounding for r:1 allocation.
  n_c_raw <- if (!is.null(td$n_c)) as.numeric(td$n_c) else NA_real_
  if (!is.na(n_c_raw)) {
    n_c_int <- ceiling(n_c_raw)
    n_t_int <- r * n_c_int
    N_int <- n_c_int + n_t_int
  } else {
    n_c_int <- NA_real_
    n_t_int <- NA_real_
    N_int <- NA_real_
  }
  list(
    tau = if (!is.null(td$tau)) as.numeric(td$tau) else NA_real_,
    E = E_int,
    n_t = n_t_int,
    n_c = n_c_int,
    N = N_int
  )
}

# Build a scenario-level table of bounds and rho* for each (example, r)
hypothetical_table <- expand.grid(
  example = hypothetical_examples$example,
  r = hypothetical_r,
  stringsAsFactors = FALSE
) |>
  dplyr::left_join(hypothetical_examples, by = "example") |>
  dplyr::arrange(.data$HR, .data$r)

hypothetical_table$hazard_c <- log(2) / hypothetical_table$mst_c
hypothetical_table$rho_lower_loose <- NA_real_
hypothetical_table$rho_lower_tight <- NA_real_
hypothetical_table$rho_numerical <- NA_real_
hypothetical_table$rho_upper <- NA_real_
hypothetical_table$tau_1 <- NA_real_
hypothetical_table$E_r <- NA_real_
hypothetical_table$n_t_r <- NA_real_
hypothetical_table$n_c_r <- NA_real_
hypothetical_table$N_r <- NA_real_

for (i in seq_len(nrow(hypothetical_table))) {
  row <- hypothetical_table[i, ]
  bnd <- calc_rho_bounds(
    HR = row$HR, hazard_c = row$hazard_c,
    d = row$d, r = row$r, time_unit = TIME_UNIT
  )
  hypothetical_table$rho_lower_loose[i] <- bnd$bounds_df$rho_lower_loose
  hypothetical_table$rho_lower_tight[i] <- bnd$bounds_df$rho_lower_tight
  hypothetical_table$rho_upper[i] <- bnd$bounds_df$rho_upper

  hypothetical_table$rho_numerical[i] <- solve_rho_for_example(
    HR = row$HR, hazard_c = row$hazard_c,
    f = row$f, d = row$d, omega1 = row$omega1,
    alpha_one_sided = row$alpha, power_lvl = row$power, r = row$r
  )

  baseline <- calc_tau_for_example(
    HR = row$HR, hazard_c = row$hazard_c,
    f = row$f, d = row$d, omega1 = row$omega1,
    alpha_one_sided = row$alpha, power_lvl = row$power,
    r = 1, rho = 1
  )
  hypothetical_table$tau_1[i] <- baseline$tau

  rstar <- calc_tau_for_example(
    HR = row$HR, hazard_c = row$hazard_c,
    f = row$f, d = row$d, omega1 = row$omega1,
    alpha_one_sided = row$alpha, power_lvl = row$power,
    r = row$r, rho = hypothetical_table$rho_numerical[i]
  )
  hypothetical_table$E_r[i] <- rstar$E
  hypothetical_table$n_t_r[i] <- rstar$n_t
  hypothetical_table$n_c_r[i] <- rstar$n_c
  hypothetical_table$N_r[i] <- rstar$N
}

utils::write.csv(
  hypothetical_table[, c("example", "HR", "mst_c", "f", "d", "omega1",
                         "alpha", "power", "r",
                         "rho_lower_loose", "rho_lower_tight",
                         "rho_numerical", "rho_upper",
                         "tau_1", "E_r", "n_t_r", "n_c_r", "N_r")],
  file = file.path(OUT_DIR, "table_hypothetical.csv"),
  row.names = FALSE
)

# Build a LaTeX table grouping by example
fmt3 <- function(x) {
  ifelse(is.na(x), "---", sprintf("%.3f", x))
}
fmt1 <- function(x) {
  ifelse(is.na(x), "---", sprintf("%.1f", x))
}

build_hypothetical_tex <- function(df, examples_meta) {
  header <- c(
    "\\begin{tabular}{lcccccccccccccccc}",
    "\\toprule",
    " & & & & & & & \\multicolumn{4}{c}{Bounds and $\\rho^{*}$} & & & & & \\\\",
    "\\cmidrule(lr){8-11}",
    "Example & HR & ${\\rm MST}_{C}$ & $f$ & $d$ & $\\omega_{1}$ & $r$ & $\\rho_{L}^{\\rm univ}$ & $\\rho_{L}^{\\rm tight}$ & $\\rho^{*}$ & $\\rho_{U}$ & $\\tau_{1}$ & $E_{r}$ & $n_{T,r}$ & $n_{C,r}$ & $N_{r}$ \\\\",
    "\\midrule"
  )
  rows <- character(0)
  ex_names <- unique(df$example)
  for (g in seq_along(ex_names)) {
    ex <- ex_names[g]
    sub <- df[df$example == ex, ]
    sub <- sub[order(sub$r), ]
    for (i in seq_len(nrow(sub))) {
      first_cols <- if (i == 1) {
        sprintf("%s & %.1f & %d & %d & %.2f & %d",
                sub$example[i], sub$HR[i],
                as.integer(sub$mst_c[i]), as.integer(sub$f[i]),
                sub$d[i], as.integer(sub$omega1[i]))
      } else {
        " & & & & & "
      }
      tau_col <- if (i == 1) fmt1(sub$tau_1[i]) else ""
      rows <- c(rows, sprintf(
        " %s & %d & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
        first_cols, as.integer(sub$r[i]),
        fmt3(sub$rho_lower_loose[i]),
        fmt3(sub$rho_lower_tight[i]),
        fmt3(sub$rho_numerical[i]),
        fmt3(sub$rho_upper[i]),
        tau_col,
        sprintf("%d", as.integer(sub$E_r[i])),
        sprintf("%d", as.integer(sub$n_t_r[i])),
        sprintf("%d", as.integer(sub$n_c_r[i])),
        sprintf("%d", as.integer(sub$N_r[i]))
      ))
    }
    if (g < length(ex_names)) {
      rows <- c(rows, "\\midrule")
    }
  }
  footer <- c(
    "\\bottomrule",
    "\\end{tabular}"
  )
  paste(c(header, rows, footer), collapse = "\n")
}

writeLines(build_hypothetical_tex(hypothetical_table, hypothetical_examples),
           con = file.path(OUT_DIR, "table_hypothetical.tex"))

# =============================================================================
# Done. Outputs:
#   fig1_re_phase2.eps            -> Phase 2 RE curves (HR = 0.6)
#   fig2_re_phase3.eps            -> Phase 3 RE curves (HR = 0.8)
#   fig3_univ_accuracy.eps        -> Signed relative gaps of all 3 bounds
#   fig4_bounds_vs_HR.eps         -> Four quantities vs r, faceted by HR x d
#   table_hypothetical.tex / .csv -> Application section table for three
#                                    hypothetical trial examples
# =============================================================================
