#' Print Method for Objects of Class "relative_efficiency"
#'
#' Displays a human-readable summary of a relative efficiency calculation
#' produced by \code{\link{calc_relative_efficiency}}, including the common
#' design parameters, the baseline (r = 1) trial characteristics, and a
#' compact table of trial characteristics and relative efficiency for each
#' (r, rho) combination.
#'
#' @param x An object of class \code{"relative_efficiency"}.
#' @param digits Integer. Number of digits to use when rounding numeric
#'   summaries shown in the output. Default is 3.
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#'
#' @return The input object \code{x}, returned invisibly. Called primarily
#'   for its side effect of printing to the console.
#'
#' @details
#' The horizontal separator lines are sized to match the longest content
#' line, so the width of the output adapts to the numeric values displayed.
#'
#' @examples
#' \dontrun{
#' result <- calc_relative_efficiency(
#'   HR = 0.7, hazard_c = log(2) / 12,
#'   a_k = c(0, 3, 6), omega_k = c(15, 20, 25),
#'   f = 12, d_t = 0.3, d_c = 0.3,
#'   r = c(2, 3, 4), rho = c(1.2, 1.3, 1.4),
#'   alpha = 0.025, power = 0.9
#' )
#' print(result)
#' }
#'
#' @seealso \code{\link{calc_relative_efficiency}}.
#'
#' @exportS3Method print relative_efficiency
print.relative_efficiency <- function(x, digits = 3, ...) {
  
  original_df <- x$original_df
  re_df <- x$re_df
  
  # Baseline row (r = 1, rho = 1)
  baseline_row <- original_df[original_df$r == 1 & original_df$rho == 1, ][1, ]
  
  unit <- baseline_row$time_unit
  title <- "Relative efficiency of unequal randomization"
  
  # Collect output lines
  lines <- character(0)
  
  # --- Common design parameters ---
  lines <- c(lines, "Common design parameters")
  lines <- c(lines, sprintf("  Hazard ratio (T vs C):    %s",
                            format(round(baseline_row$HR, digits),
                                   nsmall = digits)))
  lines <- c(lines, sprintf("  Control median survival:  %s %s",
                            format(round(baseline_row$mst_c, digits),
                                   nsmall = digits),
                            unit))
  lines <- c(lines, sprintf("  Follow-up period:         %s %s",
                            format(round(baseline_row$f, digits),
                                   nsmall = digits),
                            unit))
  lines <- c(lines, sprintf("  Dropout rate (T, C):      %s%%, %s%% per year",
                            format(round(100 * baseline_row$d_t, 1)),
                            format(round(100 * baseline_row$d_c, 1))))
  lines <- c(lines, sprintf("  Significance / Power:     %s / %s%%",
                            format(baseline_row$alpha, nsmall = 3),
                            format(round(100 * baseline_row$power, 1))))
  lines <- c(lines, sprintf("  Event formula:            %s",
                            baseline_row$formula))
  lines <- c(lines, "")
  
  # --- Baseline trial characteristics ---
  lines <- c(lines, "Baseline (r = 1, rho = 1)")
  lines <- c(lines, sprintf("  Sample size (N):          %s",
                            format(round(baseline_row$N, 1), nsmall = 1)))
  lines <- c(lines, sprintf("  Accrual duration (a_K):   %s %s",
                            format(round(baseline_row$a_K, digits),
                                   nsmall = digits),
                            unit))
  lines <- c(lines, sprintf("  Total trial duration:     %s %s",
                            format(round(baseline_row$tau, digits),
                                   nsmall = digits),
                            unit))
  lines <- c(lines, "")
  
  # --- RE comparison table (one row per r:1 scenario) ---
  # Merge per-r trial-level outputs (N, a_K, tau) with the RE table.
  non_baseline <- original_df[!(original_df$r == 1 & original_df$rho == 1), ]
  # Keep only identifier + main outputs to avoid width blow-up
  merged_tbl <- merge(
    re_df,
    non_baseline[, c("r", "rho", "N", "a_K", "tau")],
    by = c("r", "rho"),
    sort = FALSE
  )
  # Order columns for display
  display_tbl <- data.frame(
    r = merged_tbl$r,
    rho = round(merged_tbl$rho, digits),
    N = round(merged_tbl$N, 1),
    a_K = round(merged_tbl$a_K, digits),
    tau = round(merged_tbl$tau, digits),
    re_n = round(merged_tbl$re_n, digits),
    re_a_K = round(merged_tbl$re_accrual_duration, digits),
    re_tau = round(merged_tbl$re_tau, digits),
    stringsAsFactors = FALSE
  )
  lines <- c(lines, "Relative efficiency versus 1:1 randomization")
  tbl_lines <- utils::capture.output(print(display_tbl, row.names = FALSE))
  lines <- c(lines, paste0("  ", tbl_lines))
  
  # --- Footer ---
  footer <- sprintf("Time unit: %s", unit)
  
  # --- Determine separator width from the longest printed line ---
  all_for_width <- c(title, lines, footer)
  sep_width <- max(nchar(all_for_width, type = "width"))
  sep_line <- strrep("-", sep_width)
  
  # --- Emit output ---
  cat(title, "\n", sep = "")
  cat(sep_line, "\n", sep = "")
  cat(paste(lines, collapse = "\n"), "\n", sep = "")
  cat(sep_line, "\n", sep = "")
  cat(footer, "\n", sep = "")
  
  invisible(x)
}