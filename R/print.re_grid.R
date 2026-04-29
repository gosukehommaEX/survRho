#' Print Method for Objects of Class "re_grid"
#'
#' Displays a compact summary of a grid-wide relative efficiency calculation
#' produced by \code{\link{calc_relative_efficiency_grid}}, including the
#' parameter values explored, the total number of scenarios, and summary
#' statistics of the relative efficiency across the grid.
#'
#' @param x An object of class \code{"re_grid"}.
#' @param digits Integer. Number of digits to use when rounding numeric
#'   summaries shown in the output. Default is 3.
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#'
#' @return The input object \code{x}, returned invisibly. Called primarily
#'   for its side effect of printing to the console.
#'
#' @details
#' Because a grid can contain hundreds of scenarios, this method does not
#' print the full \code{re_df} table. Instead, it summarises the parameter
#' ranges and the empirical distribution of re_n, re_accrual_duration, and
#' re_tau. For full data inspection, access \code{x$re_df} directly.
#'
#' @examples
#' \dontrun{
#' results <- calc_relative_efficiency_grid(
#'   HR = c(0.7, 0.8), hazard_c = c(log(2) / 6, log(2) / 12),
#'   f = c(12, 24), d = c(0, 0.3),
#'   a_k = c(0, 3, 6), omega_k = c(15, 20, 25),
#'   r = c(2, 3), rho = c(1.2, 1.3),
#'   alpha = 0.025, power = 0.9
#' )
#' print(results)
#' }
#'
#' @seealso \code{\link{calc_relative_efficiency_grid}}.
#'
#' @exportS3Method print re_grid
print.re_grid <- function(x, digits = 3, ...) {
  
  re_df <- x$re_df
  original_df <- x$original_df
  
  # --- Extract unique values for each design dimension ---
  HR_vals <- sort(unique(re_df$HR))
  mst_c_vals <- sort(unique(re_df$mst_c))
  f_vals <- sort(unique(re_df$f))
  d_vals <- sort(unique(re_df$d))
  r_vals <- sort(unique(re_df$r))
  rho_vals <- sort(unique(re_df$rho))
  
  # --- Common-across-grid values (take from the first row of original_df) ---
  first_row <- original_df[1, ]
  unit <- first_row$time_unit
  title <- "Relative efficiency over a parameter grid"
  
  # --- Helper to format a vector of values compactly ---
  # When the number of values exceeds max_inline, show a compact summary
  # (range, count, and step if evenly spaced) instead of listing them all.
  fmt_vals <- function(v, nd = digits, max_inline = 6) {
    rounded <- round(v, nd)
    if (length(rounded) <= max_inline) {
      return(paste(format(rounded, nsmall = nd), collapse = ", "))
    }
    sorted <- sort(rounded)
    diffs <- diff(sorted)
    is_even <- length(diffs) > 0 &&
      all(abs(diffs - diffs[1]) < 10^(-nd - 2))
    range_str <- sprintf("[%s, %s]",
                         format(sorted[1], nsmall = nd),
                         format(sorted[length(sorted)], nsmall = nd))
    if (is_even) {
      sprintf("%s (step = %s)", range_str,
              format(round(diffs[1], nd + 2), nsmall = nd))
    } else {
      range_str
    }
  }
  
  # Collect output lines
  lines <- character(0)
  
  # --- Design section ---
  lines <- c(lines, "Design")
  lines <- c(lines, sprintf("  Hazard ratio (HR):          %d value%s: %s",
                            length(HR_vals),
                            if (length(HR_vals) == 1) "" else "s",
                            fmt_vals(HR_vals)))
  lines <- c(lines, sprintf("  Control median survival:    %d value%s: %s %s",
                            length(mst_c_vals),
                            if (length(mst_c_vals) == 1) "" else "s",
                            fmt_vals(mst_c_vals, 1),
                            unit))
  lines <- c(lines, sprintf("  Follow-up period (f):       %d value%s: %s %s",
                            length(f_vals),
                            if (length(f_vals) == 1) "" else "s",
                            fmt_vals(f_vals, 1),
                            unit))
  lines <- c(lines, sprintf("  Dropout rate (d):           %d value%s: %s",
                            length(d_vals),
                            if (length(d_vals) == 1) "" else "s",
                            fmt_vals(d_vals, 2)))
  lines <- c(lines, sprintf("  Allocation ratio (r):       %d value%s: %s",
                            length(r_vals),
                            if (length(r_vals) == 1) "" else "s",
                            fmt_vals(r_vals, 0)))
  lines <- c(lines, sprintf("  Accrual multiplier (rho):   %d value%s: %s",
                            length(rho_vals),
                            if (length(rho_vals) == 1) "" else "s",
                            fmt_vals(rho_vals, 2)))
  lines <- c(lines, sprintf("  Significance / Power:       %s / %s%%",
                            format(first_row$alpha, nsmall = 3),
                            format(round(100 * first_row$power, 1))))
  lines <- c(lines, sprintf("  Event formula:              %s",
                            first_row$formula))
  lines <- c(lines, "")
  
  # --- Grid summary ---
  n_scenarios <- nrow(re_df)
  n_baseline <- sum(original_df$r == 1 & original_df$rho == 1)
  lines <- c(lines, "Grid summary")
  lines <- c(lines, sprintf("  Total r:1 scenarios:        %d", n_scenarios))
  lines <- c(lines, sprintf("  Baseline (r = 1) cases:     %d", n_baseline))
  lines <- c(lines, "")
  
  # --- RE summary statistics ---
  lines <- c(lines, "Relative efficiency summary (across all r:1 scenarios)")
  re_summary <- rbind(
    re_n = summary_stats(re_df$re_n, digits),
    re_a_K = summary_stats(re_df$re_accrual_duration, digits),
    re_tau = summary_stats(re_df$re_tau, digits)
  )
  re_summary_df <- as.data.frame(re_summary)
  tbl_lines <- utils::capture.output(print(re_summary_df))
  lines <- c(lines, paste0("  ", tbl_lines))
  
  # --- Footer ---
  footer2 <- sprintf("Time unit: %s", unit)
  
  # --- Determine separator width from the longest printed line ---
  all_for_width <- c(title, lines, footer2)
  sep_width <- max(nchar(all_for_width, type = "width"))
  sep_line <- strrep("-", sep_width)
  
  # --- Emit output ---
  cat(title, "\n", sep = "")
  cat(sep_line, "\n", sep = "")
  cat(paste(lines, collapse = "\n"), "\n", sep = "")
  cat(sep_line, "\n", sep = "")
  cat(footer2, "\n", sep = "")
  
  invisible(x)
}


# --- Internal helper: five-number summary plus mean, all rounded ---
summary_stats <- function(v, digits) {
  v <- v[is.finite(v)]
  if (length(v) == 0) {
    return(c(Min = NA, Q1 = NA, Median = NA, Mean = NA, Q3 = NA, Max = NA))
  }
  q <- stats::quantile(v, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  c(
    Min = round(q[[1]], digits),
    Q1 = round(q[[2]], digits),
    Median = round(q[[3]], digits),
    Mean = round(mean(v, na.rm = TRUE), digits),
    Q3 = round(q[[4]], digits),
    Max = round(q[[5]], digits)
  )
}