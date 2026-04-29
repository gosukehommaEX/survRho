#' Print Method for Objects of Class "rho_bounds_grid"
#'
#' Displays a compact summary of a grid-wide analytical bounds calculation
#' produced by \code{\link{calc_rho_bounds_grid}}, including the parameter
#' values explored, the total number of scenarios, and summary statistics
#' of the three bounds across the grid.
#'
#' @param x An object of class \code{"rho_bounds_grid"}.
#' @param digits Integer. Number of digits to use when rounding numeric
#'   summaries shown in the output. Default is 3.
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#'
#' @return The input object \code{x}, returned invisibly. Called primarily
#'   for its side effect of printing to the console.
#'
#' @details
#' Because a grid can contain many scenarios, this method does not print
#' the full \code{bounds_df} table. Instead, it summarises the parameter
#' ranges and the empirical distribution of rho_lower_loose,
#' rho_lower_tight, and rho_upper. For full data inspection, access
#' \code{x$bounds_df} directly.
#'
#' @examples
#' \dontrun{
#' results <- calc_rho_bounds_grid(
#'   HR = c(0.6, 0.7, 0.8),
#'   hazard_c = c(log(2) / 6, log(2) / 12),
#'   d = c(0, 0.3),
#'   r = c(2, 3, 4)
#' )
#' print(results)
#' }
#'
#' @seealso \code{\link{calc_rho_bounds_grid}}.
#'
#' @exportS3Method print rho_bounds_grid
print.rho_bounds_grid <- function(x, digits = 3, ...) {
  
  bounds_df <- x$bounds_df
  unit <- x$time_unit
  title <- "Analytical bounds on the required accrual rate multiplier (grid)"
  
  # --- Extract unique values for each design dimension ---
  HR_vals <- sort(unique(bounds_df$HR))
  mst_c_vals <- sort(unique(bounds_df$mst_c))
  d_vals <- sort(unique(bounds_df$d))
  r_vals <- sort(unique(bounds_df$r))
  
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
      all(abs(diffs - diffs[1]) < 10 ^ (-nd - 2))
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
  lines <- c(lines, sprintf("  Dropout rate (d):           %d value%s: %s",
                            length(d_vals),
                            if (length(d_vals) == 1) "" else "s",
                            fmt_vals(d_vals, 2)))
  lines <- c(lines, sprintf("  Allocation ratio (r):       %d value%s: %s",
                            length(r_vals),
                            if (length(r_vals) == 1) "" else "s",
                            fmt_vals(r_vals, 0)))
  lines <- c(lines, "")
  
  # --- Grid summary ---
  n_scenarios <- nrow(bounds_df)
  lines <- c(lines, "Grid summary")
  lines <- c(lines, sprintf("  Total (HR, hazard_c, d, r)  scenarios: %d",
                            n_scenarios))
  lines <- c(lines, "")
  
  # --- Bounds summary statistics ---
  lines <- c(lines, "Bounds summary (across all scenarios)")
  bound_summary <- rbind(
    rho_lower_loose = summary_stats_bounds(bounds_df$rho_lower_loose, digits),
    rho_lower_tight = summary_stats_bounds(bounds_df$rho_lower_tight, digits),
    rho_upper = summary_stats_bounds(bounds_df$rho_upper, digits)
  )
  bound_summary_df <- as.data.frame(bound_summary)
  tbl_lines <- utils::capture.output(print(bound_summary_df))
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
# Named distinctly from summary_stats() in print.re_grid.R to avoid
# clobbering when both files are sourced.
summary_stats_bounds <- function(v, digits) {
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
