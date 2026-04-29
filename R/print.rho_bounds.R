#' Print Method for Objects of Class "rho_bounds"
#'
#' Displays a human-readable summary of analytical bounds on the required
#' accrual rate multiplier rho produced by \code{\link{calc_rho_bounds}},
#' including the common design parameters and a compact table of bounds
#' for each allocation ratio r.
#'
#' @param x An object of class \code{"rho_bounds"}.
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
#' result <- calc_rho_bounds(
#'   HR = 0.8, hazard_c = log(2) / 12,
#'   d = 0.3, r = c(2, 3, 4)
#' )
#' print(result)
#' }
#'
#' @seealso \code{\link{calc_rho_bounds}}.
#'
#' @exportS3Method print rho_bounds
print.rho_bounds <- function(x, digits = 3, ...) {
  
  bounds_df <- x$bounds_df
  unit <- x$time_unit
  title <- "Analytical bounds on the required accrual rate multiplier"
  
  # --- Derived values for display ---
  mst_c <- log(2) / x$hazard_c
  
  # Collect output lines
  lines <- character(0)
  
  # --- Common design parameters ---
  lines <- c(lines, "Common design parameters")
  lines <- c(lines, sprintf("  Hazard ratio (T vs C):    %s",
                            format(round(x$HR, digits), nsmall = digits)))
  lines <- c(lines, sprintf("  Control median survival:  %s %s",
                            format(round(mst_c, digits), nsmall = digits),
                            unit))
  lines <- c(lines, sprintf("  Dropout rate:             %s%% per year",
                            format(round(100 * x$d, 1))))
  lines <- c(lines, sprintf("  Dropout hazard (gamma):   %s per %s",
                            format(round(x$gamma, digits + 2),
                                   nsmall = digits + 2),
                            sub("s$", "", unit)))
  lines <- c(lines, sprintf("  Upper bound on theta:     %s",
                            format(round(bounds_df$theta_U[1], digits),
                                   nsmall = digits)))
  lines <- c(lines, "")
  
  # --- Bounds table (one row per r) ---
  display_tbl <- data.frame(
    r = bounds_df$r,
    rho_lower_loose = round(bounds_df$rho_lower_loose, digits),
    rho_lower_tight = round(bounds_df$rho_lower_tight, digits),
    rho_upper = round(bounds_df$rho_upper, digits),
    stringsAsFactors = FALSE
  )
  lines <- c(lines, "Bounds on rho required for RE = 1")
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