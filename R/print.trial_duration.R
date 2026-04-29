#' Print Method for Objects of Class "trial_duration"
#'
#' Displays a human-readable summary of a trial design calculation produced
#' by \code{\link{calc_trial_duration}}, including treatment effect, control
#' survival, allocation ratio, test operating characteristics, and the
#' resulting required events, sample size, accrual duration, and total trial
#' duration.
#'
#' @param x An object of class \code{"trial_duration"}.
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
#' The accrual intensities (omega_k) and interval boundaries (a_k) are
#' displayed in a compact table in the accrual section; see
#' \code{build_accrual_table}.
#'
#' @examples
#' \dontrun{
#' result <- calc_trial_duration(
#'   HR = 0.7, hazard_c = log(2) / 12,
#'   a_k = c(0, 3, 6), omega_k = c(15, 20, 25),
#'   f = 12, d_t = 0.3, d_c = 0.3, r = 2,
#'   alpha = 0.025, power = 0.9
#' )
#' print(result)
#' print(result, digits = 4)
#' }
#'
#' @seealso \code{\link{calc_trial_duration}}.
#'
#' @exportS3Method print trial_duration
print.trial_duration <- function(x, digits = 3, ...) {
  
  unit <- x$time_unit
  title <- "Trial design calculation (exponential survival)"
  
  # Collect output lines into a character vector so that separator widths
  # can be adapted to the longest line actually printed.
  lines <- character(0)
  
  # --- Treatment effect section ---
  lines <- c(lines, "Treatment effect")
  lines <- c(lines, sprintf("  Hazard ratio (T vs C):    %s (log HR = %s)",
                            format(round(x$HR, digits), nsmall = digits),
                            format(round(log(x$HR), digits), nsmall = digits)))
  lines <- c(lines, sprintf("  Median survival (T, C):   %s, %s %s",
                            format(round(x$mst_t, digits), nsmall = digits),
                            format(round(x$mst_c, digits), nsmall = digits),
                            unit))
  lines <- c(lines, sprintf("  Hazard rate (T, C):       %s, %s per %s",
                            format(round(x$hazard_t, digits + 2),
                                   nsmall = digits + 2),
                            format(round(x$hazard_c, digits + 2),
                                   nsmall = digits + 2),
                            sub("s$", "", unit)))
  lines <- c(lines, "")
  
  # --- Design section ---
  lines <- c(lines, "Design")
  lines <- c(lines, sprintf("  Allocation ratio (T:C):   %s:1",
                            format(round(x$r, digits))))
  lines <- c(lines, sprintf("  Significance (one-sided): alpha = %s",
                            format(x$alpha, nsmall = 3)))
  lines <- c(lines, sprintf("  Target power:             %s%%",
                            format(round(100 * x$power, 1))))
  lines <- c(lines, sprintf("  Dropout rate (T, C):      %s%%, %s%% per year",
                            format(round(100 * x$d_t, 1)),
                            format(round(100 * x$d_c, 1))))
  lines <- c(lines, sprintf("  Follow-up period:         %s %s",
                            format(round(x$f, digits), nsmall = digits),
                            unit))
  lines <- c(lines, "")
  
  # --- Accrual pattern section ---
  omega_cols <- grep("^omega_[0-9]+$", names(x), value = TRUE)
  a_cols <- grep("^a_[0-9]+$", names(x), value = TRUE)
  p_cols <- grep("^p_[0-9]+$", names(x), value = TRUE)
  
  accrual_lines <- NULL
  if (length(omega_cols) > 0 && length(a_cols) > 0) {
    interval_tbl <- build_accrual_table(x, omega_cols, a_cols, p_cols, digits)
    if (!is.null(interval_tbl)) {
      tbl_lines <- utils::capture.output(print(interval_tbl, row.names = FALSE))
      accrual_lines <- paste0("  ", tbl_lines)
    }
  }
  if (!is.null(accrual_lines)) {
    lines <- c(lines, "Accrual pattern", accrual_lines, "")
  }
  
  # --- Primary results section ---
  lines <- c(lines, "Primary results")
  lines <- c(lines, sprintf("  Required events (E):      %s",
                            format(round(x$E, 1), nsmall = 1)))
  lines <- c(lines, sprintf("  Required sample size:     %s (T: %s, C: %s)",
                            format(round(x$N, 1), nsmall = 1),
                            format(round(x$n_t, 1), nsmall = 1),
                            format(round(x$n_c, 1), nsmall = 1)))
  lines <- c(lines, sprintf("  Accrual duration (a_K):   %s %s",
                            format(round(x$a_K, digits), nsmall = digits),
                            unit))
  lines <- c(lines, sprintf("  Total trial duration:     %s %s (a_K + f)",
                            format(round(x$tau, digits), nsmall = digits),
                            unit))
  lines <- c(lines, "")
  
  # --- Footer line ---
  footer <- sprintf("Event formula: %s  |  Time unit: %s", x$formula, unit)
  
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