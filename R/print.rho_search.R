#' Print Method for Objects of Class "rho_search"
#'
#' Displays a human-readable summary of a root search for the accrual rate
#' multiplier rho that achieves a target relative efficiency, as produced by
#' \code{\link{find_rho_for_target_re}}. The output includes the search
#' specification (metric, target RE, allocation ratio), the result
#' (found rho, achieved RE, objective, iterations), and a final status
#' indicator.
#'
#' @param x An object of class \code{"rho_search"}.
#' @param digits Integer. Number of digits to use when rounding numeric
#'   summaries shown in the output. Default is 4.
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#'
#' @return The input object \code{x}, returned invisibly. Called primarily
#'   for its side effect of printing to the console.
#'
#' @details
#' The status line at the bottom is either "SUCCESS" (converged with small
#' objective), "FAILED" (did not converge), or "NOT AVAILABLE" (the root
#' finder returned NULL).
#'
#' @examples
#' \dontrun{
#' result <- find_rho_for_target_re(
#'   HR = 0.7, hazard_c = log(2) / 12,
#'   a_k = c(0, 3, 6), omega_k = c(15, 20, 25),
#'   f = 12, d_t = 0.3, d_c = 0.3,
#'   r = 2, alpha = 0.025, power = 0.9,
#'   metric = "tau", target_re = 1
#' )
#' print(result)
#' }
#'
#' @seealso \code{\link{find_rho_for_target_re}}.
#'
#' @exportS3Method print rho_search
print.rho_search <- function(x, digits = 4, ...) {
  
  title <- "Root search for target relative efficiency"
  
  # Human-readable description of the metric
  metric_desc <- switch(x$metric,
                        "n" = "n  (total sample size)",
                        "accrual_duration" = "accrual_duration  (a_K)",
                        "tau" = "tau  (total trial duration)",
                        x$metric)
  
  # Collect output lines
  lines <- character(0)
  
  # --- Search specification ---
  lines <- c(lines, "Search specification")
  lines <- c(lines, sprintf("  Outcome (metric):         %s", metric_desc))
  lines <- c(lines, sprintf("  Target RE:                %s",
                            format(round(x$target_re, digits),
                                   nsmall = digits)))
  lines <- c(lines, sprintf("  Allocation ratio (r):     %s",
                            format(round(x$r, digits))))
  lines <- c(lines, "")
  
  # --- Result ---
  lines <- c(lines, "Result")
  if (is.na(x$rho)) {
    lines <- c(lines, "  Rho found:                (not available)")
    lines <- c(lines, sprintf("  Converged:                %s", x$converged))
  } else {
    lines <- c(lines, sprintf("  Rho found:                %s",
                              format(round(x$rho, digits), nsmall = digits)))
    lines <- c(lines, sprintf("  Achieved RE:              %s",
                              format(round(x$re, digits), nsmall = digits)))
    if (!is.na(x$objective)) {
      lines <- c(lines, sprintf("  Objective |re - target|:  %s",
                                format(x$objective, scientific = TRUE,
                                       digits = 3)))
    }
    if (!is.na(x$iterations)) {
      lines <- c(lines, sprintf("  Iterations:               %d",
                                as.integer(x$iterations)))
    }
    lines <- c(lines, sprintf("  Converged:                %s", x$converged))
  }
  
  # --- Status ---
  status <- if (is.na(x$rho)) {
    "NOT AVAILABLE"
  } else if (isTRUE(x$converged)) {
    "SUCCESS"
  } else {
    "FAILED"
  }
  footer <- sprintf("Status: %s", status)
  
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