#' Compute Analytical Bounds on rho Over a Grid of Parameter Combinations
#'
#' Evaluates analytical bounds on the required accrual rate multiplier rho
#' over all combinations of hazard ratios, control-group hazard rates, and
#' dropout rates, for a set of allocation ratios r. For each parameter
#' combination, \code{\link{calc_rho_bounds}} is called and the results are
#' stacked into a single data frame.
#'
#' @param HR A numeric vector of hazard ratios under the alternative
#'   hypothesis. All values must satisfy \code{0 < HR < 1}.
#' @param hazard_c A numeric vector of hazard rates for the control group
#'   (per time unit).
#' @param d A numeric vector of dropout rates (annual probability) common to
#'   both groups. The underlying theory assumes a single dropout hazard
#'   shared across groups.
#' @param r A numeric vector of allocation ratios to the experimental group
#'   (T:C = r:1). All values must satisfy \code{r > 1}.
#' @param time_unit A character string specifying the time unit. One of
#'   "months" (default), "weeks", or "years".
#' @param verbose Logical. If TRUE, progress messages are printed for grids
#'   with more than 50 combinations. Default is FALSE.
#'
#' @return An object of class \code{"rho_bounds_grid"} inheriting from
#'   \code{"list"}, with components:
#' \describe{
#'   \item{bounds_df}{Bounds results with one row per (HR, hazard_c, d, r)
#'     combination. Columns include HR, hazard_c, mst_c (= log(2) / hazard_c),
#'     d, gamma, r, rho_lower_loose, rho_lower_tight, rho_upper, and
#'     theta_U.}
#'   \item{time_unit}{The input time unit (echoed).}
#' }
#' When printed, a compact summary of the grid is shown; see
#' \code{\link{print.rho_bounds_grid}}.
#'
#' @details
#' The total number of scenarios evaluated is
#' length(HR) x length(hazard_c) x length(d) x length(r).
#'
#' @note This function requires \code{\link{calc_rho_bounds}} and the
#'   \code{dplyr} package.
#'
#' @examples
#' \dontrun{
#' results <- calc_rho_bounds_grid(
#'   HR = c(0.6, 0.7, 0.8),
#'   hazard_c = c(log(2) / 6, log(2) / 12, log(2) / 24),
#'   d = c(0, 0.3),
#'   r = c(2, 3, 4),
#'   time_unit = "months",
#'   verbose = TRUE
#' )
#' print(results)
#' }
#'
#' @seealso \code{\link{calc_rho_bounds}} for single parameter set
#'   analysis; \code{\link{print.rho_bounds_grid}} for display.
#'
#' @importFrom dplyr bind_rows
#'
#' @export
calc_rho_bounds_grid <- function(HR, hazard_c, d, r,
                                 time_unit = "months",
                                 verbose = FALSE) {

  # --- Validate inputs (vector versions; per-element checks happen inside
  #     calc_rho_bounds via each call) ---
  if (!is.numeric(HR) || any(HR <= 0) || any(HR >= 1)) {
    stop("All values in HR must satisfy 0 < HR < 1")
  }
  if (!is.numeric(hazard_c) || any(hazard_c <= 0)) {
    stop("All values in hazard_c must be positive")
  }
  if (!is.numeric(d) || any(d < 0) || any(d >= 1)) {
    stop("All values in d must be in [0, 1)")
  }
  if (!is.numeric(r) || any(r <= 1)) {
    stop("All values in r must be greater than 1")
  }

  # --- Build parameter grid ---
  param_grid <- expand.grid(
    HR = HR,
    hazard_c = hazard_c,
    d = d,
    stringsAsFactors = FALSE
  )
  n_combi <- nrow(param_grid)

  # --- Iterate over parameter combinations ---
  bounds_list <- vector("list", n_combi)

  for (i in seq_len(n_combi)) {
    current_HR <- param_grid$HR[i]
    current_hazard_c <- param_grid$hazard_c[i]
    current_d <- param_grid$d[i]

    current_result <- calc_rho_bounds(
      HR = current_HR,
      hazard_c = current_hazard_c,
      d = current_d,
      r = r,
      time_unit = time_unit
    )

    # Tag each row with the parameter combination identifiers
    current_df <- current_result$bounds_df
    current_df$HR <- current_HR
    current_df$hazard_c <- current_hazard_c
    current_df$mst_c <- log(2) / current_hazard_c
    current_df$d <- current_d
    current_df$gamma <- current_result$gamma

    bounds_list[[i]] <- current_df

    if (verbose && n_combi > 50 && i %% 10 == 0) {
      message("Processed ", i, " of ", n_combi, " combinations")
    }
  }

  # --- Combine results ---
  bounds_df <- dplyr::bind_rows(bounds_list)

  # --- Reorder columns: identifiers first ---
  if (nrow(bounds_df) > 0) {
    id_cols <- c("HR", "hazard_c", "mst_c", "d", "gamma", "r")
    bound_cols <- c("rho_lower_loose", "rho_lower_tight", "rho_upper",
                    "theta_U")
    bounds_df <- bounds_df[, c(id_cols, bound_cols)]
    rownames(bounds_df) <- NULL
  }

  if (verbose) {
    message("Completed ", n_combi, " parameter combinations.")
  }

  out <- list(bounds_df = bounds_df, time_unit = time_unit)
  class(out) <- c("rho_bounds_grid", "list")
  return(out)
}
