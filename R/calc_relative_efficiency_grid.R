#' Calculate Relative Efficiency Over a Grid of Parameter Combinations
#'
#' Evaluates the relative efficiency of r:1 randomization compared to 1:1
#' randomization over all combinations of hazard ratios, control-group hazard
#' rates, follow-up periods, and common dropout rates. For each parameter
#' combination, \code{\link{calc_relative_efficiency}} is called with the
#' specified allocation ratios and accrual rate multipliers.
#'
#' @param HR A numeric vector of hazard ratios under the alternative hypothesis.
#' @param hazard_c A numeric vector of hazard rates for the control group
#'   (per time unit).
#' @param f A numeric vector of follow-up times after completion of accrual.
#' @param d A numeric vector of dropout rates (annual probability) applied to
#'   both the experimental and control groups (d_t = d_c = d).
#' @param a_k A numeric vector of accrual time interval boundaries.
#' @param omega_k A numeric vector of accrual intensities (patients per time
#'   unit) for each accrual interval.
#' @param r A numeric vector of allocation ratios to the experimental group
#'   (T:C = r:1).
#' @param rho A numeric vector of accrual rate multipliers, the same length
#'   as r.
#' @param alpha A one-sided level of significance.
#' @param power A target power.
#' @param time_unit A character string specifying the time unit. One of
#'   "months" (default), "weeks", or "years".
#' @param formula A character string specifying the formula used to calculate
#'   the required number of events. One of "Schoenfeld" (default) or "Freedman".
#' @param verbose Logical. If TRUE, progress messages are printed for grids
#'   with more than 50 combinations. Default is FALSE.
#'
#' @return An object of class \code{"re_grid"} inheriting from \code{"list"},
#'   with components:
#' \describe{
#'   \item{original_df}{Full output of calc_trial_duration() for every
#'     parameter combination, including the baseline (r = 1) case. Columns
#'     HR, hazard_c, f, and d are appended for identification.}
#'   \item{re_df}{Relative efficiency results with one row per (HR, hazard_c,
#'     f, d, r, rho) combination. Columns include HR, hazard_c, mst_c
#'     (= log(2) / hazard_c), f, d, r, rho, re_n, re_accrual_duration, and
#'     re_tau.}
#' }
#' When printed, a compact summary of the grid is shown; see
#' \code{\link{print.re_grid}}. The list components remain accessible in the
#' usual way.
#'
#' @details
#' The total number of scenarios evaluated is
#' length(HR) x length(hazard_c) x length(f) x length(d) x (1 + length(r)),
#' where "1 +" accounts for the baseline r = 1 case in each combination.
#'
#' The relative efficiency for outcome x is
#' \code{RE_x = value_{r = 1} / value_{r}},
#' so that RE_x > 1 indicates r:1 randomization is more efficient than 1:1.
#'
#' @note This function requires \code{\link{calc_relative_efficiency}} and
#'   the \code{dplyr} package.
#'
#' @examples
#' \dontrun{
#' results <- calc_relative_efficiency_grid(
#'   HR = seq(0.6, 0.9, by = 0.1),
#'   hazard_c = c(log(2) / 6, log(2) / 12, log(2) / 24),
#'   f = c(12, 18, 24),
#'   d = c(0, 0.2, 0.4),
#'   a_k = c(0, 3, 6),
#'   omega_k = c(15, 20, 25),
#'   r = c(2, 3),
#'   rho = c(1.2, 1.3),
#'   alpha = 0.025,
#'   power = 0.9,
#'   time_unit = "months",
#'   formula = "Schoenfeld",
#'   verbose = TRUE
#' )
#' print(results)
#' }
#'
#' @seealso \code{\link{calc_relative_efficiency}} for single parameter set
#'   analysis; \code{\link{print.re_grid}} for
#'   display.
#'
#' @importFrom dplyr bind_rows
#'
#' @export
calc_relative_efficiency_grid <- function(HR, hazard_c, f, d, a_k, omega_k,
                                          r, rho, alpha, power,
                                          time_unit = "months",
                                          formula = "Schoenfeld",
                                          verbose = FALSE) {
  
  # --- Validate input lengths ---
  if (length(r) != length(rho)) {
    stop("Vectors r and rho must have the same length")
  }
  
  # --- Build parameter grid ---
  param_grid <- expand.grid(
    HR = HR,
    hazard_c = hazard_c,
    f = f,
    d = d,
    stringsAsFactors = FALSE
  )
  n_combi <- nrow(param_grid)
  
  # --- Iterate over parameter combinations ---
  original_list <- vector("list", n_combi)
  re_list <- vector("list", n_combi)
  
  for (i in seq_len(n_combi)) {
    current_HR <- param_grid$HR[i]
    current_hazard_c <- param_grid$hazard_c[i]
    current_f <- param_grid$f[i]
    current_d <- param_grid$d[i]
    
    current_results <- calc_relative_efficiency(
      HR = current_HR,
      hazard_c = current_hazard_c,
      a_k = a_k,
      omega_k = omega_k,
      f = current_f,
      d_t = current_d,
      d_c = current_d,
      r = r,
      rho = rho,
      alpha = alpha,
      power = power,
      time_unit = time_unit,
      formula = formula
    )
    
    # Tag each row with the parameter combination identifiers
    current_original <- current_results$original_df
    current_original$HR <- current_HR
    current_original$hazard_c <- current_hazard_c
    current_original$f <- current_f
    current_original$d <- current_d
    
    current_re <- current_results$re_df
    current_re$HR <- current_HR
    current_re$hazard_c <- current_hazard_c
    current_re$mst_c <- log(2) / current_hazard_c
    current_re$f <- current_f
    current_re$d <- current_d
    
    original_list[[i]] <- current_original
    re_list[[i]] <- current_re
    
    if (verbose && n_combi > 50 && i %% 10 == 0) {
      message("Processed ", i, " of ", n_combi, " combinations")
    }
  }
  
  # --- Combine results ---
  original_df <- dplyr::bind_rows(original_list)
  re_df <- dplyr::bind_rows(re_list)
  
  # --- Reorder columns: identifiers first ---
  if (nrow(original_df) > 0) {
    id_cols <- c("HR", "hazard_c", "f", "d")
    other_cols <- setdiff(names(original_df), id_cols)
    original_df <- original_df[, c(id_cols, other_cols)]
    rownames(original_df) <- NULL
  }
  
  if (nrow(re_df) > 0) {
    id_cols <- c("HR", "hazard_c", "mst_c", "f", "d", "r", "rho")
    other_cols <- setdiff(names(re_df), id_cols)
    re_df <- re_df[, c(id_cols, other_cols)]
    rownames(re_df) <- NULL
  }
  
  if (verbose) {
    message("Completed ", n_combi, " parameter combinations.")
  }
  
  out <- list(original_df = original_df, re_df = re_df)
  class(out) <- c("re_grid", "list")
  return(out)
}