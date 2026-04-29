#' Find Rho for Target Relative Efficiency Over a Grid of Parameter Combinations
#'
#' Searches for the accrual rate multiplier rho that achieves a specified
#' target relative efficiency for each combination of trial design parameters
#' (hazard ratio, control-group hazard rate, follow-up period, dropout rate,
#' allocation ratio, and target RE). For each combination,
#' \code{\link{find_rho_for_target_re}} is called internally.
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
#' @param target_re A numeric vector of target relative efficiency values.
#' @param metric A character string specifying the outcome for which the
#'   target RE is sought. One of "n" (total sample size),
#'   "accrual_duration" (accrual duration a_K), or "tau" (total trial
#'   duration). Default is "n".
#' @param alpha A one-sided level of significance.
#' @param power A target power.
#' @param time_unit A character string specifying the time unit. One of
#'   "months" (default), "weeks", or "years".
#' @param rho_range A numeric vector of length 2 giving the search range for
#'   rho. Default is c(0.1, 5).
#' @param tolerance A numeric value specifying the tolerance for the root
#'   finder. Default is 1e-6.
#' @param max_iter A maximum number of iterations for the root finder.
#'   Default is 1000.
#' @param formula A character string specifying the formula used to calculate
#'   the required number of events. One of "Schoenfeld" (default) or "Freedman".
#' @param verbose Logical. If TRUE, progress messages are printed for grids
#'   with more than 50 combinations, and a summary is printed at the end.
#'   Default is FALSE.
#'
#' @return A data frame with one row per parameter combination, with columns:
#' \describe{
#'   \item{HR}{Hazard ratio}
#'   \item{hazard_c}{Control-group hazard rate}
#'   \item{mst_c}{Median survival time for the control group
#'     (= log(2) / hazard_c)}
#'   \item{f}{Follow-up time}
#'   \item{d}{Dropout rate (annual probability)}
#'   \item{r}{Allocation ratio (experimental:control)}
#'   \item{target_re}{Target relative efficiency}
#'   \item{metric}{Outcome used for optimization}
#'   \item{rho}{Rho value that achieves the target RE}
#'   \item{achieved_re}{Relative efficiency actually achieved at rho}
#'   \item{converged}{Logical indicating whether the root finder converged}
#'   \item{objective}{Absolute difference |achieved_re - target_re|}
#'   \item{iterations}{Number of iterations used by the root finder}
#'   \item{status}{"Success", "Failed", or "Error"}
#'   \item{error_message}{Error message if status = "Error", otherwise NA}
#' }
#'
#' @details
#' The total number of scenarios evaluated is
#' length(HR) x length(hazard_c) x length(f) x length(d) x length(r) x length(target_re).
#' Failures in individual combinations (e.g., target_re outside the achievable
#' range) are caught and reported in the status column; the remaining
#' combinations are still processed.
#'
#' @note This function requires \code{\link{find_rho_for_target_re}}.
#'
#' @examples
#' \dontrun{
#' results <- find_rho_for_target_re_grid(
#'   HR = c(0.6, 0.7, 0.8),
#'   hazard_c = c(log(2) / 6, log(2) / 12, log(2) / 24),
#'   f = c(12, 18, 24),
#'   d = c(0, 0.2, 0.4),
#'   a_k = c(0, 3, 6),
#'   omega_k = c(15, 20, 25),
#'   r = c(2, 3),
#'   target_re = c(0.9, 1, 1.1),
#'   metric = "n",
#'   alpha = 0.025,
#'   power = 0.9,
#'   time_unit = "months",
#'   formula = "Schoenfeld",
#'   verbose = TRUE
#' )
#' head(subset(results, status == "Success"))
#' }
#'
#' @seealso \code{\link{find_rho_for_target_re}} for single parameter set
#'   optimization, \code{\link{calc_relative_efficiency_grid}} for grid-wide
#'   relative efficiency analysis.
#'
#' @export
find_rho_for_target_re_grid <- function(HR, hazard_c, f, d, a_k, omega_k, r,
                                        target_re, metric = "n",
                                        alpha, power, time_unit = "months",
                                        rho_range = c(0.1, 5),
                                        tolerance = 1e-6, max_iter = 1000,
                                        formula = "Schoenfeld",
                                        verbose = FALSE) {

  # --- Validate inputs ---
  valid_metrics <- c("n", "accrual_duration", "tau")
  if (!metric %in% valid_metrics) {
    stop("metric must be one of: ", paste(valid_metrics, collapse = ", "))
  }

  if (!is.numeric(target_re) || any(target_re <= 0)) {
    stop("All target_re values must be positive")
  }

  # --- Build parameter grid ---
  param_grid <- expand.grid(
    HR = HR,
    hazard_c = hazard_c,
    f = f,
    d = d,
    r = r,
    target_re = target_re,
    stringsAsFactors = FALSE
  )
  n_combi <- nrow(param_grid)

  # --- Pre-allocate result vectors ---
  rho_vec <- rep(NA_real_, n_combi)
  achieved_re_vec <- rep(NA_real_, n_combi)
  converged_vec <- rep(NA, n_combi)
  objective_vec <- rep(NA_real_, n_combi)
  iterations_vec <- rep(NA_integer_, n_combi)
  status_vec <- rep(NA_character_, n_combi)
  error_message_vec <- rep(NA_character_, n_combi)

  # --- Iterate over parameter combinations ---
  for (i in seq_len(n_combi)) {
    current_result <- tryCatch(
      find_rho_for_target_re(
        HR = param_grid$HR[i],
        hazard_c = param_grid$hazard_c[i],
        a_k = a_k,
        omega_k = omega_k,
        f = param_grid$f[i],
        d_t = param_grid$d[i],
        d_c = param_grid$d[i],
        r = param_grid$r[i],
        alpha = alpha,
        power = power,
        time_unit = time_unit,
        metric = metric,
        target_re = param_grid$target_re[i],
        rho_range = rho_range,
        tolerance = tolerance,
        max_iter = max_iter,
        formula = formula
      ),
      error = function(e) {
        list(error = TRUE, message = conditionMessage(e))
      }
    )

    if (!is.null(current_result$error) && current_result$error) {
      status_vec[i] <- "Error"
      error_message_vec[i] <- current_result$message
      converged_vec[i] <- FALSE
    } else {
      rho_vec[i] <- current_result$rho
      achieved_re_vec[i] <- current_result$re
      converged_vec[i] <- current_result$converged
      objective_vec[i] <- current_result$objective
      iterations_vec[i] <- current_result$iterations

      if (isTRUE(current_result$converged) &&
          !is.na(current_result$objective) &&
          current_result$objective < tolerance * 10) {
        status_vec[i] <- "Success"
      } else {
        status_vec[i] <- "Failed"
      }
    }

    if (verbose && n_combi > 50 && i %% 10 == 0) {
      message("Processed ", i, " of ", n_combi, " combinations")
    }
  }

  # --- Assemble results data frame ---
  results_df <- data.frame(
    HR = param_grid$HR,
    hazard_c = param_grid$hazard_c,
    mst_c = log(2) / param_grid$hazard_c,
    f = param_grid$f,
    d = param_grid$d,
    r = param_grid$r,
    target_re = param_grid$target_re,
    metric = metric,
    rho = rho_vec,
    achieved_re = achieved_re_vec,
    converged = converged_vec,
    objective = objective_vec,
    iterations = iterations_vec,
    status = status_vec,
    error_message = error_message_vec,
    stringsAsFactors = FALSE
  )

  # --- Print summary if verbose ---
  if (verbose) {
    n_success <- sum(results_df$status == "Success", na.rm = TRUE)
    n_failed <- sum(results_df$status == "Failed", na.rm = TRUE)
    n_error <- sum(results_df$status == "Error", na.rm = TRUE)
    message("Summary of find_rho_for_target_re_grid:")
    message("  Total scenarios evaluated: ", n_combi)
    message("  Success: ", n_success)
    message("  Failed:  ", n_failed)
    message("  Error:   ", n_error)
    message("  Metric:  ", metric)
    message("  Formula: ", formula)
    if (n_success > 0) {
      successful_rhos <- results_df$rho[results_df$status == "Success"]
      message("  Rho range (successful): [",
              round(min(successful_rhos, na.rm = TRUE), 3), ", ",
              round(max(successful_rhos, na.rm = TRUE), 3), "]")
      message("  Mean rho (successful):  ",
              round(mean(successful_rhos, na.rm = TRUE), 3))
    }
  }

  return(results_df)
}
