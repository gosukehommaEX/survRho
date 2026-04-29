#' Calculate Relative Efficiency of Unequal Randomization in Survival Trials
#'
#' Calculates the relative efficiency (RE) of r:1 randomization compared to
#' 1:1 randomization for a given set of trial design parameters. The function
#' evaluates efficiency across three outcomes: total sample size (N), accrual
#' duration (a_K), and total trial duration (tau). For each r, the accrual
#' rates are multiplied by the corresponding rho value.
#'
#' @param HR A hazard ratio under the alternative hypothesis.
#' @param hazard_c A hazard rate for the control group (per time unit).
#' @param a_k A numeric vector of accrual time interval boundaries (in the
#'   specified time unit).
#' @param omega_k A numeric vector of accrual intensities (patients per time
#'   unit) for each accrual interval.
#' @param f A follow-up time after completion of accrual.
#' @param d_t A dropout rate for the experimental group (annual probability).
#' @param d_c A dropout rate for the control group (annual probability).
#' @param r A numeric vector of allocation ratios to the experimental group
#'   (i.e., T:C = r:1).
#' @param rho A numeric vector of accrual rate multipliers corresponding to
#'   each r value. Must have the same length as r.
#' @param alpha A one-sided level of significance.
#' @param power A target power.
#' @param time_unit A character string specifying the time unit. One of
#'   "months" (default), "weeks", or "years".
#' @param formula A character string specifying the formula used to calculate
#'   the required number of events. One of "Schoenfeld" (default) or "Freedman".
#'
#' @return An object of class \code{"relative_efficiency"} inheriting from
#'   \code{"list"}, with components:
#' \describe{
#'   \item{original_df}{Full output of calc_trial_duration() for the baseline
#'     (r = 1, rho = 1) case and each (r, rho) combination. A rho column is
#'     appended for bookkeeping.}
#'   \item{re_df}{Relative efficiency results with columns r, rho, re_n,
#'     re_accrual_duration, and re_tau, where re_x is defined as the ratio of
#'     the baseline (r = 1) value to the r:1 value for outcome x.}
#' }
#' When printed, a human-readable summary is shown instead of the raw list;
#' see \code{\link{print.relative_efficiency}}. The list components
#' remain accessible in the usual way (e.g., \code{result$re_df}).
#'
#' @details
#' The relative efficiency for outcome x is defined as
#' \code{RE_x = value_{r = 1} / value_{r}}.
#' A value of RE_x greater than 1 indicates that r:1 randomization is more
#' efficient than 1:1 randomization for outcome x; a value less than 1
#' indicates that 1:1 randomization is more efficient. For each r value, the
#' accrual rates omega_k are scaled by the corresponding rho before calling
#' \code{\link{calc_trial_duration}}.
#'
#' @note The vectors r and rho must have the same length. This function
#'   requires \code{\link{calc_trial_duration}} and the \code{dplyr} package.
#'
#' @examples
#' \dontrun{
#' result <- calc_relative_efficiency(
#'   HR = 0.7,
#'   hazard_c = log(2) / 12,
#'   a_k = c(0, 3, 6),
#'   omega_k = c(15, 20, 25),
#'   f = 12,
#'   d_t = 0.3,
#'   d_c = 0.3,
#'   r = c(2, 3, 4),
#'   rho = c(1.2, 1.3, 1.4),
#'   alpha = 0.025,
#'   power = 0.9,
#'   time_unit = "months",
#'   formula = "Schoenfeld"
#' )
#' print(result)             # Human-readable summary
#' result$re_df              # Access the RE data frame
#' result$original_df        # Access the full trial-level output
#' }
#'
#' @seealso \code{\link{calc_trial_duration}} for the underlying trial design
#'   calculation; \code{\link{print.relative_efficiency}} for the print method.
#'
#' @importFrom dplyr bind_rows
#'
#' @export
calc_relative_efficiency <- function(HR, hazard_c, a_k, omega_k, f, d_t, d_c,
                                     r, rho, alpha, power,
                                     time_unit = "months",
                                     formula = "Schoenfeld") {
  
  # --- Validate input lengths ---
  if (length(r) != length(rho)) {
    stop("Vectors r and rho must have the same length")
  }
  
  # --- Baseline case (r = 1, rho = 1) ---
  baseline_result <- calc_trial_duration(
    HR = HR,
    hazard_c = hazard_c,
    a_k = a_k,
    omega_k = omega_k,
    f = f,
    d_t = d_t,
    d_c = d_c,
    r = 1,
    alpha = alpha,
    power = power,
    time_unit = time_unit,
    formula = formula
  )
  baseline_result$rho <- 1
  
  baseline_N <- baseline_result$N
  baseline_a_K <- baseline_result$a_K
  baseline_tau <- baseline_result$tau
  
  # --- Loop over (r, rho) combinations ---
  results_list <- vector("list", length(r) + 1)
  results_list[[1]] <- baseline_result
  
  re_n <- numeric(length(r))
  re_accrual_duration <- numeric(length(r))
  re_tau <- numeric(length(r))
  
  for (i in seq_along(r)) {
    current_result <- calc_trial_duration(
      HR = HR,
      hazard_c = hazard_c,
      a_k = a_k,
      omega_k = omega_k * rho[i],
      f = f,
      d_t = d_t,
      d_c = d_c,
      r = r[i],
      alpha = alpha,
      power = power,
      time_unit = time_unit,
      formula = formula
    )
    current_result$rho <- rho[i]
    results_list[[i + 1]] <- current_result
    
    re_n[i] <- baseline_N / current_result$N
    re_accrual_duration[i] <- baseline_a_K / current_result$a_K
    re_tau[i] <- baseline_tau / current_result$tau
  }
  
  # --- Combine results ---
  # bind_rows handles differing column sets (e.g., number of omega_k columns)
  # by filling missing columns with NA.
  original_df <- dplyr::bind_rows(results_list)
  
  re_df <- data.frame(
    r = r,
    rho = rho,
    re_n = re_n,
    re_accrual_duration = re_accrual_duration,
    re_tau = re_tau,
    stringsAsFactors = FALSE
  )
  
  out <- list(original_df = original_df, re_df = re_df)
  class(out) <- c("relative_efficiency", "list")
  return(out)
}