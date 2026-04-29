#' Find Accrual Rate Multiplier for a Target Relative Efficiency
#'
#' Finds the value of rho (accrual rate multiplier) at which the relative
#' efficiency (RE) of r:1 randomization compared to 1:1 randomization equals
#' a specified target value for a chosen outcome. The root is found with
#' \code{\link[stats]{uniroot}} (Brent's method).
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
#' @param r A scalar allocation ratio to the experimental group (T:C = r:1).
#' @param alpha A one-sided level of significance.
#' @param power A target power.
#' @param time_unit A character string specifying the time unit. One of
#'   "months" (default), "weeks", or "years".
#' @param metric A character string specifying the outcome for which the
#'   target RE is sought. One of "n" (total sample size),
#'   "accrual_duration" (accrual duration a_K), or "tau" (total trial
#'   duration). Default is "n".
#' @param target_re A target relative efficiency value. Default is 1.
#' @param rho_range A numeric vector of length 2 giving the search range for
#'   rho. Default is c(0.1, 5).
#' @param tolerance A numeric value specifying the tolerance for the root
#'   finder. Default is 1e-6.
#' @param max_iter A maximum number of iterations for the root finder.
#'   Default is 1000.
#' @param formula A character string specifying the formula used to calculate
#'   the required number of events. One of "Schoenfeld" (default) or "Freedman".
#'
#' @return An object of class \code{"rho_search"} inheriting from
#'   \code{"list"}, with components:
#' \describe{
#'   \item{rho}{The rho value at which RE equals target_re.}
#'   \item{re}{The relative efficiency achieved at the returned rho.}
#'   \item{target_re}{The target relative efficiency.}
#'   \item{metric}{The outcome used for optimization.}
#'   \item{r}{The allocation ratio used.}
#'   \item{converged}{Logical indicating whether the root finder converged.}
#'   \item{objective}{The absolute difference |re - target_re| at the returned rho.}
#'   \item{iterations}{The number of iterations used by the root finder.}
#' }
#' When printed, a human-readable summary is shown; see
#' \code{\link{print.rho_search}}. The list components remain accessible in
#' the usual way (e.g., \code{result$rho}).
#'
#' @details
#' The relative efficiency is defined as \code{RE = value_{r = 1} / value_{r}} for
#' the chosen outcome. The function solves RE(rho) - target_re = 0 using
#' \code{\link[stats]{uniroot}}. Before the root search, the function checks
#' whether target_re lies within the interval [RE(rho_range[1]), RE(rho_range[2])]
#' and issues an informative error message if not.
#'
#' @note This function requires \code{\link{calc_relative_efficiency}}.
#'
#' @examples
#' \dontrun{
#' result <- find_rho_for_target_re(
#'   HR = 0.7,
#'   hazard_c = log(2) / 12,
#'   a_k = c(0, 3, 6),
#'   omega_k = c(15, 20, 25),
#'   f = 12,
#'   d_t = 0.3,
#'   d_c = 0.3,
#'   r = 2,
#'   alpha = 0.025,
#'   power = 0.9,
#'   metric = "tau",
#'   target_re = 1,
#'   formula = "Schoenfeld"
#' )
#' print(result)
#' result$rho
#' }
#'
#' @seealso \code{\link{calc_relative_efficiency}} for the underlying
#'   relative efficiency calculation; \code{\link{print.rho_search}} for the
#'   print method.
#'
#' @importFrom stats uniroot
#' @export
find_rho_for_target_re <- function(HR, hazard_c, a_k, omega_k, f, d_t, d_c,
                                   r, alpha, power, time_unit = "months",
                                   metric = "n", target_re = 1,
                                   rho_range = c(0.1, 5),
                                   tolerance = 1e-6, max_iter = 1000,
                                   formula = "Schoenfeld") {
  
  # --- Validate inputs ---
  valid_metrics <- c("n", "accrual_duration", "tau")
  if (!metric %in% valid_metrics) {
    stop("metric must be one of: ", paste(valid_metrics, collapse = ", "))
  }
  
  if (!is.numeric(target_re) || length(target_re) != 1 || target_re <= 0) {
    stop("target_re must be a single positive numeric value")
  }
  
  if (length(rho_range) != 2 || rho_range[1] >= rho_range[2]) {
    stop("rho_range must be a numeric vector of length 2 with rho_range[1] < rho_range[2]")
  }
  if (any(rho_range <= 0)) stop("rho_range values must be positive")
  
  if (tolerance <= 0) stop("tolerance must be positive")
  if (max_iter <= 0) stop("max_iter must be a positive integer")
  
  # --- Map metric name to re_df column name ---
  re_column <- switch(metric,
                      "n" = "re_n",
                      "accrual_duration" = "re_accrual_duration",
                      "tau" = "re_tau")
  
  # --- Helper: calculate RE at a given rho ---
  calc_re_at_rho <- function(rho_val) {
    if (!is.finite(rho_val) || rho_val <= 0) return(NA_real_)
    result <- calc_relative_efficiency(
      HR = HR,
      hazard_c = hazard_c,
      a_k = a_k,
      omega_k = omega_k,
      f = f,
      d_t = d_t,
      d_c = d_c,
      r = r,
      rho = rho_val,
      alpha = alpha,
      power = power,
      time_unit = time_unit,
      formula = formula
    )
    result$re_df[[re_column]]
  }
  
  # --- Feasibility check at endpoints ---
  re_lower <- calc_re_at_rho(rho_range[1])
  re_upper <- calc_re_at_rho(rho_range[2])
  
  if (is.na(re_lower) || is.na(re_upper)) {
    stop("Cannot evaluate relative efficiency at the given rho_range. ",
         "Please adjust rho_range.")
  }
  
  if (target_re < min(re_lower, re_upper) || target_re > max(re_lower, re_upper)) {
    stop(sprintf(paste("Target RE (%.3f) is not achievable within rho_range [%.3f, %.3f].",
                       "Achievable range: [%.3f, %.3f].",
                       "Consider: (1) adjusting target_re,",
                       "(2) expanding rho_range,",
                       "(3) changing other trial parameters."),
                 target_re, rho_range[1], rho_range[2],
                 min(re_lower, re_upper), max(re_lower, re_upper)))
  }
  
  # --- Objective function: RE(rho) - target_re ---
  objective_fun <- function(rho_val) calc_re_at_rho(rho_val) - target_re
  
  # --- Solve for rho using Brent's method ---
  root_result <- tryCatch(
    uniroot(objective_fun,
            lower = rho_range[1],
            upper = rho_range[2],
            tol = tolerance,
            maxiter = max_iter),
    error = function(e) NULL
  )
  
  if (is.null(root_result)) {
    out <- list(
      rho = NA_real_,
      re = NA_real_,
      target_re = target_re,
      metric = metric,
      r = r,
      converged = FALSE,
      objective = NA_real_,
      iterations = NA_integer_
    )
    class(out) <- c("rho_search", "list")
    return(out)
  }
  
  optimal_rho <- root_result$root
  final_re <- calc_re_at_rho(optimal_rho)
  final_objective <- abs(final_re - target_re)
  # uniroot considers convergence successful if it returns without error
  # and the estimated precision is within tolerance
  converged <- !is.null(root_result$estim.prec) &&
    root_result$estim.prec <= tolerance * 10
  
  out <- list(
    rho = optimal_rho,
    re = final_re,
    target_re = target_re,
    metric = metric,
    r = r,
    converged = converged,
    objective = final_objective,
    iterations = root_result$iter
  )
  class(out) <- c("rho_search", "list")
  return(out)
}