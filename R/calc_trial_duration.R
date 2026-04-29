#' Calculate Accrual Duration and Total Trial Duration for Survival Trials
#'
#' Calculates the required accrual duration and total trial duration for a
#' survival trial with a time-to-event endpoint, comparing an experimental
#' group (T) against a control group (C). The function determines the accrual
#' duration (a_K) that yields the target number of events and then computes
#' the total trial duration (tau = a_K + f). Required sample size and number
#' of events are obtained from the Lachin and Foulkes (1986) formula, assuming
#' piecewise uniform accrual and exponential survival distributions.
#'
#' @param HR A hazard ratio under the alternative hypothesis.
#' @param hazard_c A hazard rate for the control group (per time unit).
#' @param a_k A numeric vector of accrual time interval boundaries (in the
#'   specified time unit). The first element is typically 0.
#' @param omega_k A numeric vector of accrual intensities (patients per time
#'   unit) for each accrual interval. If NULL, p_k must be specified.
#' @param p_k A numeric vector of proportions of total enrollment for each
#'   accrual interval (must sum to 1). If specified, omega_k must be NULL.
#' @param f A follow-up time after completion of accrual (in the specified
#'   time unit).
#' @param d_t A dropout rate for the experimental group (annual probability).
#' @param d_c A dropout rate for the control group (annual probability).
#' @param r An allocation ratio to the experimental group (i.e., T:C = r:1).
#' @param alpha A one-sided level of significance.
#' @param power A target power.
#' @param time_unit A character string specifying the time unit. One of
#'   "months" (default), "weeks", or "years".
#' @param formula A character string specifying the formula used to calculate
#'   the required number of events. One of "Schoenfeld" (default) or "Freedman".
#'
#' @return An object of class \code{"trial_duration"} inheriting from
#'   \code{"data.frame"}, with one row and the following columns:
#' \describe{
#'   \item{HR}{Hazard ratio}
#'   \item{mst_t, mst_c}{Median survival times for experimental and control groups}
#'   \item{hazard_t, hazard_c}{Hazard rates for experimental and control groups}
#'   \item{d_t, d_c}{Dropout rates for experimental and control groups (annual probability)}
#'   \item{f}{Follow-up time after completion of accrual}
#'   \item{r}{Allocation ratio (experimental:control)}
#'   \item{alpha}{One-sided significance level}
#'   \item{power}{Target power}
#'   \item{E}{Required number of events}
#'   \item{n_t, n_c}{Required sample sizes for experimental and control groups}
#'   \item{N}{Total required sample size}
#'   \item{a_K}{Total accrual duration (primary output)}
#'   \item{tau}{Total trial duration, tau = a_K + f (primary output)}
#'   \item{time_unit}{Time unit used for the calculations}
#'   \item{formula}{Formula used for event calculation}
#'   \item{omega_1, omega_2, ...}{Accrual intensities for each accrual interval}
#'   \item{a_1, a_2, ...}{Accrual time boundaries including the final a_K}
#'   \item{p_1, p_2, ...}{Proportions for each accrual interval (only if p_k was specified)}
#' }
#' The object can be used as a regular data frame (e.g., \code{result$N}).
#' When printed, a human-readable summary is shown instead of the raw data
#' frame; see \code{\link{print.trial_duration}}.
#'
#' @details
#' The function can be used in two modes depending on whether omega_k or p_k
#' is specified.
#'
#' Mode 1: Specify omega_k directly. The function uses the given accrual
#' intensities and solves for the accrual duration a_K that achieves the
#' target number of events.
#'
#' Mode 2: Specify p_k (proportions). The function first uses provisional
#' uniform accrual rates to find a_K, then back-calculates omega_k so that
#' the specified proportions are exactly achieved, and finally recomputes the
#' sample sizes.
#'
#' The accrual duration a_K is obtained as the root of the equation
#' \code{sum_{k} omega_k * (a_k - a_{k-1}) * (r * Phi_T(a_K) + Phi_C(a_K)) / (1 + r) = E},
#' where Phi_j(a_K) is the per-interval event occurrence probability derived
#' from Lachin and Foulkes (1986). The root is obtained with
#' \code{\link[stats]{uniroot}} (Brent's method).
#'
#' This function assumes one-sided tests and exponential survival
#' distributions. Dropout rates are specified as annual probabilities and are
#' internally converted to the specified time unit (months: 12, weeks: 365.25/7,
#' years: 1).
#'
#' @note
#' Either omega_k or p_k must be specified, but not both.
#'
#' @examples
#' \dontrun{
#' # Mode 1: Specify omega_k directly
#' result1 <- calc_trial_duration(
#'   HR = 0.7,
#'   hazard_c = log(2) / 12,
#'   a_k = c(0, 3, 6),
#'   omega_k = c(15, 20, 25),
#'   f = 12,
#'   d_t = 0.3,
#'   d_c = 0.3,
#'   r = 1,
#'   alpha = 0.025,
#'   power = 0.9,
#'   time_unit = "months",
#'   formula = "Schoenfeld"
#' )
#' print(result1)      # Human-readable summary
#' result1$N           # Access as data frame
#'
#' # Mode 2: Specify proportions p_k
#' result2 <- calc_trial_duration(
#'   HR = 0.7,
#'   hazard_c = log(2) / 12,
#'   a_k = c(0, 3, 6),
#'   omega_k = NULL,
#'   p_k = c(0.2, 0.3, 0.5),
#'   f = 12,
#'   d_t = 0.3,
#'   d_c = 0.3,
#'   r = 1,
#'   alpha = 0.025,
#'   power = 0.9,
#'   time_unit = "months",
#'   formula = "Schoenfeld"
#' )
#' }
#'
#' @references
#' Lachin, J. M., and Foulkes, M. A. (1986). Evaluation of sample size and
#' power for analyses of survival with allowance for nonuniform patient entry,
#' losses to follow-up, noncompliance, and stratification. Biometrics, 42(3),
#' 507-519.
#'
#' Freedman, L. S. (1982). Tables of the number of patients required in
#' clinical trials using the logrank test. Statistics in Medicine, 1(2),
#' 121-129.
#'
#' Brent, R. P. (1973). Algorithms for Minimization without Derivatives.
#' Prentice-Hall.
#'
#' @seealso \code{\link{print.trial_duration}} for the print method.
#'
#' @importFrom stats qnorm setNames uniroot
#' @export
calc_trial_duration <- function(HR, hazard_c, a_k, omega_k = NULL, p_k = NULL,
                                f, d_t, d_c, r, alpha, power,
                                time_unit = "months", formula = "Schoenfeld") {
  
  # --- Validate scalar inputs ---
  valid_units <- c("months", "weeks", "years")
  if (!time_unit %in% valid_units) {
    stop("time_unit must be one of: ", paste(valid_units, collapse = ", "))
  }
  
  valid_formulas <- c("Schoenfeld", "Freedman")
  if (!formula %in% valid_formulas) {
    stop("formula must be one of: ", paste(valid_formulas, collapse = ", "))
  }
  
  if (HR <= 0) stop("HR must be positive")
  if (hazard_c <= 0) stop("hazard_c must be positive")
  if (f < 0) stop("f must be non-negative")
  if (d_t < 0 || d_t >= 1) stop("d_t must be in [0, 1)")
  if (d_c < 0 || d_c >= 1) stop("d_c must be in [0, 1)")
  if (r <= 0) stop("r must be positive")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be in (0, 1)")
  if (power <= 0 || power >= 1) stop("power must be in (0, 1)")
  
  # --- Validate a_k ---
  if (any(diff(a_k) <= 0)) stop("a_k must be strictly increasing")
  
  # --- Validate omega_k / p_k exclusivity and set mode ---
  if (is.null(omega_k) && is.null(p_k)) {
    stop("Either omega_k or p_k must be specified")
  }
  if (!is.null(omega_k) && !is.null(p_k)) {
    stop("Only one of omega_k or p_k can be specified, not both")
  }
  
  original_p_k <- NULL
  if (!is.null(omega_k)) {
    if (length(a_k) != length(omega_k)) {
      stop("Length of a_k and omega_k must be equal")
    }
    if (any(omega_k <= 0)) stop("All values in omega_k must be positive")
    use_proportions <- FALSE
  } else {
    if (length(a_k) != length(p_k)) {
      stop("Length of a_k and p_k must be equal")
    }
    if (any(p_k <= 0)) stop("All values in p_k must be positive")
    if (abs(sum(p_k) - 1) > 1e-10) stop("Values in p_k must sum to 1")
    
    original_p_k <- p_k
    # Provisional uniform rates; scale roughly with the expected sample size
    initial_rate <- if (HR <= 0.4) 5
    else if (HR <= 0.6) 15
    else if (HR <= 0.8) 50
    else 100
    omega_k <- rep(initial_rate, length(a_k))
    use_proportions <- TRUE
  }
  
  # --- Derived quantities ---
  time_conversion <- switch(time_unit,
                            "months" = 12,
                            "weeks"  = 365.25 / 7,
                            "years"  = 1)
  
  # Dropout hazard rate per time unit (j = T, C)
  d_j <- c(d_t, d_c)
  gamma_j <- -log(1 - d_j) / time_conversion
  
  # Event hazard rate per time unit (j = T, C)
  hazard_t <- hazard_c * HR
  lambda_j <- c(hazard_t, hazard_c)
  
  # Median survival time
  mst_j <- log(2) / lambda_j
  
  # Required number of events
  z_sum <- qnorm(power) + qnorm(1 - alpha)
  E <- if (formula == "Schoenfeld") {
    (z_sum / log(HR))^2 * ((1 + r)^2 / r)
  } else {
    (z_sum * (1 + r * HR) / (1 - HR))^2 * (1 / r)
  }
  
  # --- Event occurrence probability function ---
  # Returns a list with elements "T" and "C"; each is a vector of length K
  # giving Phi_j on each accrual interval (k = 1, ..., K).
  phi_fun <- function(a_K_val) {
    a_k_ext <- c(a_k, a_K_val)
    lapply(setNames(seq_len(2), c("T", "C")), function(j) {
      lam <- lambda_j[j] + gamma_j[j]
      sapply(seq_along(a_k_ext)[-1], function(k) {
        xi_k <- exp(-lam * (a_K_val + f - a_k_ext[k])) -
          exp(-lam * (a_K_val + f - a_k_ext[k - 1]))
        (lambda_j[j] / lam) *
          (1 - xi_k / (lam * (a_k_ext[k] - a_k_ext[k - 1])))
      })
    })
  }
  
  # Objective function: expected total events minus target E.
  # The root a_K satisfies event_fun(a_K) = 0.
  event_fun <- function(a_K_val) {
    phi_vals <- phi_fun(a_K_val)
    sum(omega_k * diff(c(a_k, a_K_val)) *
          (r * phi_vals[["T"]] + phi_vals[["C"]]) / (1 + r)) - E
  }
  
  # --- Feasibility check ---
  a_K_lower <- max(a_k) + 1e-4
  a_K_upper <- 1e+4
  
  val_upper <- event_fun(a_K_upper)
  if (val_upper < 0) {
    stop(sprintf(paste("Cannot achieve the required number of events (E = %.1f)",
                       "even with very long accrual duration.",
                       "Maximum achievable events is approximately %.1f.",
                       "Consider: (1) reducing power or increasing alpha,",
                       "(2) increasing accrual rates (omega_k),",
                       "(3) extending follow-up period (f),",
                       "(4) adjusting hazard ratio (HR)."),
                 E, val_upper + E))
  }
  val_lower <- event_fun(a_K_lower)
  if (val_lower > 0) {
    stop("The target events are already achieved at the minimum accrual duration. ",
         "Check input parameters.")
  }
  
  # --- Solve for a_K using Brent's method ---
  root_result <- uniroot(event_fun,
                         lower = a_K_lower,
                         upper = a_K_upper,
                         tol = .Machine$double.eps^0.25,
                         maxiter = 1000)
  a_K <- root_result$root
  
  if (a_K <= max(a_k)) {
    stop(sprintf("Calculated a_K (%.3f) is not larger than max(a_k) (%.3f). ",
                 a_K, max(a_k)),
         "Check input parameters.")
  }
  
  # --- In proportions mode, back-calculate omega_k from p_k and a_K ---
  if (use_proportions) {
    current_N <- sum(omega_k * diff(c(a_k, a_K)))
    K <- length(a_k)
    new_omega_k <- numeric(K)
    for (k in seq_len(K - 1)) {
      n_k <- current_N * original_p_k[k]
      duration_k <- a_k[k + 1] - a_k[k]
      new_omega_k[k] <- n_k / duration_k
    }
    n_last <- current_N * original_p_k[K]
    duration_last <- a_K - a_k[K]
    new_omega_k[K] <- n_last / duration_last
    omega_k <- new_omega_k
    p_k <- original_p_k
  }
  
  # --- Final sample sizes ---
  n_j <- c(r, 1) / (1 + r) * sum(omega_k * diff(c(a_k, a_K)))
  N <- sum(n_j)
  
  if (any(n_j < 0) || N < 0) {
    stop("Negative sample sizes detected. Check input parameters.")
  }
  
  # --- Assemble output ---
  output_data <- data.frame(
    HR = HR,
    mst_t = mst_j[1],
    mst_c = mst_j[2],
    hazard_t = lambda_j[1],
    hazard_c = lambda_j[2],
    d_t = d_j[1],
    d_c = d_j[2],
    f = f,
    r = r,
    alpha = alpha,
    power = power,
    E = E,
    n_t = n_j[1],
    n_c = n_j[2],
    N = N,
    a_K = a_K,
    tau = a_K + f,
    time_unit = time_unit,
    formula = formula,
    stringsAsFactors = FALSE
  )
  
  # Add omega_k columns
  for (k in seq_along(omega_k)) {
    output_data[[paste0("omega_", k)]] <- omega_k[k]
  }
  
  # Add a_k columns (include final a_K boundary)
  final_a_k <- c(a_k[a_k < a_K], a_K)
  for (k in seq_along(final_a_k)) {
    output_data[[paste0("a_", k)]] <- final_a_k[k]
  }
  
  # Add p_k columns if proportions mode was used
  if (use_proportions) {
    for (k in seq_along(original_p_k)) {
      output_data[[paste0("p_", k)]] <- original_p_k[k]
    }
  }
  
  # Attach S3 class; the object remains usable as a data frame
  class(output_data) <- c("trial_duration", class(output_data))
  return(output_data)
}