#' Analytical Bounds on the Required Accrual Rate Multiplier
#'
#' Computes theoretical bounds on the accrual rate multiplier rho that is
#' required to maintain equal total trial duration (i.e., RE = 1) when
#' switching from 1:1 to r:1 randomization in a two-group survival trial.
#' Under the exponential survival / exponential dropout framework of
#' Lachin and Foulkes (1986) and the event-count formula of Schoenfeld
#' (1981), rho satisfies
#'
#' \deqn{\rho = \frac{(1 + r)^3}{8 r} \cdot \frac{1 + \theta(\tau)}{1 + r \theta(\tau)}}
#'
#' where \eqn{\theta(\tau) = \phi_T(\tau) / \phi_C(\tau)} is the ratio of
#' event probabilities in the two groups. This function returns three
#' analytical bounds on rho that do not require numerical root-finding:
#' (i) a tight lower bound and an upper bound that depend on
#' (\eqn{HR}, \eqn{\lambda_C}, \eqn{\gamma}, \eqn{r}); and (ii) a loose
#' lower bound \eqn{(r - 1)^2 / (4 r) + 1} that depends on \eqn{r} only.
#'
#' @param HR A hazard ratio under the alternative hypothesis. Must satisfy
#'   \code{HR < 1} (experimental better than control).
#' @param hazard_c A hazard rate for the control group (per time unit).
#' @param d A dropout rate common to both groups (annual probability). The
#'   underlying theory assumes a single dropout hazard gamma shared across
#'   groups, so only one value is accepted here (unlike
#'   \code{\link{calc_trial_duration}} which allows group-specific dropout).
#' @param r A numeric vector of allocation ratios to the experimental group
#'   (i.e., T:C = r:1). All values must satisfy \code{r > 1}.
#' @param time_unit A character string specifying the time unit. One of
#'   "months" (default), "weeks", or "years".
#'
#' @return An object of class \code{"rho_bounds"} inheriting from
#'   \code{"list"}, with components:
#' \describe{
#'   \item{bounds_df}{A data frame with columns r, rho_lower_loose,
#'     rho_lower_tight, rho_upper, and theta_U, giving the three bounds
#'     and the corresponding upper bound on theta for each r.}
#'   \item{HR}{The input hazard ratio (echoed).}
#'   \item{hazard_c}{The input control-group hazard rate (echoed).}
#'   \item{d}{The input annual dropout probability (echoed).}
#'   \item{gamma}{The dropout hazard per time unit derived from d.}
#'   \item{time_unit}{The input time unit (echoed).}
#' }
#' When printed, a human-readable summary is shown instead of the raw list;
#' see \code{\link{print.rho_bounds}}. The list components remain accessible
#' in the usual way (e.g., \code{result$bounds_df}).
#'
#' @details
#' The three bounds returned are:
#' \itemize{
#'   \item \strong{Loose lower bound} (depends on r only):
#'     \eqn{\rho \ge (r - 1)^2 / (4 r) + 1}. This universal bound holds for
#'     any HR, any control hazard, any dropout rate, any accrual duration,
#'     and any follow-up duration. It equals the Schoenfeld event-count
#'     ratio \eqn{E_r / E_1}.
#'   \item \strong{Tight lower bound}:
#'     \eqn{\rho \ge \frac{(1 + r)^3}{8 r} \cdot \frac{1 + \theta_U}{1 + r \theta_U}},
#'     where
#'     \eqn{\theta_U = HR (\lambda_C + \gamma) / (HR \lambda_C + \gamma)}.
#'   \item \strong{Upper bound}:
#'     \eqn{\rho \le \frac{(1 + r)^3}{8 r} \cdot \frac{1 + HR}{1 + r HR}}.
#' }
#' The bounds follow from monotonicity of the map
#' \eqn{\theta \mapsto \rho(\theta)} (decreasing for r > 1) combined with
#' the inequality \eqn{HR \le \theta(\tau) \le \theta_U}. The loose lower
#' bound follows from the additional observation that
#' \eqn{(1 + \theta_U) / (1 + r \theta_U) \ge 2 / (1 + r)} whenever
#' \eqn{\theta_U \le 1}.
#'
#' @note The annual dropout probability \code{d} is converted to a hazard
#'   rate \code{gamma} per time unit using \code{gamma = -log(1 - d) / C},
#'   where \code{C} equals 12 for months, 365.25 / 7 for weeks, and 1 for
#'   years. This matches the convention used by
#'   \code{\link{calc_trial_duration}}.
#'
#' @examples
#' \dontrun{
#' # Phase 3 oncology setting (HR = 0.8, control median survival = 12 months,
#' # 30% annual dropout)
#' result <- calc_rho_bounds(
#'   HR = 0.8,
#'   hazard_c = log(2) / 12,
#'   d = 0.3,
#'   r = c(2, 3, 4),
#'   time_unit = "months"
#' )
#' print(result)              # Human-readable summary
#' result$bounds_df           # Access the bounds data frame
#' }
#'
#' @seealso \code{\link{calc_relative_efficiency}} for the numerical
#'   calculation of rho via root-finding; \code{\link{print.rho_bounds}}
#'   for the print method.
#'
#' @export
calc_rho_bounds <- function(HR, hazard_c, d, r, time_unit = "months") {

  # --- Validate inputs ---
  valid_units <- c("months", "weeks", "years")
  if (!time_unit %in% valid_units) {
    stop("time_unit must be one of: ", paste(valid_units, collapse = ", "))
  }
  if (length(HR) != 1 || !is.numeric(HR)) {
    stop("HR must be a single numeric value")
  }
  if (HR <= 0 || HR >= 1) {
    stop("HR must be in (0, 1); theoretical bounds assume HR < 1")
  }
  if (length(hazard_c) != 1 || !is.numeric(hazard_c) || hazard_c <= 0) {
    stop("hazard_c must be a single positive numeric value")
  }
  if (length(d) != 1 || !is.numeric(d) || d < 0 || d >= 1) {
    stop("d must be a single numeric value in [0, 1)")
  }
  if (!is.numeric(r) || any(r <= 1)) {
    stop("All values in r must be numeric and greater than 1")
  }

  # --- Convert annual dropout probability to hazard per time unit ---
  # Matches the convention of calc_trial_duration()
  time_conversion <- switch(time_unit,
                            "months" = 12,
                            "weeks"  = 365.25 / 7,
                            "years"  = 1)
  gamma <- -log(1 - d) / time_conversion

  # --- Derived quantities that do not depend on r ---
  # theta_U = upper bound on theta(tau)
  theta_U <- HR * (hazard_c + gamma) / (HR * hazard_c + gamma)

  # --- Compute bounds for each r ---
  # rho_lower_loose : (r - 1)^2 / (4 r) + 1        [depends on r only]
  # rho_lower_tight : (1 + r)^3 / (8 r) * (1 + theta_U) / (1 + r * theta_U)
  # rho_upper       : (1 + r)^3 / (8 r) * (1 + HR) / (1 + r * HR)
  coef <- (1 + r) ^ 3 / (8 * r)
  rho_lower_loose <- (r - 1) ^ 2 / (4 * r) + 1
  rho_lower_tight <- coef * (1 + theta_U) / (1 + r * theta_U)
  rho_upper       <- coef * (1 + HR)      / (1 + r * HR)

  bounds_df <- data.frame(
    r = r,
    rho_lower_loose = rho_lower_loose,
    rho_lower_tight = rho_lower_tight,
    rho_upper = rho_upper,
    theta_U = rep(theta_U, length(r)),
    stringsAsFactors = FALSE
  )

  # --- Assemble return value ---
  out <- list(
    bounds_df = bounds_df,
    HR = HR,
    hazard_c = hazard_c,
    d = d,
    gamma = gamma,
    time_unit = time_unit
  )
  class(out) <- c("rho_bounds", "list")
  return(out)
}
