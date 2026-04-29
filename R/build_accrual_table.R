#' Build a Compact Accrual Table from a trial_duration Object
#'
#' Internal helper that extracts the accrual intervals (omega_k and interval
#' boundaries) from a \code{\link{calc_trial_duration}} result and returns a
#' compact data frame with one row per interval.
#'
#' @param x An object of class \code{"trial_duration"} (or a data frame with
#'   the same structure).
#' @param omega_cols A character vector of column names for omega_k values
#'   (matching \code{^omega_[0-9]+$}).
#' @param a_cols A character vector of column names for a_k boundaries
#'   (matching \code{^a_[0-9]+$}).
#' @param p_cols A character vector of column names for p_k proportions
#'   (matching \code{^p_[0-9]+$}); may be empty if p_k was not specified.
#' @param digits Integer. Number of digits to round to in the returned table.
#'
#' @return A data frame with columns \code{interval}, \code{start},
#'   \code{end}, \code{omega_k}, and optionally \code{p_k}. Returns
#'   \code{NULL} if the lengths of omega_cols and a_cols are inconsistent.
#'
#' @details
#' In calc_trial_duration()'s output, there are K values of omega_k (one per
#' interval) and typically K+1 values of a_k (the interval boundaries
#' a_1, ..., a_K, a_{K+1} = a_K[final]). The k-th interval runs from
#' a_k[k] to a_k[k+1]. When length(a_vec) == K (a fallback case), the first
#' start is taken to be 0.
#'
#' @keywords internal
#' @noRd
build_accrual_table <- function(x, omega_cols, a_cols, p_cols, digits) {
  omega_vec <- as.numeric(x[1, omega_cols])
  a_vec <- as.numeric(x[1, a_cols])
  K <- length(omega_vec)

  # Determine start / end for each of the K intervals
  if (length(a_vec) == K + 1) {
    start_vec <- a_vec[seq_len(K)]
    end_vec <- a_vec[seq_len(K) + 1]
  } else if (length(a_vec) == K) {
    start_vec <- c(0, a_vec[-K])
    end_vec <- a_vec
  } else {
    return(NULL)
  }

  tbl <- data.frame(
    interval = seq_len(K),
    start = round(start_vec, digits),
    end = round(end_vec, digits),
    omega_k = round(omega_vec, digits),
    stringsAsFactors = FALSE
  )

  if (length(p_cols) == K) {
    p_vec <- as.numeric(x[1, p_cols])
    tbl$p_k <- round(p_vec, digits)
  }

  tbl
}
