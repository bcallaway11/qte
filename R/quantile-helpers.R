# =============================================================================
# Title: Weighted quantile helper
# Description: Shared weighted type-1 quantile helper used by ddid_gt() and
#   qdid_gt(). Not exported.
# Author: Brant Callaway
# Last update: 2026-07-21
# Date created: 2026-07-21
# =============================================================================

#' @title Weighted type-1 quantile
#'
#' @description Smallest \code{y} with weighted CDF \code{>= p}. Uses a small
#'   numerical tolerance on the comparison so that a \code{p} landing exactly
#'   on a floating-point rounding boundary (common with uniform weights, e.g.
#'   \code{1/n}) picks the intended order statistic instead of flipping based
#'   on last-bit differences in summation across platforms.
#'
#' @param y_vec numeric vector of values.
#' @param w_vec numeric vector of weights, same length as \code{y_vec}.
#' @param probs numeric vector of probabilities in \code{[0, 1]}.
#'
#' @return numeric vector, same length as \code{probs}.
#'
#' @keywords internal
wquant <- function(y_vec, w_vec, probs) {
  ord <- order(y_vec)
  y_s <- y_vec[ord]
  cdf <- cumsum(w_vec[ord])
  tol <- sqrt(.Machine$double.eps)
  y_s[vapply(probs, function(p) {
    idx <- which(cdf >= p - tol)
    if (length(idx) == 0L) length(y_s) else idx[1L]
  }, integer(1))]
}
