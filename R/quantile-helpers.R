# =============================================================================
# Title: Quantile helpers
# Description: Shared, unexported quantile helpers: a weighted type-1
#   quantile (used by ddid_gt() and qdid_gt()), and conditional CDF/quantile
#   evaluation from quantile regressions with rearrangement (used by the
#   covariate-adjusted CiC, DDID and panel QTT estimators).
# Author: Brant Callaway
# Last update: 2026-11-14
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

# --- Conditional quantiles from quantile regression ---------------------------

#' @title Rearranged conditional quantiles from a quantile regression
#'
#' @description Predicted conditional quantiles at each row of
#'   \code{newdata}, with each row sorted. Sorting is the rearrangement of
#'   Chernozhukov, Fernandez-Val and Galichon (2010): it makes the conditional
#'   quantile function monotone when the quantile regressions cross, and
#'   leaves it unchanged when they don't.
#'
#'   Used instead of \code{predict(fit, type = "Fhat"/"Qhat", stepfun = TRUE)},
#'   which fails when \code{newdata} has a single row and returns a
#'   non-monotone \code{Fhat} when quantiles cross (issue #9).
#'
#' @param fit an \code{rqs} object (\code{rq()} with a vector of \code{tau}).
#' @param newdata data frame of covariate values; if missing, the data used
#'   to fit the model.
#'
#' @return n x length(fit$tau) matrix; row i holds the sorted predicted
#'   quantiles for row i of \code{newdata}.
#'
#' @references Chernozhukov, V., I. Fernandez-Val and A. Galichon (2010).
#'   ``Quantile and Probability Curves Without Crossing.'' Econometrica
#'   78(3), 1093-1125.
#'
#' @keywords internal
qr_sorted_quantiles <- function(fit, newdata) {
  q <- if (missing(newdata)) predict(fit) else predict(fit, newdata = newdata)
  q <- matrix(q, ncol = length(fit$tau))
  # apply() returns one column per row of q; transpose back to n x K
  matrix(t(apply(q, 1, sort)), ncol = ncol(q))
}

#' @title Conditional CDF and quantile function from sorted QR predictions
#'
#' @description Step-function evaluation matching quantreg's
#'   \code{predict(..., stepfun = TRUE)} conventions:
#'   \code{qr_cdf} returns the largest \code{tau_k} whose predicted quantile
#'   is \code{<= y} (\code{tau_1} below the first quantile), and
#'   \code{qr_quantile} returns the predicted quantile at the largest
#'   \code{tau_k <= p} (the first quantile below \code{tau_1}).
#'
#' @param q_sorted matrix from \code{qr_sorted_quantiles()}.
#' @param taus quantile levels, \code{fit$tau}.
#' @param y,p numeric vectors with one value per row of \code{q_sorted}.
#'
#' @return numeric vector, one value per row of \code{q_sorted}.
#'
#' @keywords internal
qr_cdf <- function(q_sorted, taus, y) {
  vapply(seq_along(y), function(i) {
    taus[max(1L, findInterval(y[i], q_sorted[i, ]))]
  }, numeric(1))
}

#' @rdname qr_cdf
#' @keywords internal
qr_quantile <- function(q_sorted, taus, p) {
  vapply(seq_along(p), function(i) {
    q_sorted[i, max(1L, findInterval(p[i], taus))]
  }, numeric(1))
}
