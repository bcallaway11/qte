# =============================================================================
# Title: Quantile helpers
# Description: Shared, unexported quantile helpers: a weighted type-1
#   quantile (used by ddid_gt() and qdid_gt()), fast weighted ranks and
#   quantiles, reuse of identical quantile regressions across (g,t) cells,
#   and conditional CDF/quantile evaluation from quantile regressions with
#   rearrangement (used by the covariate-adjusted CiC, DDID and panel QTT
#   estimators).
# Author: Brant Callaway
# Last update: 2026-11-20
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
  y_s[first_geq(cdf, probs - tol)]
}

# --- Fast weighted ranks and quantiles ----------------------------------------
# These replace loops of the form vapply(x, function(v) <scan all of y>),
# which cost O(length(x) * length(y)), with one sort plus findInterval().
# They make the same comparisons, so they return the same values. The one
# exception is wcdf_at() with unequal weights, where the summation order
# differs and a value can differ in the last bit. Missing values fall back
# to the original loops, which define the NA behavior.

# Index of the first element of the non-decreasing vector `cdf` that is
# >= p, for each p; length(cdf) if there is none (or p is NA). Same as
# vapply(p, function(pp) { idx <- which(cdf >= pp); if (!length(idx))
# length(cdf) else idx[1] }, integer(1)).
first_geq <- function(cdf, p) {
  n <- length(cdf)
  if (anyNA(cdf)) {
    return(vapply(p, function(pp) {
      idx <- which(cdf >= pp)
      if (length(idx) == 0L) n else idx[1L]
    }, integer(1)))
  }
  # number of elements < p, plus one
  idx <- findInterval(p, cdf, left.open = TRUE) + 1L
  idx[is.na(idx) | idx > n] <- n
  idx
}

# Weighted CDF of `y` evaluated at each value of `x`:
# sum(w * (y <= x[i])) for each i.
wcdf_at <- function(y, w, x) {
  if (anyNA(y) || anyNA(w) || anyNA(x)) {
    return(vapply(x, function(v) sum(w * (y <= v)), numeric(1)))
  }
  o <- order(y)
  cw <- c(0, cumsum(w[o]))
  cw[findInterval(x, y[o]) + 1L]
}

# --- Reusing a quantile regression across (g,t) cells -------------------------
# For post-treatment cells, the comparison group in period t (G == 0 or
# G > t) does not depend on the treated cohort g, so the quantile regression
# of control outcomes in period t is identical for every cohort already
# treated by t. ptetools loops over t and then g, so those cells come one
# after another. cached_rq() keeps the most recent coefficients and returns
# them when the outcome, design matrix (kept columns, see cond_fit()),
# weights and quantiles are identical; otherwise it fits and stores the new
# ones. A reused fit is identical to refitting.
# options(qte.disable_qr_cache = TRUE) turns this off.
.qr_cache <- new.env(parent = emptyenv())

cached_rq <- function(y, xk, w, tau) {
  # row names differ between cells with the same data; drop them from the key
  rownames(xk) <- NULL
  y <- unname(y)
  w <- unname(w)
  use_cache <- !isTRUE(getOption("qte.disable_qr_cache"))
  key <- list(y = y, xk = xk, w = w, tau = tau)
  if (use_cache && !is.null(.qr_cache$fit) && identical(.qr_cache$key, key)) {
    return(.qr_cache$fit)
  }
  fit <- rq_coef(y, xk, w, tau)
  if (use_cache) {
    .qr_cache$key <- key
    .qr_cache$fit <- fit
  }
  fit
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
#' @description Evaluate the conditional CDF (\code{qr_cdf}) or quantile
#'   function (\code{qr_quantile}) for each row of \code{q_sorted} by linear
#'   interpolation between the estimated quantiles \code{(tau_k, Q(tau_k))}.
#'   Interpolating, rather than using a step function, avoids rounding ranks
#'   down to the quantile grid, which biases CiC-type counterfactuals.
#'
#'   Beyond the grid the functions are flat: \code{qr_cdf} returns
#'   \code{tau_1} below the first quantile and \code{tau_K} above the last,
#'   and \code{qr_quantile} returns the first (last) quantile for \code{p}
#'   below \code{tau_1} (above \code{tau_K}). In simulations, extrapolating
#'   with exponential tails reduced bias only slightly under normal errors
#'   and badly hurt heavy-tailed (t3) designs, and estimating extra extreme
#'   quantiles added more variance than it removed bias.
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
    q <- q_sorted[i, ]
    # degenerate row (all quantiles equal): approx() needs two distinct x's
    if (q[1] == q[length(q)]) return(if (y[i] < q[1]) taus[1] else taus[length(taus)])
    # ties (equal adjacent quantiles after rearrangement) get the mean tau
    stats::approx(q, taus, xout = y[i], rule = 2,
                  ties = list("ordered", mean))$y
  }, numeric(1))
}

#' @rdname qr_cdf
#' @keywords internal
qr_quantile <- function(q_sorted, taus, p) {
  vapply(seq_along(p), function(i) {
    stats::approx(taus, q_sorted[i, ], xout = p[i], rule = 2)$y
  }, numeric(1))
}
