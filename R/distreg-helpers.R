# =============================================================================
# Title: Distribution regression helpers
# Description: Conditional distributions by distribution regression (DR): one
#   logit of 1{Y <= y_j} on the covariates for each threshold y_j on a grid
#   (Foresi and Peracchi 1995; Chernozhukov, Fernandez-Val and Melly 2013).
#   This is the cond_dist = "dr" alternative to the quantile regression helpers
#   in quantile-helpers.R, with the same interface: a per-unit conditional CDF
#   at a value (dr_cdf, cf. qr_cdf) and a per-unit conditional quantile at a
#   probability (dr_quantile, cf. qr_quantile).
# Author: Brant Callaway
# Last update: 2026-11-19
# Date created: 2026-11-15
# =============================================================================

# --- Threshold grid ----------------------------------------------------------

# Thresholds for the DR fits. Continuous outcomes: the weighted percentiles
# 1%, ..., 99% of the fitting sample (the same grid as the QR estimators'
# u_seq), deduplicated, so a mass point (e.g. many zeros) becomes a single
# threshold at that value. Discrete outcomes (at most n_grid distinct
# values): every support point except the largest, at which the CDF is 1.
#
# Args: y outcome vector; w weights; n_grid number of percentiles.
# Returns: list(grid, support, discrete).
dr_grid <- function(y, w, n_grid = getOption("qte.dr_n_grid", 99L)) {
  support <- sort(unique(y))
  if (length(support) <= n_grid + 1L) {
    return(list(grid = support[-length(support)], support = support,
                discrete = TRUE))
  }
  probs <- seq_len(n_grid) / (n_grid + 1)
  grid <- unique(wquant(y, w, probs))
  # P(Y <= max) = 1 trivially; a logit there has no zeros to fit
  list(grid = grid[grid < max(y)], support = NULL, discrete = FALSE)
}

# --- Penalized weighted logit ------------------------------------------------

# Newton's method for a weighted logit with a small ridge penalty on the
# non-intercept coefficients, on standardized covariates. The penalty keeps
# the coefficients finite when a threshold separates the data (common in the
# tails with sparse dummies) at a negligible cost in bias; with no
# covariates the fit is the exact weighted share. Step halving guards
# against overshooting.
#
# Args: xs n x p standardized design (column 1 = intercept); z 0/1 response;
#   w weights summing to 1; b starting values; lambda penalty; maxit, tol.
# Returns: coefficient vector (standardized scale).
dr_logit <- function(xs, z, w, b, lambda, maxit = 50L, tol = 1e-10) {
  p <- ncol(xs)
  pen <- c(0, rep(lambda, p - 1L))
  obj <- function(bb) {
    eta <- drop(xs %*% bb)
    sum(w * (z * eta - log1p(exp(-abs(eta))) - pmax(eta, 0))) -
      0.5 * sum(pen * bb^2)
  }
  f_old <- obj(b)
  for (it in seq_len(maxit)) {
    eta <- drop(xs %*% b)
    mu <- stats::plogis(eta)
    grad <- drop(crossprod(xs, w * (z - mu))) - pen * b
    hess <- crossprod(xs, xs * (w * mu * (1 - mu)))
    diag(hess) <- diag(hess) + pen
    step <- tryCatch(solve(hess, grad),
                     error = function(e) qr.solve(hess + diag(1e-8, p), grad))
    # step halving: accept the first step that doesn't lower the objective
    s <- 1
    repeat {
      b_new <- b + s * step
      f_new <- obj(b_new)
      if (f_new >= f_old - 1e-14 || s < 1e-4) break
      s <- s / 2
    }
    converged <- max(abs(b_new - b)) < tol || abs(f_new - f_old) < tol
    b <- b_new
    f_old <- f_new
    if (converged) break
  }
  b
}

# --- Fit ----------------------------------------------------------------------

# Fit DR for `formula` (response Y, covariates on the right) on `data`, with
# sampling weights from the `.w` column when `weighted` (as the QR fits use).
# One penalized logit per threshold, warm-started from the previous one.
#
# Returns an object of class "qte_dr": coef (n_thr x p, original covariate
# scale), grid, support, discrete flag, sample range, and the terms needed to
# predict.
dr_fit <- function(formula, data, weighted = TRUE,
                   lambda = getOption("qte.dr_lambda", 1e-4)) {
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  tt <- stats::delete.response(stats::terms(mf))
  xmat <- stats::model.matrix(tt, mf)
  w <- if (!weighted || is.null(data$.w)) rep(1, nrow(xmat)) else data$.w
  w <- w / sum(w)
  g <- dr_grid(y, w)

  # Standardize covariates (column 1 is the intercept). Columns with no
  # variation in this sample (e.g. a dummy for a group with no comparison
  # units) carry no information and get a zero coefficient.
  m <- colSums(xmat * w)
  s <- sqrt(colSums(w * sweep(xmat, 2, m)^2))
  keep <- c(TRUE, s[-1] > 1e-12)
  xs <- sweep(sweep(xmat[, keep, drop = FALSE], 2, c(0, m[keep][-1])), 2,
              c(1, s[keep][-1]), "/")

  n_thr <- length(g$grid)
  coef_mat <- matrix(0, n_thr, ncol(xmat), dimnames = list(NULL, colnames(xmat)))
  b <- c(stats::qlogis(min(max(sum(w * (y <= g$grid[1])), 1e-6), 1 - 1e-6)),
         rep(0, sum(keep) - 1L))
  for (j in seq_len(n_thr)) {
    b <- dr_logit(xs, as.numeric(y <= g$grid[j]), w, b, lambda)
    # back to the original scale: eta = b0 + sum b_k (x_k - m_k) / s_k
    beta <- b / c(1, s[keep][-1])
    beta[1] <- b[1] - sum(beta[-1] * m[keep][-1])
    coef_mat[j, keep] <- beta
  }
  structure(list(coef = coef_mat, grid = g$grid, support = g$support,
                 discrete = g$discrete, terms = tt,
                 range = range(y)),
            class = "qte_dr")
}

# One-slot cache for the post-period control fit, which is identical across
# the cohorts already treated by that period (the same idea as cached_rq()).
# options(qte.disable_qr_cache = TRUE) turns it off too.
.dr_cache <- new.env(parent = emptyenv())

cached_dr <- function(formula, data) {
  key_data <- data[, c(all.vars(formula), ".w"), drop = FALSE]
  rownames(key_data) <- NULL
  key <- list(paste(deparse(formula), collapse = " "), key_data,
              getOption("qte.dr_n_grid", 99L), getOption("qte.dr_lambda", 1e-4))
  use_cache <- !isTRUE(getOption("qte.disable_qr_cache"))
  if (use_cache && !is.null(.dr_cache$fit) && identical(.dr_cache$key, key)) {
    return(.dr_cache$fit)
  }
  fit <- dr_fit(formula, data)
  if (use_cache) {
    .dr_cache$key <- key
    .dr_cache$fit <- fit
  }
  fit
}

# --- Predict ------------------------------------------------------------------

# n x n_thr matrix of conditional CDFs F(grid_j | x_i) at the rows of `newdata`,
# with each row sorted: the rearrangement of Chernozhukov, Fernandez-Val and
# Galichon (2010), which makes each CDF monotone when the logits cross.
dr_cdf_matrix <- function(fit, newdata) {
  xmat <- stats::model.matrix(fit$terms,
                              stats::model.frame(fit$terms, newdata))
  cdf_mat <- stats::plogis(xmat %*% t(fit$coef))
  if (ncol(cdf_mat) > 1) cdf_mat <- matrix(t(apply(cdf_mat, 1, sort)), nrow = nrow(cdf_mat))
  cdf_mat
}

# Tail anchors for continuous outcomes. The thresholds are *marginal*
# percentiles of Y, so for units with extreme covariates much of the
# conditional distribution lies outside [grid_1, grid_J]; a flat tail (the QR
# convention, which is right there because QR's grid is in probability space)
# would cap their ranks well inside (0, 1). Instead each conditional CDF is
# anchored at 0 at the sample minimum and at 1 at the sample maximum and
# interpolated linearly in between. An anchor is skipped when it coincides
# with a threshold (a mass point at the minimum, e.g. many zeros, is then a
# jump: F = 0 below it).
#
# Returns list(grid, cdf_mat) with the anchors added as extra columns.
dr_anchor <- function(fit, cdf_mat) {
  g <- fit$grid
  lo <- fit$range[1] < g[1]
  hi <- fit$range[2] > g[length(g)]
  if (lo) {
    g <- c(fit$range[1], g)
    cdf_mat <- cbind(0, cdf_mat)
  }
  if (hi) {
    g <- c(g, fit$range[2])
    cdf_mat <- cbind(cdf_mat, 1)
  }
  list(grid = g, cdf_mat = cdf_mat, lo = lo)
}

# Conditional CDF of unit i at y[i]. Continuous outcomes: linear
# interpolation between the anchored thresholds (see dr_anchor()); 0 below
# the sample minimum, 1 above the maximum. Discrete outcomes: the
# step-function CDF at the support points.
dr_cdf <- function(fit, newdata, y) {
  dr_cdf_rows(fit, dr_cdf_matrix(fit, newdata), y)
}

# Weighted average over the rows of `newdata` of the conditional CDF,
# evaluated at each point of y: sum_i w_i F(y | x_i). Averaging the
# threshold rows first and interpolating once is exact (both are linear) and
# avoids an n x length(y) matrix. Used by cdf_did_gt().
#
# Args: fit (qte_dr), newdata, w (weights summing to 1), y (points).
dr_avg_cdf <- function(fit, newdata, w, y) {
  avg <- colSums(dr_cdf_matrix(fit, newdata) * w)
  dr_cdf_rows(fit, matrix(avg, nrow = length(y), ncol = length(avg),
                          byrow = TRUE), y)
}

# dr_cdf() given the threshold matrix: row i of cdf_mat (from
# dr_cdf_matrix() or an average of its rows) evaluated at y[i].
dr_cdf_rows <- function(fit, cdf_mat, y) {
  n <- nrow(cdf_mat)
  n_thr <- ncol(cdf_mat)
  g <- fit$grid
  k <- findInterval(y, g)  # number of thresholds <= y
  if (fit$discrete) {
    # y below the smallest support point: 0; at or above the largest: 1
    out <- numeric(n)
    mid <- k >= 1 & y < max(fit$support)
    out[mid] <- cdf_mat[cbind(which(mid), k[mid])]
    out[y >= max(fit$support)] <- 1
    return(out)
  }
  a <- dr_anchor(fit, cdf_mat)
  g <- a$grid
  cdf_mat <- a$cdf_mat
  n_thr <- ncol(cdf_mat)
  k <- findInterval(y, g)
  out <- numeric(n)
  # below the first point: 0 if it is the sample-minimum anchor, else (mass
  # point at the minimum) 0 as well, since the CDF jumps there
  hi <- k == n_thr
  mid <- k >= 1 & !hi
  out[hi] <- cdf_mat[hi, n_thr]
  if (any(mid)) {
    i <- which(mid)
    kk <- k[mid]
    f0 <- cdf_mat[cbind(i, kk)]
    f1 <- cdf_mat[cbind(i, kk + 1L)]
    out[mid] <- f0 + (f1 - f0) * (y[mid] - g[kk]) / (g[kk + 1L] - g[kk])
  }
  out
}

# Invert sorted conditional CDFs: for each row i of `cdf_mat` (from
# dr_cdf_matrix()), the quantile at p[i]. Continuous: linear interpolation
# between the anchored thresholds (see dr_anchor()). Discrete: the smallest
# support point whose CDF is at least p.
dr_invert <- function(fit, cdf_mat, p) {
  n <- nrow(cdf_mat)
  if (fit$discrete) {
    k <- rowSums(cdf_mat < p)  # thresholds with CDF below p
    return(fit$support[k + 1L])
  }
  a <- dr_anchor(fit, cdf_mat)
  g <- a$grid
  cdf_mat <- a$cdf_mat
  n_thr <- ncol(cdf_mat)
  k <- rowSums(cdf_mat <= p)  # rows are sorted, so F[i, k] <= p < F[i, k + 1]
  out <- numeric(n)
  lo <- k == 0           # p below F at the first point (mass point at min)
  hi <- k == n_thr
  mid <- !lo & !hi
  out[lo] <- g[1]
  out[hi] <- g[n_thr]
  if (any(mid)) {
    i <- which(mid)
    kk <- k[mid]
    f0 <- cdf_mat[cbind(i, kk)]
    f1 <- cdf_mat[cbind(i, kk + 1L)]
    out[mid] <- g[kk] + (p[mid] - f0) / (f1 - f0) * (g[kk + 1L] - g[kk])
  }
  out
}

# Conditional quantile of unit i at p[i]: the inverse of dr_cdf().
dr_quantile <- function(fit, newdata, p) {
  dr_invert(fit, dr_cdf_matrix(fit, newdata), p)
}

# n x K matrix of conditional quantiles at the probabilities `taus` (one row
# per row of `newdata`), the DR counterpart of predicting a quantile
# regression at a vector of tau. The CDF matrix is computed once.
dr_quantile_matrix <- function(fit, newdata, taus) {
  cdf_mat <- dr_cdf_matrix(fit, newdata)
  n <- nrow(cdf_mat)
  vapply(taus, function(tau) dr_invert(fit, cdf_mat, rep(tau, n)), numeric(n))
}

# =============================================================================
# Dispatch: one interface for cond_dist = "qr" and "dr"
# =============================================================================

# The covariate-adjusted estimators call these instead of quantreg/DR
# directly, so each estimator branches on cond_dist in one place. The "qr"
# paths reproduce the previous code exactly (same rq() calls, rearrangement
# and interpolation).

# Fit a conditional distribution model of `formula` on `data`.
#   cond_dist "qr" (quantile regression at `taus`) or "dr".
#   cache     reuse the last identical fit (cached_rq() / cached_dr()); for
#             the post-period control fit shared across cohorts.
#   weighted  use the `.w` column as sampling weights.
cond_fit <- function(formula, data, cond_dist = c("qr", "dr"),
                     taus = seq(0.01, 0.99, 0.01), cache = FALSE,
                     weighted = TRUE) {
  cond_dist <- match.arg(cond_dist)
  if (cond_dist == "dr") {
    if (cache && weighted) return(cached_dr(formula, data))
    return(dr_fit(formula, data, weighted = weighted))
  }
  fit <- if (cache && weighted) {
    cached_rq(formula, data, taus)
  } else if (weighted) {
    suppressWarnings(quantreg::rq(formula, data = data, tau = taus,
                                  weights = .w))
  } else {
    suppressWarnings(quantreg::rq(formula, tau = taus, data = data))
  }
  structure(list(fit = fit, taus = taus), class = "qte_qr")
}

# Conditional CDF of unit i (row i of `newdata`) at y[i].
cond_cdf <- function(model, newdata, y) {
  if (inherits(model, "qte_dr")) return(dr_cdf(model, newdata, y))
  qr_cdf(qr_sorted_quantiles(model$fit, newdata), model$taus, y)
}

# Conditional quantile of unit i (row i of `newdata`) at p[i].
cond_quantile <- function(model, newdata, p) {
  if (inherits(model, "qte_dr")) return(dr_quantile(model, newdata, p))
  qr_quantile(qr_sorted_quantiles(model$fit, newdata), model$taus, p)
}

# n x K matrix of conditional quantiles at `taus` (rows sorted), for the
# estimators that pool predicted quantiles to integrate out the covariates.
cond_quantile_matrix <- function(model, newdata, taus = model$taus) {
  if (inherits(model, "qte_dr")) {
    return(dr_quantile_matrix(model, newdata, taus))
  }
  t(apply(as.matrix(stats::predict(model$fit, newdata = newdata)), 1, sort))
}
