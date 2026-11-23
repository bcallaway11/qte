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
# Last update: 2026-11-22
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
# `xlev`: factor levels for the covariates (see cond_xlevels()), so that the
# design has the same columns as for the units the fit is evaluated at.
#
# Returns an object of class "qte_dr": coef (n_thr x p, original covariate
# scale), grid, support, discrete flag, sample range, and the terms and
# factor levels needed to predict.
dr_fit <- function(formula, data, weighted = TRUE,
                   lambda = getOption("qte.dr_lambda", 1e-4), xlev = NULL) {
  mf <- stats::model.frame(formula, data, xlev = xlev)
  y <- stats::model.response(mf)
  tt <- stats::delete.response(stats::terms(mf))
  xmat <- stats::model.matrix(tt, mf)
  w <- if (!weighted || is.null(data$.w)) rep(1, nrow(xmat)) else data$.w
  w <- w / sum(w)
  g <- dr_grid(y, w)

  # Standardize covariates (column 1 is the intercept). Columns that are
  # collinear in this sample (constant, e.g. a dummy for a group with no
  # comparison units, or a combination of other columns) are not identified
  # and get a zero coefficient; check_cond_support() warns when that matters
  # for the units the fit is evaluated at.
  m <- colSums(xmat * w)
  s <- sqrt(colSums(w * sweep(xmat, 2, m)^2))
  keep <- seq_len(ncol(xmat)) %in% cond_keep(xmat)
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
                 discrete = g$discrete, terms = tt, xlev = xlev,
                 range = range(y)),
            class = "qte_dr")
}

# One-slot cache for the post-period control fit, which is identical across
# the cohorts already treated by that period (the same idea as cached_rq()).
# options(qte.disable_qr_cache = TRUE) turns it off too.
.dr_cache <- new.env(parent = emptyenv())

cached_dr <- function(formula, data, xlev = NULL) {
  key_data <- data[, c(all.vars(formula), ".w"), drop = FALSE]
  rownames(key_data) <- NULL
  key <- list(paste(deparse(formula), collapse = " "), key_data, xlev,
              getOption("qte.dr_n_grid", 99L), getOption("qte.dr_lambda", 1e-4))
  use_cache <- !isTRUE(getOption("qte.disable_qr_cache"))
  if (use_cache && !is.null(.dr_cache$fit) && identical(.dr_cache$key, key)) {
    return(.dr_cache$fit)
  }
  fit <- dr_fit(formula, data, xlev = xlev)
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
  xmat <- cond_design(fit$terms, newdata, fit$xlev)
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

# --- Design matrices and covariate support ----------------------------------
#
# A conditional model fit on one sample (say comparison units) and
# evaluated at another (treated units) is only identified at covariate
# values the fitting sample spans. Columns that are collinear in the
# fitting sample are dropped (a zero coefficient); that is harmless when the
# evaluation units satisfy the same linear relation (a duplicate column, a
# factor level no unit has) and an extrapolation otherwise (a dummy for a
# region with treated units only). check_cond_support() reports the second
# case, with ptetools::unsupported_covariates() (qte#16, ptetools#34).

# Factor levels of the covariates in `formula` over all rows of `data`
# (typically a whole (g,t) cell). Design matrices built for subsets with
# these levels all have the same columns, even when a subset lacks a level:
# that column is then 0 there, and dropped as collinear.
cond_xlevels <- function(formula, data) {
  tt <- stats::delete.response(stats::terms(formula))
  mf <- stats::model.frame(tt, data, na.action = stats::na.pass)
  stats::.getXlevels(tt, mf)
}

# Design matrix of terms `tt` (no response) at the rows of `data`.
cond_design <- function(tt, data, xlev = NULL) {
  mf <- stats::model.frame(tt, data, xlev = xlev, na.action = stats::na.pass)
  stats::model.matrix(tt, mf)
}

# Columns of a design matrix kept for a fit: those that are not collinear,
# by the pivoted QR rule of ptetools::pte_attgt() (and of
# ptetools::unsupported_covariates() with keep_covs = NULL), so the fit and
# the support check drop the same columns. With full rank, all columns.
cond_keep <- function(xmat) {
  qr0 <- qr(t(xmat) %*% xmat / nrow(xmat))
  sort(qr0$pivot[seq_len(qr0$rank)])
}

# Warn (once, for this (g,t) cell) when a conditional model fit on the
# samples in `fits` is evaluated at units of `target` that the fitting
# samples do not span: the dropped-column case described above.
#
# Args: formula (outcome model formula; its covariates are checked); fits
#   (list of data frames the models are fit on; with several, a unit is
#   flagged if any of them fails); target (data frame of the units the
#   models are evaluated at; shares are weighted by its .w column);
#   fitted_on, predicted_for (who these are, for the message).
# Returns: NULL invisibly; signals ptetools::pte_cell_warning(), which
#   ptetools::pte() reports once per call, listing the cells.
check_cond_support <- function(formula, fits, target,
                               fitted_on = "comparison units",
                               predicted_for = "treated units") {
  if (length(BMisc::rhs_vars(formula)) == 0L || nrow(target) == 0L) {
    return(invisible(NULL))
  }
  tt <- stats::delete.response(stats::terms(formula))
  # one set of factor levels for all the samples, so the columns line up
  xlev <- cond_xlevels(formula, do.call(rbind, c(fits, list(target))))
  x1 <- cond_design(tt, target, xlev)

  flagged <- list()
  for (f in fits) {
    x0 <- cond_design(tt, f, xlev)
    x0 <- x0[stats::complete.cases(x0), , drop = FALSE]
    if (nrow(x0) == 0L) next
    fl <- ptetools::unsupported_covariates(
      rbind(x0, x1), rep(c(0, 1), c(nrow(x0), nrow(x1)))
    )
    fl[is.na(fl)] <- FALSE
    for (nm in colnames(fl)) {
      flagged[[nm]] <- if (is.null(flagged[[nm]])) fl[, nm] else
        flagged[[nm]] | fl[, nm]
    }
  }
  if (length(flagged) == 0L) return(invisible(NULL))

  w <- if (is.null(target$.w)) rep(1, nrow(target)) else target$.w
  share <- vapply(flagged, function(f) sum(w[f]) / sum(w), numeric(1))
  ptetools::pte_cell_warning(
    paste0("Covariates dropped from the conditional outcome models because ",
           "they do not vary among ", fitted_on, " (given the other ",
           "covariates) but do among ", predicted_for, "; the estimates ",
           "extrapolate for those units (see ptetools::check_overlap())"),
    paste0(names(share), " (", round(100 * share, 1), "% of ",
           predicted_for, ")", collapse = ", ")
  )
}

# --- Conditional distributions -----------------------------------------------

# Fit a conditional distribution model of `formula` on `data`.
#   cond_dist "qr" (quantile regression at `taus`) or "dr".
#   cache     reuse the last identical fit (cached_rq() / cached_dr()); for
#             the post-period control fit shared across cohorts.
#   weighted  use the `.w` column as sampling weights.
#   xlev      factor levels (cond_xlevels() of the whole cell), so the fit
#             can be evaluated at units with levels this sample lacks.
# Collinear columns are dropped (cond_keep()); QR is fit on the kept
# columns, so a dummy with no variation no longer stops rq() with
# "Singular design matrix". With full rank the fit is unchanged.
cond_fit <- function(formula, data, cond_dist = c("qr", "dr"),
                     taus = seq(0.01, 0.99, 0.01), cache = FALSE,
                     weighted = TRUE, xlev = NULL) {
  cond_dist <- match.arg(cond_dist)
  if (cond_dist == "dr") {
    if (cache && weighted) return(cached_dr(formula, data, xlev))
    return(dr_fit(formula, data, weighted = weighted, xlev = xlev))
  }
  mf <- if (weighted) {
    stats::model.frame(formula, data, weights = .w, xlev = xlev)
  } else {
    stats::model.frame(formula, data, xlev = xlev)
  }
  y <- stats::model.response(mf)
  tt <- stats::delete.response(stats::terms(mf))
  xmat <- stats::model.matrix(tt, mf)
  keep <- cond_keep(xmat)
  xk <- xmat[, keep, drop = FALSE]
  w <- if (weighted) stats::model.weights(mf) else NULL
  coef <- if (cache && weighted) cached_rq(y, xk, w, taus) else
    rq_coef(y, xk, w, taus)
  structure(list(coef = coef, keep = keep, terms = tt, xlev = xlev,
                 taus = taus),
            class = "qte_qr")
}

# Coefficients (p x length(taus)) of a quantile regression of y on the
# design xk (which includes the intercept column), weighted by w unless
# NULL.
#
# Uses the Frisch-Newton interior point method (method = "fn"; Portnoy and
# Koenker 1997), not rq()'s default Barrodale-Roberts simplex ("br"). With
# a mass point in y (e.g. earnings with many zeros) and tau inside it, the
# linear program is degenerate and "br" can cycle forever; the loop is in
# Fortran, so it cannot be interrupted (this hung a bootstrap draw of
# unc_qte() on lalonde.psid). "fn" does not pivot and cannot cycle, is as
# fast or faster here, and qte uses none of the "br"-only features (rank-
# inversion inference, the full quantile process): only coef() at a fixed
# grid of taus, with bootstrap inference. Estimates agree with "br" up to
# the solver tolerance when the solution is unique.
rq_coef <- function(y, xk, w, taus) {
  fit <- if (is.null(w)) {
    suppressWarnings(quantreg::rq(y ~ xk - 1, tau = taus, method = "fn"))
  } else {
    suppressWarnings(quantreg::rq(y ~ xk - 1, tau = taus, weights = w,
                                  method = "fn"))
  }
  cf <- as.matrix(stats::coef(fit))
  if (length(taus) == 1L) cf <- matrix(cf, ncol = 1L)
  cf
}

# n x K matrix of predicted conditional quantiles at the rows of `newdata`,
# each row sorted (the rearrangement in qr_sorted_quantiles()).
qr_model_quantiles <- function(model, newdata) {
  x <- cond_design(model$terms, newdata, model$xlev)[, model$keep,
                                                     drop = FALSE]
  q <- x %*% model$coef
  matrix(t(apply(q, 1, sort)), ncol = ncol(q))
}

# Conditional CDF of unit i (row i of `newdata`) at y[i].
cond_cdf <- function(model, newdata, y) {
  if (inherits(model, "qte_dr")) return(dr_cdf(model, newdata, y))
  qr_cdf(qr_model_quantiles(model, newdata), model$taus, y)
}

# Conditional quantile of unit i (row i of `newdata`) at p[i].
cond_quantile <- function(model, newdata, p) {
  if (inherits(model, "qte_dr")) return(dr_quantile(model, newdata, p))
  qr_quantile(qr_model_quantiles(model, newdata), model$taus, p)
}

# n x K matrix of conditional quantiles at `taus` (rows sorted), for the
# estimators that pool predicted quantiles to integrate out the covariates.
cond_quantile_matrix <- function(model, newdata, taus = model$taus) {
  if (inherits(model, "qte_dr")) {
    return(dr_quantile_matrix(model, newdata, taus))
  }
  qr_model_quantiles(model, newdata)
}

# --- Conditional means (mdid) ------------------------------------------------

# Weighted least squares of `formula` on `data` with collinear columns
# dropped by the same rule (lm() would give them NA coefficients and
# predict() a rank-deficiency warning). Returns the coefficients and what
# cond_mean_predict() needs. With full rank, the same fit as lm().
cond_mean_fit <- function(formula, data, xlev = NULL) {
  mf <- stats::model.frame(formula, data, weights = .w, xlev = xlev)
  tt <- stats::delete.response(stats::terms(mf))
  xmat <- stats::model.matrix(tt, mf)
  keep <- cond_keep(xmat)
  fit <- stats::lm.wfit(xmat[, keep, drop = FALSE], stats::model.response(mf),
                        stats::model.weights(mf))
  list(coef = fit$coefficients, keep = keep, terms = tt, xlev = xlev)
}

cond_mean_predict <- function(model, newdata) {
  x <- cond_design(model$terms, newdata, model$xlev)[, model$keep,
                                                     drop = FALSE]
  drop(x %*% model$coef)
}
