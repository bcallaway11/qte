# =============================================================================
# Title: Covariates are taken from the base period
# Description: With panel data, every estimator conditions on each unit's
#   covariates in one pre-treatment period: the base period, or pre2 (the
#   earliest period used) for panel_qtt. Time-varying covariates must then
#   give exactly the estimates of the same data with the covariates replaced
#   by their reference-period values. Before 2.1.0, cic, qdid, mdid and
#   cdf_did fit the comparison group's post-period distribution on
#   post-period covariates, and lou_qtt used post-treatment covariates.
# Author: Brant Callaway
# Last update: 2026-11-20
# Date created: 2026-11-20
# =============================================================================

# panel with time-varying covariate x; `ref_period` is the period whose x
# every row gets in the "fixed" version
tv_covariate_data <- function(n_periods, g, ref_period, seed = 1) {
  set.seed(seed)
  n <- 1200
  G <- ifelse(runif(n) < 0.4, g, 0)
  x0 <- rnorm(n)
  d <- expand.grid(id = seq_len(n), period = seq_len(n_periods))
  d$G <- G[d$id]
  # x drifts over time, differently for each unit
  d$x <- x0[d$id] + d$period * rnorm(nrow(d), sd = 0.7)
  d$y <- x0[d$id] + 0.5 * d$period + rnorm(nrow(d)) +
    ifelse(d$G > 0 & d$period >= d$G, 1, 0)
  fixed <- d
  fixed$x <- d$x[d$period == ref_period][fixed$id]
  list(tv = d, fixed = fixed)
}

fit_tv <- function(f, data, ...) {
  set.seed(1)
  suppressWarnings(suppressMessages(f(
    yname = "y", gname = "G", tname = "period", idname = "id", data = data,
    xformula = ~x, gt_type = "qtt", probs = c(0.25, 0.5, 0.75), biters = 2,
    ...
  )))
}

expect_same_estimates <- function(f, d, label, ...) {
  a <- fit_tv(f, d$tv, ...)
  b <- fit_tv(f, d$fixed, ...)
  expect_equal(a$overall$qtt, b$overall$qtt, tolerance = 1e-10,
               label = paste(label, "QTT"))
  expect_equal(a$att$overall$att, b$att$overall$att, tolerance = 1e-10,
               label = paste(label, "ATT"))
}

test_that("two-period estimators use base-period covariates", {
  d <- tv_covariate_data(n_periods = 2, g = 2, ref_period = 1)
  # the covariate really does vary over time
  expect_false(isTRUE(all.equal(d$tv$x, d$fixed$x)))
  expect_same_estimates(cic, d, "cic")
  expect_same_estimates(cic, d, "cic dr", cond_dist = "dr")
  expect_same_estimates(qdid, d, "qdid")
  expect_same_estimates(mdid, d, "mdid")
  expect_same_estimates(ddid, d, "ddid")
  expect_same_estimates(cdf_did, d, "cdf_did")
  expect_same_estimates(lou_qtt, d, "lou_qtt")
  expect_same_estimates(lou_qtt, d, "lou_qtt reg", est_method = "reg")
})

test_that("panel_qtt uses covariates from pre2, the earliest period used", {
  # g = 3 with periods 1-3: pre2 = 1, pre1 = 2, post = 3 (long and short)
  d <- tv_covariate_data(n_periods = 3, g = 3, ref_period = 1)
  expect_same_estimates(panel_qtt, d, "panel_qtt")
  expect_same_estimates(panel_qtt, d, "panel_qtt short", pre_copula = "short")
})
