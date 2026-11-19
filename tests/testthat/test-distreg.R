# =============================================================================
# Title: Tests for the distribution regression (cond_dist = "dr") helpers
# Description: Exact properties of dr_fit()/dr_cdf()/dr_quantile() (no
#   covariates and saturated covariates reproduce weighted ECDFs; quantile
#   inverts CDF; discrete outcomes stay on the support; tail anchors), and
#   that cic_gt(cond_dist = "dr") recovers a known ATT.
# Author: Brant Callaway
# Last update: 2026-11-15
# Date created: 2026-11-15
# =============================================================================

sim_cdf_data <- function(n = 2000, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = rbinom(n, 1, 0.4), .w = runif(n, 0.5, 2))
  d$Y <- 1 + d$x + (1 + d$x) * rnorm(n)
  d
}

wecdf_at <- function(y, w, grid) {
  vapply(grid, function(c) sum(w * (y <= c)) / sum(w), numeric(1))
}

test_that("without covariates DR reproduces the weighted ECDF at the thresholds", {
  d <- sim_cdf_data()
  fit <- dr_fit(Y ~ 1, d)
  expect_false(fit$discrete)
  expect_equal(drop(dr_cdf_matrix(fit, d[1, , drop = FALSE])),
               wecdf_at(d$Y, d$.w, fit$grid), tolerance = 1e-8)
})

test_that("with a saturated discrete covariate DR reproduces cell ECDFs", {
  d <- sim_cdf_data()
  op <- options(qte.dr_lambda = 0)
  fit <- dr_fit(Y ~ x, d)
  options(op)
  Fm <- dr_cdf_matrix(fit, data.frame(x = c(0, 1)))
  for (v in 0:1) {
    cell <- d[d$x == v, ]
    expect_equal(Fm[v + 1, ], wecdf_at(cell$Y, cell$.w, fit$grid),
                 tolerance = 1e-6)
  }
})

test_that("conditional CDFs are monotone, in [0, 1], and inverted by dr_quantile", {
  d <- sim_cdf_data()
  fit <- dr_fit(Y ~ x, d)
  Fm <- dr_cdf_matrix(fit, d[1:100, ])
  expect_true(all(apply(Fm, 1, function(r) all(diff(r) >= 0))))
  nd <- d[1:300, ]
  yy <- runif(300, min(d$Y), max(d$Y))
  u <- dr_cdf(fit, nd, yy)
  expect_true(all(u >= 0 & u <= 1))
  expect_equal(dr_quantile(fit, nd, u), yy, tolerance = 1e-8)
  # tail anchors: 0 below the sample minimum, 1 at and above the maximum
  expect_equal(dr_cdf(fit, nd[1:2, ], c(min(d$Y) - 1, max(d$Y) + 1)), c(0, 1))
})

test_that("discrete outcomes give step CDFs on the support", {
  set.seed(2)
  d <- data.frame(x = rnorm(1000), .w = 1)
  d$Y <- rpois(1000, exp(0.5 + 0.3 * d$x))
  fit <- dr_fit(Y ~ x, d)
  expect_true(fit$discrete)
  q <- dr_quantile(fit, d[1:50, ], runif(50))
  expect_true(all(q %in% fit$support))
  # CDF is constant between support points
  s <- fit$support[2]
  expect_equal(dr_cdf(fit, d[1:5, ], rep(s, 5)), dr_cdf(fit, d[1:5, ], rep(s + 0.5, 5)))
})

test_that("cic_gt with cond_dist = 'dr' recovers a known ATT", {
  # rank depends on D and x and is stable over time; effect = +1
  set.seed(3)
  n <- 3000
  x <- rnorm(n)
  D <- rbinom(n, 1, plogis(0.5 * x))
  u <- rnorm(n, 0.4 * D + 0.6 * x)
  gt <- rbind(
    data.frame(id = 1:n, D = D, name = "pre",  Y = u + x, x = x, .w = 1),
    data.frame(id = 1:n, D = D, name = "post", Y = 1 + 1.5 * u + x + D, x = x, .w = 1)
  )
  op <- options(qte.disable_qr_cache = TRUE)
  att_dr <- suppressWarnings(cic_gt(gt, ~x, cond_dist = "dr"))$attgt
  options(op)
  expect_equal(att_dr, 1, tolerance = 0.05)
  expect_error(cic_gt(gt, ~x, cond_dist = "nope"))
})

test_that("the cond_fit/cond_cdf/cond_quantile dispatch reproduces the QR helpers", {
  d <- sim_cdf_data()
  taus <- seq(0.01, 0.99, 0.01)
  m <- cond_fit(Y ~ x, d, "qr", taus)
  direct <- suppressWarnings(quantreg::rq(Y ~ x, data = d, tau = taus, weights = .w))
  q <- qr_sorted_quantiles(direct, d)
  expect_identical(cond_cdf(m, d, d$Y), qr_cdf(q, taus, d$Y))
  p <- runif(nrow(d))
  expect_identical(cond_quantile(m, d, p), qr_quantile(q, taus, p))
})

test_that("qdid_gt with cond_dist = 'dr' recovers a known ATT", {
  # time trend depends on x (unconditional PT fails); effect = +1
  set.seed(4)
  n <- 2000
  x <- rnorm(n)
  D <- rbinom(n, 1, plogis(0.5 * x))
  u <- rnorm(n, 0.4 * D + 0.6 * x)
  gt <- rbind(
    data.frame(id = 1:n, D = D, name = "pre",  Y = x + u, x = x, .w = 1),
    data.frame(id = 1:n, D = D, name = "post", Y = 0.5 + 1.2 * x + u + D, x = x, .w = 1)
  )
  op <- options(qte.disable_qr_cache = TRUE)
  expect_equal(suppressWarnings(qdid_gt(gt, ~x, cond_dist = "dr"))$attgt, 1,
               tolerance = 0.02)
  options(op)
})

test_that("cond_dist = 'dr' runs end to end in every covariate estimator", {
  data(mpdta, package = "did")
  small <- mpdta[mpdta$countyreal %in% unique(mpdta$countyreal)[1:200], ]
  args <- list(yname = "lemp", gname = "first.treat", tname = "year",
               idname = "countyreal", data = small, xformula = ~lpop,
               cond_dist = "dr", biters = 2, cl = 1, cband = FALSE)
  att_of <- function(r) r$overall_results$att
  set.seed(1)
  for (f in list(cic, qdid, ddid, panel_qtt)) {
    expect_true(is.finite(att_of(suppressWarnings(do.call(f, args)))))
  }
  expect_true(is.finite(att_of(suppressWarnings(
    do.call(lou_qtt, c(args, est_method = "reg"))))))
  data(lalonde, package = "qte")
  u <- suppressWarnings(unc_qte(yname = "re78", dname = "treat",
                                data = lalonde.psid, xformula = ~ age + re75,
                                est_method = "reg", target = "qtt",
                                cond_dist = "dr", biters = 2,
                                probs = c(0.25, 0.5, 0.75)))
  expect_true(all(is.finite(u$qte)))
})
