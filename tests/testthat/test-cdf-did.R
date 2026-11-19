# =============================================================================
# Title: Tests for cdf_did() / cdf_did_gt()
# Description: Without covariates the ATT is the DiD ATT; the counterfactual
#   distribution and QTT are recovered when parallel trends in CDFs holds,
#   for a binary and a continuous outcome; with covariates, all three
#   methods (aipw, ipw, reg) recover the truth in a design where
#   only conditional parallel trends holds; the rearrangement produces a
#   valid CDF and reports how much it changed.
# Author: Brant Callaway
# Last update: 2026-11-19
# Date created: 2026-11-19
# =============================================================================

# one 2x2 cell from draws of the four (period x group) outcome distributions
make_cell <- function(trt_pre, trt_post, ctrl_pre, ctrl_post, x = NULL) {
  part <- function(y, name, D, xx) {
    out <- data.frame(name = name, D = D, Y = y, .w = 1)
    if (!is.null(xx)) out$x <- xx
    out
  }
  cell <- rbind(part(trt_pre$y, "pre", 1, trt_pre$x),
                part(trt_post$y, "post", 1, trt_post$x),
                part(ctrl_pre$y, "pre", 0, ctrl_pre$x),
                part(ctrl_post$y, "post", 0, ctrl_post$x))
  cell$id <- seq_len(nrow(cell))
  cell
}

# quantile() of BMisc distribution functions is exact as of BMisc 1.4.11
# (BMisc #10; before, it was wrong with ties or unequal weights)
qtt_from <- function(res, probs) {
  egr <- res$extra_gt_returns
  quantile(egr$F1, probs, type = 1) - quantile(egr$F0, probs, type = 1)
}

test_that("without covariates the ATT is the DiD ATT (panel and RCS)", {
  data(mpdta, package = "did")
  fit <- function(f, ...) {
    set.seed(1)
    suppressMessages(f(yname = "lemp", gname = "first.treat", tname = "year",
                       idname = "countyreal", data = mpdta, biters = 2, ...))
  }
  expect_equal(fit(cdf_did)$overall_results$att,
               fit(mdid)$overall_results$att, tolerance = 1e-10)
  expect_equal(fit(cdf_did, panel = FALSE)$overall_results$att,
               fit(mdid)$overall_results$att, tolerance = 1e-10)
})

test_that("binary outcome: counterfactual distribution and QTT recovered", {
  # P(Y = 1): treated 0.3 -> 0.4 untreated (0.6 observed, effect +0.2),
  # comparison 0.5 -> 0.6. So F0(0) = 0.7 + 0.4 - 0.5 = 0.6.
  set.seed(1)
  n <- 20000
  b <- function(p) list(y = rbinom(n, 1, p))
  res <- cdf_did_gt(make_cell(b(0.3), b(0.6), b(0.5), b(0.6)))
  expect_equal(res$attgt, 0.2, tolerance = 0.03)
  expect_equal(res$extra_gt_returns$F0(0), 0.6, tolerance = 0.02)
  expect_equal(unname(qtt_from(res, c(0.25, 0.5, 0.75))), c(0, 1, 0))
})

# Mixtures of N(0, 1) and N(3, 1): the CDF is linear in the mixing weight,
# so parallel trends in mixing weights is parallel trends in CDFs. Treated
# post-period outcomes are untreated ones + 1, so the QTT is 1 everywhere.
mix <- function(n, pi, x = NULL) {
  list(y = ifelse(runif(n) < pi, rnorm(n), rnorm(n, 3)), x = x)
}

test_that("continuous outcome: QTT recovered under CDF parallel trends", {
  set.seed(2)
  n <- 20000
  trt_post <- mix(n, 0.3)
  trt_post$y <- trt_post$y + 1
  res <- cdf_did_gt(make_cell(mix(n, 0.5), trt_post, mix(n, 0.6),
                              mix(n, 0.4)))
  expect_equal(res$attgt, 1, tolerance = 0.05)
  expect_equal(unname(qtt_from(res, c(0.25, 0.5, 0.75))), rep(1, 3),
               tolerance = 0.1)
})

test_that("covariates: aipw, ipw and reg recover the truth", {
  # binary x, more common among the treated (0.7 vs 0.3). Mixing weights
  # follow parallel trends within x, with a steeper trend at x = 1, so
  # unconditional parallel trends fails. Truth: ATT = 1, QTT = 1.
  set.seed(3)
  n <- 20000
  draw <- function(px1, pi0, pi1, shift = 0) {
    x <- rbinom(n, 1, px1)
    out <- mix(n, ifelse(x == 1, pi1, pi0), x)
    out$y <- out$y + shift
    out
  }
  cell <- make_cell(trt_pre = draw(0.7, 0.5, 0.5),
                    trt_post = draw(0.7, 0.4, 0.1, shift = 1),
                    ctrl_pre = draw(0.3, 0.6, 0.6),
                    ctrl_post = draw(0.3, 0.5, 0.2))

  unconditional <- cdf_did_gt(cell)
  expect_gt(abs(unconditional$attgt - 1), 0.2)

  for (em in c("aipw", "ipw", "reg")) {
    res <- cdf_did_gt(cell, xformula = ~x, est_method = em)
    expect_equal(res$attgt, 1, tolerance = 0.06, label = em)
    expect_equal(unname(qtt_from(res, c(0.25, 0.5, 0.75))), rep(1, 3),
                 tolerance = 0.15, label = em)
  }
})

test_that("rearrangement gives a valid CDF and reports the change", {
  # treated pre all 0, comparison 1 -> 0: raw F0 = (2, 1) at (0, 1)
  one <- function(v) list(y = rep(v, 50))
  res <- cdf_did_gt(make_cell(one(0), one(0), one(1), one(0)))
  F0 <- res$extra_gt_returns$F0
  expect_equal(F0(c(-1, 0, 1)), c(0, 1, 1))
  expect_equal(res$extra_gt_returns$cdf_violation, 1)

  # identical distributions everywhere: nothing to rearrange
  set.seed(4)
  same <- function() list(y = rnorm(200))
  res2 <- cdf_did_gt(make_cell(same(), same(), same(), same()))
  expect_gte(res2$extra_gt_returns$cdf_violation, 0)
  expect_true(all(diff(res2$extra_gt_returns$F0(seq(-3, 3, 0.1))) >= 0))
})

test_that("cdf_did QTT runs; only QTTs warn for discrete outcomes", {
  data(mpdta, package = "did")
  set.seed(1)
  q <- suppressMessages(cdf_did(
    yname = "lemp", gname = "first.treat", tname = "year",
    idname = "countyreal", data = mpdta, gt_type = "qtt",
    probs = c(0.25, 0.5, 0.75), biters = 5
  ))
  expect_s3_class(q, "pte_qtt")
  expect_false(anyNA(q$overall$qtt))

  d <- mpdta
  d$y <- as.numeric(d$lemp > median(d$lemp))
  run <- function(type) {
    suppressMessages(cdf_did(yname = "y", gname = "first.treat",
                             tname = "year", idname = "countyreal", data = d,
                             gt_type = type, probs = 0.5, biters = 2))
  }
  expect_no_warning(run("att"))
  expect_warning(run("qtt"), class = "qte_discrete_outcome")
})
