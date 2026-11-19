# =============================================================================
# Title: Tests for the discrete-outcome check (#17)
# Description: outcome_discreteness() summarizes mass points and distinct
#   values (with weights); the estimators warn exactly once per call, with
#   class qte_discrete_outcome, when the outcome looks discrete, and not for
#   continuous or finely rounded outcomes or for ATTs that don't depend on
#   the outcome being continuous.
# Author: Brant Callaway
# Last update: 2026-11-19
# Date created: 2026-11-19
# =============================================================================

test_that("outcome_discreteness() finds the largest mass point, with weights", {
  d <- outcome_discreteness(c(0, 0, 0, 1, 2, NA))
  expect_equal(d$n, 5)
  expect_equal(d$n_unique, 3)
  expect_equal(d$max_mass, 0.6)
  expect_equal(d$mass_at, 0)

  # weights move the largest mass point
  dw <- outcome_discreteness(c(0, 0, 1), w = c(1, 1, 5))
  expect_equal(dw$max_mass, 5 / 7)
  expect_equal(dw$mass_at, 1)
})

test_that("binary, zero-inflated flagged; rounded, log counts not", {
  set.seed(1)
  n <- 5000
  flag <- function(y) {
    check_discrete_outcome(data.frame(y = y), "y", warn = FALSE)$discrete
  }
  expect_true(flag(rbinom(n, 1, 0.1)))
  expect_true(flag(ifelse(runif(n) < 0.5, 0, rlnorm(n))))
  expect_true(flag(sample(1:7, n, TRUE)))
  expect_false(flag(round(rnorm(n), 2)))
  expect_false(flag(log(rpois(n, 20) + 1)))
})

# small two-cohort panel with a binary outcome
binary_panel <- function() {
  set.seed(2)
  n <- 400
  g <- sample(c(0, 3), n, replace = TRUE)
  d <- expand.grid(id = seq_len(n), period = 1:3)
  d$g <- g[d$id]
  d$y <- as.numeric(rnorm(nrow(d)) + 0.3 * d$period > 1.2)
  d
}

count_discrete_warnings <- function(expr) {
  n <- 0L
  withCallingHandlers(
    expr,
    qte_discrete_outcome = function(w) {
      n <<- n + 1L
      invokeRestart("muffleWarning")
    },
    warning = function(w) invokeRestart("muffleWarning")
  )
  n
}

run_binary <- function(estimator, ...) {
  suppressMessages(estimator(
    yname = "y", gname = "g", tname = "period", idname = "id",
    data = binary_panel(), biters = 3, ...
  ))
}

test_that("estimators warn exactly once per call for a binary outcome", {
  # rank/quantile-based: ATT and QTT runs both warn
  expect_equal(count_discrete_warnings(run_binary(cic)), 1L)
  expect_equal(count_discrete_warnings(
    run_binary(cic, gt_type = "qtt", probs = c(0.25, 0.5, 0.75))
  ), 1L)
  # mdid's ATT is ordinary DiD: only quantile results warn
  expect_equal(count_discrete_warnings(run_binary(mdid)), 0L)
  expect_equal(count_discrete_warnings(
    run_binary(mdid, gt_type = "qtt", probs = c(0.25, 0.5, 0.75))
  ), 1L)
})

test_that("results carry the diagnostics; continuous outcomes don't warn", {
  res <- suppressWarnings(run_binary(cic))
  expect_true(res$outcome_discreteness$discrete)
  expect_equal(res$outcome_discreteness$n_unique, 2)

  data(mpdta, package = "did")
  n_warn <- count_discrete_warnings(suppressMessages(cic(
    yname = "lemp", gname = "first.treat", tname = "year",
    idname = "countyreal", data = mpdta, biters = 2
  )))
  expect_equal(n_warn, 0L)
})

test_that("unc_qte warns for lalonde's zero earnings and notes it in summary", {
  # 12% of re78 is exactly 0
  run <- function() {
    set.seed(1)
    unc_qte(yname = "re78", dname = "treat", data = lalonde.psid,
            probs = c(0.25, 0.5, 0.75), biters = 4)
  }
  expect_equal(count_discrete_warnings(run()), 1L)
  res <- suppressWarnings(run())
  expect_equal(res$outcome_discreteness$mass_at, 0)
  expect_output(print(summary(res)), "looks discrete")

  # turning the check off via the options
  op <- options(qte.mass_point_tol = 1, qte.discrete_n_unique = 0)
  on.exit(options(op))
  expect_equal(count_discrete_warnings(run()), 0L)
})
