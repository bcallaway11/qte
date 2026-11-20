# =============================================================================
# Title: Multi-cohort QTT aggregation tests
# Description: With several treated cohorts of unequal size and effects that
#   grow with event time, the dynamic QTT at each event time is a known
#   constant (every cell at that event time has the same location shift).
#   Aggregating with weights attached to the wrong cells (ptetools 8ede89a,
#   qte#19) mixes in pre-treatment cells or other event times and misses by
#   whole units. Also checks the ATT and QTT - ATT that come with QTT runs
#   (ptetools #27), including panel_qtt()'s own long-copula aggregation.
# Author: Brant Callaway
# Last update: 2026-11-19
# Date created: 2026-11-19
# =============================================================================

# never treated and cohorts 3, 4, 5 (unequal sizes), periods 1-6. Untreated
# outcomes: unit effect + period + noise, so every estimator's assumptions
# hold. Effect at event time e >= 0 is 1 + 3e, a location shift.
multicohort_data <- function(seed = 1) {
  set.seed(seed)
  g <- rep(c(0, 3, 4, 5), times = c(1500, 1500, 800, 400))
  d <- expand.grid(id = seq_along(g), period = 1:6)
  d$group <- g[d$id]
  e <- d$period - d$group
  treated <- d$group > 0 & e >= 0
  d$y <- rnorm(length(g))[d$id] + d$period + rnorm(nrow(d)) +
    ifelse(treated, 1 + 3 * e, 0)
  d
}

fit_multicohort <- function(f, ...) {
  set.seed(1)
  suppressMessages(f(
    yname = "y", gname = "group", tname = "period", idname = "id",
    data = multicohort_data(), gt_type = "qtt",
    probs = c(0.25, 0.5, 0.75), biters = 2, ...
  ))
}

test_that("dynamic QTTs and ATTs match the effect at each event time", {
  ests <- list(cic = cic, qdid = qdid, mdid = mdid, panel_qtt = panel_qtt,
               panel_qtt_short = function(...) {
                 panel_qtt(..., pre_copula = "short")
               })
  for (nm in names(ests)) {
    res <- fit_multicohort(ests[[nm]])
    dyn <- res$dynamic
    truth <- ifelse(dyn$e >= 0, 1 + 3 * dyn$e, 0)
    expect_equal(dyn$qtt, truth, tolerance = 0.25, label = paste(nm, "QTT"))

    att <- res$att$dynamic
    expect_false(is.null(att), label = paste(nm, "ATT returned"))
    expect_equal(att$att, ifelse(att$e >= 0, 1 + 3 * att$e, 0),
                 tolerance = 0.15, label = paste(nm, "ATT"))
    # location shifts: QTT(tau) = ATT at every tau
    expect_equal(res$qtt_minus_att$dynamic$qtt_minus_att,
                 rep(0, nrow(res$qtt_minus_att$dynamic)), tolerance = 0.25,
                 label = paste(nm, "QTT - ATT"))
  }
})

test_that("panel_qtt's long copula aggregates the ATT with the QTT weights", {
  res <- fit_multicohort(panel_qtt)
  # overall QTT - ATT is the overall QTT minus the overall ATT, same draws
  expect_equal(res$qtt_minus_att$overall$qtt_minus_att,
               res$overall$qtt - res$att$overall$att, tolerance = 1e-12)
  # only event times 0 and 1 have a long-copula window here; the overall
  # ATT averages cells at those event times, so it lies between 1 and 4
  expect_true(res$att$overall$att > 1 && res$att$overall$att < 4)
  expect_false(anyNA(res$att$overall$se))
})
