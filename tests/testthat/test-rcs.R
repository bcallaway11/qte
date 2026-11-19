# =============================================================================
# Title: Tests for repeated cross sections (#18)
# Description: cic(), qdid(), mdid(), cdf_did() and lou_qtt() (with
#   lagged_outcome_cov = FALSE) accept panel = FALSE. Without covariates they
#   only use each period's outcome distributions, so on panel data the
#   repeated-cross-section estimates equal the panel ones, for the ATT and
#   the QTT. They also run on genuine repeated cross sections (new units
#   every period, no idname). ddid() and panel_qtt() need panel data and
#   have no panel argument; lou_qtt()'s lagged outcome and gt_type = "qott"
#   need panel data too, and say so.
# Author: Brant Callaway
# Last update: 2026-11-19
# Date created: 2026-11-19
# =============================================================================

# balanced panel: never treated, and cohorts 2 and 3; effect +1
rcs_panel_data <- function(seed = 42) {
  set.seed(seed)
  n <- 60
  g <- rep(c(0, 2, 3), each = n)
  d <- expand.grid(id = seq_along(g), period = 1:4)
  d$group <- g[d$id]
  d$y <- rnorm(length(g))[d$id] + 0.5 * d$period +
    (d$group > 0 & d$period >= d$group) + rnorm(nrow(d))
  d
}

# the same design, but a fresh sample of units every period
rcs_fresh_data <- function(seed = 7) {
  set.seed(seed)
  n <- 150
  d <- do.call(rbind, lapply(1:4, function(t) {
    data.frame(period = t, group = rep(c(0, 2, 3), each = n))
  }))
  d$y <- rnorm(nrow(d)) + 0.5 * d$period +
    (d$group > 0 & d$period >= d$group) + rnorm(nrow(d))
  d
}

rcs_estimators <- list(
  cic = list(f = cic, args = list()),
  qdid = list(f = qdid, args = list()),
  mdid = list(f = mdid, args = list()),
  cdf_did = list(f = cdf_did, args = list()),
  lou_qtt = list(f = lou_qtt, args = list(lagged_outcome_cov = FALSE))
)

fit_rcs <- function(e, data, ...) {
  set.seed(1)
  suppressMessages(do.call(e$f, c(
    list(yname = "y", gname = "group", tname = "period", data = data,
         biters = 5, ...),
    e$args
  )))
}

test_that("panel = FALSE gives the panel estimates on panel data", {
  d <- rcs_panel_data()
  for (nm in names(rcs_estimators)) {
    e <- rcs_estimators[[nm]]
    att_p <- fit_rcs(e, d, idname = "id")
    att_r <- fit_rcs(e, d, panel = FALSE)
    expect_equal(att_r$overall_results$att, att_p$overall_results$att,
                 tolerance = 1e-10, label = paste(nm, "ATT"))
    expect_equal(att_r$dyn_results$att.e, att_p$dyn_results$att.e,
                 tolerance = 1e-10, label = paste(nm, "event study"))

    qtt <- list(gt_type = "qtt", probs = c(0.25, 0.5, 0.75))
    qtt_p <- do.call(fit_rcs, c(list(e, d, idname = "id"), qtt))
    qtt_r <- do.call(fit_rcs, c(list(e, d, panel = FALSE), qtt))
    expect_equal(qtt_r$overall$qtt, qtt_p$overall$qtt, tolerance = 1e-10,
                 label = paste(nm, "QTT"))
  }
})

test_that("genuine repeated cross sections: ATT and QTT near the truth", {
  d <- rcs_fresh_data()
  for (nm in names(rcs_estimators)) {
    e <- rcs_estimators[[nm]]
    att <- fit_rcs(e, d, panel = FALSE)
    expect_false(anyNA(att$attgt_results$att), label = nm)
    expect_equal(att$overall_results$att, 1, tolerance = 0.3, label = nm)

    q <- fit_rcs(e, d, panel = FALSE, gt_type = "qtt",
                 probs = c(0.25, 0.5, 0.75))
    expect_s3_class(q, "pte_qtt")
    expect_false(anyNA(q$overall$se), label = nm)
    expect_equal(q$overall$qtt, rep(1, 3), tolerance = 0.5, label = nm)
  }
})

test_that("panel-only features error clearly with panel = FALSE", {
  d <- rcs_fresh_data()
  base <- list(yname = "y", gname = "group", tname = "period", data = d,
               panel = FALSE, biters = 2)
  # lou_qtt's default lagged outcome needs the same units in both periods
  expect_error(do.call(lou_qtt, base), "lagged_outcome_cov = FALSE")
  # the distribution of individual effects needs panel data
  for (f in list(cic, qdid, mdid)) {
    expect_error(do.call(f, c(base, gt_type = "qott")), "panel")
  }
  # ddid() and panel_qtt() have no repeated-cross-section version
  expect_error(do.call(ddid, base), "unused argument")
  expect_error(do.call(panel_qtt, base), "unused argument")
})
