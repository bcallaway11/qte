# =============================================================================
# Title: Fast smoke tests — ddid() distributional DiD estimator
# Description: Regression checks for ddid() covering panel ATT and QTT modes.
#   Uses the same make_staggered_panel() helper as test-cic.R.  biters kept
#   low (20) so the suite runs quickly.  Reference values computed from the
#   first verified run.
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

# Re-uses the same helper as test-cic.R (loaded automatically from helper-*.R).
# If this file is run in isolation, define it here.
if (!exists("make_staggered_panel")) {
  make_staggered_panel <- function(seed = 42) {
    set.seed(seed)
    n <- 30L; periods <- 1:4; groups <- c(0, 2, 3)
    do.call(rbind, lapply(groups, function(g) {
      do.call(rbind, lapply(seq_len(n), function(i) {
        id <- (match(g, groups) - 1L) * n + i
        data.frame(id = id, period = periods, group = g,
                   y = rnorm(4L) + ifelse(g > 0 & periods >= g, 1, 0))
      }))
    }))
  }
}

# Panel ATT: pte_emp_boot class, 6 ATT(g,t) cells (2 groups x 3 post periods),
# no NAs, overall ATT near 1 (the true treatment effect in the synthetic data).
test_that("ddid panel ATT returns valid pte_emp_boot with stable estimate", {
  dat <- make_staggered_panel()
  set.seed(42)
  res <- ddid(yname = "y", gname = "group", tname = "period", idname = "id",
              data = dat, biters = 20)
  expect_s3_class(res, "pte_emp_boot")
  expect_true(is.numeric(res$overall_results$att))
  expect_false(anyNA(res$attgt_results$att))
  expect_equal(nrow(res$attgt_results), 6L)
  expect_equal(res$overall_results$att, 1.134847, tolerance = 1e-4)
})

# QTT mode: pte_qtt class, overall curve has one row per prob, all CIs valid,
# median QTT near 1.
test_that("ddid panel QTT returns valid pte_qtt with correct structure", {
  dat <- make_staggered_panel()
  set.seed(42)
  res <- ddid(yname = "y", gname = "group", tname = "period", idname = "id",
              data = dat, biters = 20, gt_type = "qtt", probs = c(0.25, 0.5, 0.75))
  expect_s3_class(res, "pte_qtt")
  expect_true(all(c("probs", "qtt", "se", "lower", "upper") %in% names(res$overall)))
  expect_equal(nrow(res$overall), 3L)
  expect_false(anyNA(res$overall$qtt))
  expect_true(all(res$overall$lower < res$overall$upper))
  expect_equal(res$overall$qtt[res$overall$probs == 0.5], 0.9470309, tolerance = 1e-4)
})
