# =============================================================================
# Title: Fast smoke tests — panel_qtt() distributional panel QTT estimator
# Description: Regression checks for panel_qtt() covering panel ATT and QTT
#   modes. Uses a synthetic 4-period staggered panel (1 group treated at t=3)
#   so that pre2 = period 1 and pre1 = period 2 are always available.
#   biters kept low (20) so the suite runs quickly. Reference values computed
#   from the first verified run.
# Author: Brant Callaway
# Last update: 2026-05-20
# Date created: 2026-05-20
# =============================================================================

# Synthetic 4-period panel: never-treated (group 0) + one group treated at t=3.
# Periods 1..4 give pre2=1, pre1=2, post={3,4} for the treated group.
make_panel_qtt_data <- function(seed = 42) {
  set.seed(seed)
  n <- 30L; periods <- 1:4; groups <- c(0, 3)
  do.call(rbind, lapply(groups, function(g) {
    do.call(rbind, lapply(seq_len(n), function(i) {
      id <- (match(g, groups) - 1L) * n + i
      data.frame(id = id, period = periods, group = g,
                 y = rnorm(4L) + ifelse(g > 0 & periods >= g, 1, 0))
    }))
  }))
}

# Panel ATT (pre_copula = "short"): 2 ATT(g,t) cells (1 group x 2 post periods),
# no NAs, overall ATT near 1 (true treatment effect in the synthetic data).
# Using "short" here so both post-treatment cells are valid with a 4-period dataset.
test_that("panel_qtt panel ATT returns valid pte_emp_boot with stable estimate", {
  dat <- make_panel_qtt_data()
  set.seed(42)
  res <- panel_qtt(yname = "y", gname = "group", tname = "period", idname = "id",
                   data = dat, biters = 20, pre_copula = "short")
  expect_s3_class(res, "pte_emp_boot")
  expect_true(is.numeric(res$overall_results$att))
  expect_false(anyNA(res$attgt_results$att))
  expect_equal(nrow(res$attgt_results), 2L)
  expect_equal(res$overall_results$att, 1.057125, tolerance = 1e-4)
})

# QTT mode (pre_copula = "long", the default): only e=0 is valid for a 4-period
# dataset with group=3, so aggregation is over a single cell. pte_qtt class,
# overall curve has one row per prob, all CIs valid, median QTT near 1.
test_that("panel_qtt panel QTT returns valid pte_qtt with correct structure", {
  dat <- make_panel_qtt_data()
  set.seed(42)
  res <- panel_qtt(yname = "y", gname = "group", tname = "period", idname = "id",
                   data = dat, biters = 20, gt_type = "qtt", probs = c(0.25, 0.5, 0.75))
  expect_s3_class(res, "pte_qtt")
  expect_true(all(c("probs", "qtt", "se", "lower", "upper") %in% names(res$overall)))
  expect_equal(nrow(res$overall), 3L)
  expect_false(anyNA(res$overall$qtt))
  expect_true(all(res$overall$lower < res$overall$upper))
  expect_equal(res$overall$qtt[res$overall$probs == 0.5], 0.9732522, tolerance = 1e-4)
})
