# =============================================================================
# Title: Fast smoke tests — cic() staggered adoption estimator
# Description: Regression checks for cic() covering panel ATT, repeated cross
#   sections ATT, and panel QTT. Uses biters=20 to keep runtime low.
#   Reference values were computed from the first verified run.
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

# Minimal balanced panel: 3 groups (never-treated, treated from t=2, from t=3),
# 30 units each, 4 periods. Outcome = noise + 1 unit treatment effect post-adoption.
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

# Panel + staggered adoption: pte_emp_boot, 6 ATT(g,t) cells (2 groups x 3 post
# periods), overall ATT near 1 (the true treatment effect in the synthetic data).
test_that("cic panel ATT returns valid pte_emp_boot with stable estimate", {
  dat <- make_staggered_panel()
  set.seed(42)
  res <- cic(yname = "y", gname = "group", tname = "period", idname = "id",
             data = dat, biters = 20)
  expect_s3_class(res, "pte_emp_boot")
  expect_true(is.numeric(res$overall_results$att))
  expect_false(anyNA(res$attgt_results$att))
  expect_equal(nrow(res$attgt_results), 6L)
  expect_equal(res$overall_results$att, 1.208897, tolerance = 1e-4)
})

# Repeated cross sections: same data with panel=FALSE; observations treated as
# independent within each period. Point estimate matches panel (same data).
test_that("cic RCS ATT returns valid pte_emp_boot with stable estimate", {
  dat <- make_staggered_panel()
  set.seed(42)
  res <- cic(yname = "y", gname = "group", tname = "period",
             data = dat, panel = FALSE, biters = 20)
  expect_s3_class(res, "pte_emp_boot")
  expect_true(is.numeric(res$overall_results$att))
  expect_false(anyNA(res$attgt_results$att))
  expect_equal(res$overall_results$att, 1.208897, tolerance = 1e-4)
})

# QTT mode: pte_qtt class, overall curve has one row per prob, all CIs valid,
# median QTT near 1. Checks that the CDF-mixing aggregation ran end-to-end.
test_that("cic panel QTT returns valid pte_qtt with correct structure", {
  dat <- make_staggered_panel()
  set.seed(42)
  res <- cic(yname = "y", gname = "group", tname = "period", idname = "id",
             data = dat, biters = 20, gt_type = "qtt", probs = c(0.25, 0.5, 0.75))
  expect_s3_class(res, "pte_qtt")
  expect_true(all(c("probs", "qtt", "se", "lower_pw", "upper_pw") %in% names(res$overall)))
  expect_equal(nrow(res$overall), 3L)
  expect_false(anyNA(res$overall$qtt))
  expect_true(all(res$overall$lower_pw < res$overall$upper_pw))
  expect_equal(res$overall$qtt[res$overall$probs == 0.5], 0.8210850, tolerance = 1e-4)
})
