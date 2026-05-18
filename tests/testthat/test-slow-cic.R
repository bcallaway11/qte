# =============================================================================
# Title: Slow smoke tests — cic() staggered adoption estimator
# Description: Bootstrap SE checks for cic() at biters=100. Run before
#   releases with: R_SLOW_TESTS=true Rscript -e "devtools::test()"
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

# Shared helper — mirrors the one in test-cic.R
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

# Panel ATT with full bootstrap: SE is positive and CIs are non-degenerate.
test_that("cic panel ATT SE path produces valid intervals", {
  skip_slow()
  dat <- make_staggered_panel()
  set.seed(42)
  res <- cic(yname = "y", gname = "group", tname = "period", idname = "id",
             data = dat, biters = 100)
  expect_s3_class(res, "pte_emp_boot")
  expect_true(res$overall_results$se > 0)
  expect_true(res$overall_results$att - res$overall_results$se <
                res$overall_results$att)
})

# RCS ATT SE: same check for the repeated cross sections path.
test_that("cic RCS ATT SE path produces valid intervals", {
  skip_slow()
  dat <- make_staggered_panel()
  set.seed(42)
  res <- cic(yname = "y", gname = "group", tname = "period",
             data = dat, panel = FALSE, biters = 100)
  expect_s3_class(res, "pte_emp_boot")
  expect_true(res$overall_results$se > 0)
})

# Panel QTT with full bootstrap: SEs are positive, CIs bracket the estimates.
test_that("cic panel QTT SE path produces valid intervals", {
  skip_slow()
  dat <- make_staggered_panel()
  set.seed(42)
  res <- cic(yname = "y", gname = "group", tname = "period", idname = "id",
             data = dat, biters = 100, gt_type = "qtt", probs = c(0.25, 0.5, 0.75))
  expect_s3_class(res, "pte_qtt")
  expect_true(all(res$overall$se > 0))
  expect_true(all(res$overall$lower < res$overall$qtt))
  expect_true(all(res$overall$upper > res$overall$qtt))
  expect_true(all(res$dynamic$se > 0))
  expect_true(all(res$group$se > 0))
})
