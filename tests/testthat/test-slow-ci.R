# =============================================================================
# Title: Slow tests — CI estimators with bootstrap SE
# Description: Verifies SE paths for ci.qte and ci.qtet.
#   Uses iters=5 to keep runtime manageable. Gate with R_SLOW_TESTS=true.
# Author: Brant Callaway
# Last update: 2026-05-22
# Date created: 2026-05-18
# =============================================================================

test_that("ci.qte SE path produces valid intervals", {
  skip_slow()
  set.seed(42)
  cq1 <- ci.qte(re78 ~ treat,
    data = lalonde.psid, se = TRUE, iters = 5,
    probs = seq(0.25, 0.75, 0.25)
  )
  expect_s3_class(cq1, "QTE")
  expect_true(all(cq1$qte.se > 0))
  expect_true(all(cq1$qte.lower < cq1$qte))
  expect_true(all(cq1$qte.upper > cq1$qte))
})

test_that("ci.qtet SE path produces valid intervals", {
  skip_slow()
  set.seed(42)
  cqt1 <- ci.qtet(re78 ~ treat,
    data = lalonde.psid, se = TRUE, iters = 5,
    probs = seq(0.25, 0.75, 0.25)
  )
  expect_s3_class(cqt1, "QTE")
  expect_true(all(cqt1$qte.se > 0))
  expect_true(all(cqt1$qte.lower < cqt1$qte))
  expect_true(all(cqt1$qte.upper > cqt1$qte))
})
