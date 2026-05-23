# =============================================================================
# Title: Slow tests — DiD estimators with bootstrap SE
# Description: Verifies SE paths for CiC, QDiD, MDiD, ddid2, panel.qtet.
#   Uses iters=5 to keep runtime manageable. Gate with R_SLOW_TESTS=true.
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

test_that("CiC SE path produces valid intervals", {
  skip_slow()
  set.seed(42)
  cic1 <- CiC(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    se = TRUE, iters = 5, probs = seq(0.25, 0.75, 0.25)
  )
  expect_s3_class(cic1, "QTE")
  expect_true(all(cic1$qte.se > 0))
  expect_true(all(cic1$qte.lower < cic1$qte))
  expect_true(all(cic1$qte.upper > cic1$qte))
  expect_true(cic1$ate.se > 0)
})

test_that("QDiD SE path produces valid intervals", {
  skip_slow()
  set.seed(42)
  qdid1 <- QDiD(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    se = TRUE, iters = 5, probs = seq(0.25, 0.75, 0.25)
  )
  expect_s3_class(qdid1, "QTE")
  expect_true(all(qdid1$qte.se > 0))
  expect_true(all(qdid1$qte.lower < qdid1$qte))
  expect_true(all(qdid1$qte.upper > qdid1$qte))
  expect_true(qdid1$ate.se > 0)
})

test_that("MDiD SE path produces valid intervals", {
  skip_slow()
  set.seed(42)
  md1 <- MDiD(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    se = TRUE, iters = 5, probs = seq(0.25, 0.75, 0.25)
  )
  expect_s3_class(md1, "QTE")
  expect_true(all(md1$qte.se > 0))
  expect_true(md1$ate.se > 0)
})

test_that("ddid2 SE path produces valid intervals", {
  skip_slow()
  set.seed(42)
  dd1 <- ddid2(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    se = TRUE, iters = 5, probs = seq(0.25, 0.75, 0.25)
  )
  expect_s3_class(dd1, "QTE")
  expect_true(all(dd1$qte.se > 0))
  expect_true(all(dd1$qte.lower < dd1$qte))
  expect_true(all(dd1$qte.upper > dd1$qte))
  expect_true(dd1$ate.se > 0)
})

test_that("panel.qtet SE path produces valid intervals", {
  skip_slow()
  set.seed(42)
  suppressWarnings({
    pq1 <- panel.qtet(re ~ treat,
      t = 1978, tmin1 = 1975, tmin2 = 1974, tname = "year",
      data = lalonde.psid.panel, idname = "id",
      se = TRUE, iters = 5, probs = seq(0.25, 0.75, 0.25)
    )
  })
  expect_s3_class(pq1, "QTE")
  expect_true(all(pq1$qte.se > 0))
  expect_true(all(pq1$qte.lower < pq1$qte))
  expect_true(all(pq1$qte.upper > pq1$qte))
  expect_true(pq1$ate.se > 0)
})
