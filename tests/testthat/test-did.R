# =============================================================================
# Title: Fast smoke tests — DiD estimators
# Description: Regression checks for CiC, QDiD, MDiD, ddid2, panel.qtet.
#   Uses se=FALSE to skip bootstrap. Reference values from first verified run.
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

# CiC and QDiD: canonical two-period DiD estimators; ATT values are stable
# reference points for the PSID observational data.
test_that("CiC and QDiD return correct structure and stable ATT", {
  cic1 <- CiC(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    se = FALSE, probs = seq(0.25, 0.75, 0.25)
  )
  expect_s3_class(cic1, "QTE")
  expect_length(cic1$qte, 3)
  expect_true(is.numeric(cic1$ate))
  expect_false(anyNA(cic1$qte))
  expect_equal(cic1$ate, 5089.643, tolerance = 1e-3)

  qdid1 <- QDiD(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    se = FALSE, probs = seq(0.25, 0.75, 0.25)
  )
  expect_s3_class(qdid1, "QTE")
  expect_length(qdid1$qte, 3)
  expect_equal(qdid1$ate, 1684.959, tolerance = 1e-3)
})

# MDiD and ddid2: alternative DiD estimators; both should recover the mean DiD
# ATT (2326.506) on this dataset, which serves as a cross-method sanity check.
test_that("MDiD and ddid2 return correct structure and stable ATT", {
  md1 <- MDiD(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    se = FALSE, probs = seq(0.25, 0.75, 0.25)
  )
  expect_s3_class(md1, "QTE")
  expect_true(is.numeric(md1$ate))
  expect_equal(md1$ate, 2326.506, tolerance = 1e-3)

  suppressWarnings({
    dd1 <- ddid2(re ~ treat,
      t = 1978, tmin1 = 1975, tname = "year",
      data = lalonde.psid.panel, idname = "id",
      se = FALSE, probs = seq(0.25, 0.75, 0.25)
    )
  })
  expect_s3_class(dd1, "QTE")
  expect_length(dd1$qte, 3)
  expect_equal(dd1$ate, 2326.506, tolerance = 1e-3)
})

# panel.qtet uses three time periods; tmin2=1974 is the earliest available
# period in lalonde.psid.panel (1974, 1975, 1978).
test_that("panel.qtet returns correct structure and stable ATT", {
  suppressWarnings({
    pq1 <- panel.qtet(re ~ treat,
      t = 1978, tmin1 = 1975, tmin2 = 1974, tname = "year",
      data = lalonde.psid.panel, idname = "id",
      se = FALSE, probs = seq(0.25, 0.75, 0.25)
    )
  })
  expect_s3_class(pq1, "QTE")
  expect_length(pq1$qte, 3)
  expect_true(is.numeric(pq1$ate))
  expect_false(anyNA(pq1$qte))
  expect_equal(pq1$ate, 2326.506, tolerance = 1e-3)
})
