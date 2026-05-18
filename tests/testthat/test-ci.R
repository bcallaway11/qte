# =============================================================================
# Title: Fast smoke tests — CI, spatt, bounds estimators
# Description: Regression checks for ci.qte, ci.qtet, spatt, bounds.
#   Uses se=FALSE to skip bootstrap. Reference values from first verified run.
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

# ci.qte and ci.qtet: cross-sectional randomization-based estimators.
# The large negative values reflect the PSID observational comparison group —
# a regression check that these stay stable, not a causal claim.
test_that("ci.qte and ci.qtet return correct structure and stable ATE", {
  cq1 <- ci.qte(re78 ~ treat,
    data = lalonde.psid, se = FALSE,
    probs = seq(0.25, 0.75, 0.25)
  )
  expect_s3_class(cq1, "QTE")
  expect_length(cq1$qte, 3)
  expect_true(is.numeric(cq1$ate))
  expect_false(anyNA(cq1$qte))
  expect_equal(cq1$ate, -15204.78, tolerance = 1e-2)

  cqt1 <- ci.qtet(re78 ~ treat,
    data = lalonde.psid, se = FALSE,
    probs = seq(0.25, 0.75, 0.25)
  )
  expect_s3_class(cqt1, "QTE")
  expect_length(cqt1$qte, 3)
  expect_equal(cqt1$ate, -15204.78, tolerance = 1e-2)
})

# spatt: Abadie (2005) semiparametric DiD; should match mean DiD ATT without
# covariates. bounds returns BoundsObj, not QTE.
test_that("spatt and bounds return correct structure and stable ATT", {
  suppressWarnings({
    sp1 <- spatt(re ~ treat,
      t = 1978, tmin1 = 1975, tname = "year",
      data = lalonde.psid.panel, idname = "id", se = FALSE
    )
  })
  expect_s3_class(sp1, "QTE")
  expect_true(is.numeric(sp1$ate))
  expect_equal(sp1$ate, 2326.506, tolerance = 1e-3)

  b1 <- bounds(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    probs = seq(0.25, 0.75, 0.25)
  )
  expect_s3_class(b1, "BoundsObj")
  expect_true(is.numeric(b1$att))
  expect_length(b1$lb.qte, 3)
  expect_length(b1$ub.qte, 3)
  # lower bound < upper bound at each quantile
  expect_true(all(b1$lb.qte < b1$ub.qte))
})
