# =============================================================================
# Title: Compatibility smoke tests — legacy DiD wrappers
# Description: Minimal checks that CiC, QDiD, MDiD, ddid2, panel.qtet emit
#   deprecation warnings and return a QTE object. Numerical regression tests
#   for the modern replacements live in test-cic.R, test-qdid.R, test-mdid.R,
#   test-ddid.R, and test-panel-qtt.R.
# Author: Brant Callaway
# Last update: 2026-05-22
# Date created: 2026-05-18
# =============================================================================

# CiC and QDiD: deprecated wrappers should warn "deprecated" and return QTE.
test_that("CiC and QDiD emit deprecation warnings and return QTE class", {
  expect_warning(
    cic1 <- CiC(re ~ treat,
      t = 1978, tmin1 = 1975, tname = "year",
      data = lalonde.psid.panel, idname = "id",
      se = FALSE, probs = seq(0.25, 0.75, 0.25)
    ),
    regexp = "deprecated", ignore.case = TRUE
  )
  expect_s3_class(cic1, "QTE")

  expect_warning(
    qdid1 <- QDiD(re ~ treat,
      t = 1978, tmin1 = 1975, tname = "year",
      data = lalonde.psid.panel, idname = "id",
      se = FALSE, probs = seq(0.25, 0.75, 0.25)
    ),
    regexp = "deprecated", ignore.case = TRUE
  )
  expect_s3_class(qdid1, "QTE")
})

# MDiD and ddid2: same pattern.
test_that("MDiD and ddid2 emit deprecation warnings and return QTE class", {
  expect_warning(
    md1 <- MDiD(re ~ treat,
      t = 1978, tmin1 = 1975, tname = "year",
      data = lalonde.psid.panel, idname = "id",
      se = FALSE, probs = seq(0.25, 0.75, 0.25)
    ),
    regexp = "deprecated", ignore.case = TRUE
  )
  expect_s3_class(md1, "QTE")

  expect_warning(
    dd1 <- ddid2(re ~ treat,
      t = 1978, tmin1 = 1975, tname = "year",
      data = lalonde.psid.panel, idname = "id",
      se = FALSE, probs = seq(0.25, 0.75, 0.25)
    ),
    regexp = "deprecated", ignore.case = TRUE
  )
  expect_s3_class(dd1, "QTE")
})

# panel.qtet: three-period legacy wrapper.
test_that("panel.qtet emits deprecation warning and returns QTE class", {
  expect_warning(
    pq1 <- panel.qtet(re ~ treat,
      t = 1978, tmin1 = 1975, tmin2 = 1974, tname = "year",
      data = lalonde.psid.panel, idname = "id",
      se = FALSE, probs = seq(0.25, 0.75, 0.25)
    ),
    regexp = "deprecated", ignore.case = TRUE
  )
  expect_s3_class(pq1, "QTE")
})
