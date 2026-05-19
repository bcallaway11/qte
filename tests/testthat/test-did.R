# =============================================================================
# Title: Fast smoke tests — DiD estimators
# Description: Regression checks for CiC, QDiD, MDiD, ddid2, panel.qtet.
#   Uses se=FALSE to skip bootstrap. Reference values from first verified run.
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

# Shared synthetic staggered panel used by qdid tests below.
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

# qdid panel ATT: pte_emp_boot, 6 ATT(g,t) cells, overall ATT near 1.
test_that("qdid panel ATT returns valid pte_emp_boot with stable estimate", {
  dat <- make_staggered_panel()
  set.seed(42)
  res <- qdid(yname = "y", gname = "group", tname = "period", idname = "id",
              data = dat, biters = 20)
  expect_s3_class(res, "pte_emp_boot")
  expect_true(is.numeric(res$overall_results$att))
  expect_false(anyNA(res$attgt_results$att))
  expect_equal(nrow(res$attgt_results), 6L)
  expect_equal(res$overall_results$att, 1.131815, tolerance = 1e-4)
})

# qdid QTT mode: pte_qtt class, correct structure, stable median.
test_that("qdid panel QTT returns valid pte_qtt with correct structure", {
  dat <- make_staggered_panel()
  set.seed(42)
  res <- qdid(yname = "y", gname = "group", tname = "period", idname = "id",
              data = dat, biters = 20, gt_type = "qtt", probs = c(0.25, 0.5, 0.75))
  expect_s3_class(res, "pte_qtt")
  expect_equal(nrow(res$overall), 3L)
  expect_false(anyNA(res$overall$qtt))
  expect_true(all(res$overall$lower < res$overall$upper))
  expect_equal(res$overall$qtt[res$overall$probs == 0.5], 0.8957124, tolerance = 1e-4)
})

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
