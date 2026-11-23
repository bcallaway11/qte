# =============================================================================
# Title: Regression test for degenerate quantile regressions
# Description: A bootstrap resample of lalonde.psid (re78 has a mass point at
#   0) on which quantreg's default Barrodale-Roberts simplex (method = "br")
#   cycles forever at tau = 0.06 -- it hung unc_qte(est_method = "aipw") in
#   the intro vignette. rq_coef() uses Frisch-Newton (method = "fn"), which
#   solves it in well under a second. A regression here shows up as a hang
#   rather than a failure (the "br" loop is in Fortran and ignores time
#   limits). The fixture holds the 2481 resampled row indices.
# Author: Brant Callaway
# Last update: 2026-11-22
# Date created: 2026-11-22
# =============================================================================

test_that("rq_coef() solves a degenerate problem on which rq()'s default cycles", {
  idx <- readRDS(test_path("fixtures", "rq-degenerate-rows.rds"))
  d <- lalonde.psid[idx, ]
  xk <- stats::model.matrix(
    ~ age + I(age^2) + education + black + hispanic + married + nodegree, d
  )
  taus <- seq(0.01, 0.99, 0.01)
  cf <- rq_coef(d$re78, xk, NULL, taus)
  expect_equal(dim(cf), c(ncol(xk), length(taus)))
  expect_true(all(is.finite(cf)))
  # inside the mass point at zero, the fitted quantile at the mean covariates
  # is (essentially) zero
  q_low <- drop(colMeans(xk) %*% cf[, taus <= 0.05])
  expect_true(all(abs(q_low) < 1))
})
