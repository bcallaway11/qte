# =============================================================================
# Title: Fast smoke tests — unconfoundedness QTE estimator
# Description: Regression checks for unc_qte (Firpo 2007 IPW estimator).
#   Uses a low biters count to keep the test fast. Reference ATE from first
#   verified run. The large negative value reflects the PSID observational
#   comparison group — a regression check that estimates stay stable, not a
#   causal claim. Also checks that the AIPW ATT is doubly robust (issue #11).
# Author: Brant Callaway
# Last update: 2026-11-15
# Date created: 2026-05-18
# =============================================================================

# Smoke test: no covariates (simple quantile differences with bootstrap SE)
test_that("unc_qte returns correct structure and stable ATE", {
  set.seed(42)
  cq1 <- unc_qte(
    yname = "re78", dname = "treat",
    data  = lalonde.psid,
    probs = seq(0.25, 0.75, 0.25),
    biters = 20
  )
  expect_s3_class(cq1, "QTE")
  expect_length(cq1$qte, 3)
  expect_true(is.numeric(cq1$ate))
  expect_false(anyNA(cq1$qte))
  expect_equal(cq1$ate, -15204.78, tolerance = 1e-2)
})

# Smoke test: with covariates (propensity score reweighting); checks pscore fitted
test_that("unc_qte with xformla returns valid IPW estimates", {
  set.seed(42)
  cq2 <- unc_qte(
    yname   = "re78", dname = "treat",
    data    = lalonde.psid,
    xformla = ~ age + education + black + hispanic + married + nodegree,
    probs   = seq(0.25, 0.75, 0.25),
    biters  = 20
  )
  expect_s3_class(cq2, "QTE")
  expect_length(cq2$qte, 3)
  expect_false(anyNA(cq2$qte))
  expect_true(!is.null(cq2$pscore.reg))
  expect_equal(cq2$ate, -13115.59, tolerance = 1e-2)
})

test_that("the AIPW ATT is doubly robust: right with a wrong propensity model (#11)", {
  # selection through sin(2z) (the logit on z, z^2 is wrong); outcome linear
  # in (z, z^2) (right). IPW is badly biased here; AIPW should not be.
  set.seed(12)
  n <- 3000
  z <- rnorm(n)
  D <- rbinom(n, 1, plogis(-0.3 + 1.5 * sin(2 * z) + z))
  d <- data.frame(z = z, z2 = z^2, D = D, Y = 1 + z + 0.8 * z^2 + rnorm(n) + D)
  att <- function(em) {
    suppressWarnings(unc_qte(yname = "Y", dname = "D", data = d,
                             xformla = ~ z + z2, est_method = em,
                             target = "qtt", biters = 2, probs = 0.5))$ate
  }
  expect_gt(abs(att("ipw") - 1), 0.5)
  expect_equal(att("aipw"), 1, tolerance = 0.15)
})
