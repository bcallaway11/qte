# =============================================================================
# Title: Fast smoke tests — unconfoundedness QTE estimator
# Description: Regression checks for unc_qte (Firpo 2007 IPW estimator).
#   Uses a low biters count to keep the test fast. Reference ATE from first
#   verified run. The large negative value reflects the PSID observational
#   comparison group — a regression check that estimates stay stable, not a
#   causal claim. Also checks that the AIPW ATT is doubly robust (issue #11)
#   that a non-0/1 treatment errors clearly (issue #1), and that rows with
#   missing values are dropped.
# Author: Brant Callaway
# Last update: 2026-11-19
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

test_that("sampling weights act like replicated rows in every unc_qte branch", {
  # integer weights must give the same estimates as the unweighted fit on
  # row-replicated data (distribution regression: unique solutions, so exact)
  set.seed(5)
  n <- 600
  z <- rnorm(n)
  D <- rbinom(n, 1, plogis(0.3 + 0.8 * z))
  wi <- sample(1:3, n, replace = TRUE)
  d <- data.frame(z = z, D = D, Y = 1 + z + 0.5 * z^2 + rnorm(n) + D,
                  wt = wi / mean(wi))
  big <- d[rep(seq_len(n), wi), ]
  for (em in c("ipw", "or", "aipw")) {
    for (tg in c("qte", "qtt")) {
      est <- function(dat, wn) {
        r <- suppressWarnings(unc_qte(yname = "Y", dname = "D", data = dat,
                                      xformla = ~z, weightsname = wn,
                                      est_method = em, target = tg,
                                      cond_dist = "dr", biters = 2,
                                      probs = c(0.25, 0.5, 0.75)))
        c(r$qte, r$ate)
      }
      expect_equal(est(d, "wt"), est(big, NULL), tolerance = 1e-8)
    }
  }
})

test_that("unc_qte errors informatively when dname is not coded 0/1 (#1)", {
  set.seed(42)
  df <- data.frame(y = rnorm(50), d = rep(1:2, each = 25))
  expect_error(unc_qte(yname = "y", dname = "d", data = df, biters = 2),
               "coded 0 \\(untreated\\) / 1 \\(treated\\); found values: 1, 2")
  expect_error(unc_qte(yname = "y", dname = "nope", data = df, biters = 2),
               "is not a column of `data`")
  df$d <- df$d == 2
  expect_no_error(suppressWarnings(
    unc_qte(yname = "y", dname = "d", data = df, biters = 2)
  ))
})

test_that("unc_qte drops rows with missing values, with a warning", {
  d <- lalonde.psid
  d$re78[1:3] <- NA
  set.seed(1)
  expect_warning(
    r_na <- unc_qte(yname = "re78", dname = "treat", data = d,
                    probs = 0.5, biters = 2),
    "dropping 3 row"
  )
  set.seed(1)
  r_cc <- unc_qte(yname = "re78", dname = "treat", data = d[-(1:3), ],
                  probs = 0.5, biters = 2)
  expect_equal(r_na$qte, r_cc$qte)
  expect_equal(r_na$qte.se, r_cc$qte.se)
})
