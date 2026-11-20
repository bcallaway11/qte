# =============================================================================
# Title: Fast smoke tests — unconfoundedness QTE estimator
# Description: Regression checks for unc_qte (Firpo 2007 IPW estimator).
#   Uses a low biters count to keep the test fast. Reference ATE from first
#   verified run. The large negative value reflects the PSID observational
#   comparison group — a regression check that estimates stay stable, not a
#   causal claim. Also checks that the AIPW ATT is doubly robust (issue #11)
#   that a non-0/1 treatment errors clearly (issue #1), the uniform vs
#   pointwise bands (cband), that rows with missing values are dropped,
#   clustervars, the deprecated xformla alias, the slim pscore.reg, and the
#   bootstrap scale at quantiles inside a mass point.
# Author: Brant Callaway
# Last update: 2026-11-19
# Date created: 2026-05-18
# =============================================================================

# Smoke test: no covariates (simple quantile differences with bootstrap SE)
test_that("unc_qte returns correct structure and stable ATE", {
  set.seed(42)
  cq1 <- quiet_discrete(unc_qte(
    yname = "re78", dname = "treat",
    data  = lalonde.psid,
    probs = seq(0.25, 0.75, 0.25),
    biters = 20
  ))
  expect_s3_class(cq1, "QTE")
  expect_length(cq1$qte, 3)
  expect_true(is.numeric(cq1$ate))
  expect_false(anyNA(cq1$qte))
  expect_equal(cq1$ate, -15204.78, tolerance = 1e-2)
})

# Smoke test: with covariates (propensity score reweighting); checks pscore fitted
test_that("unc_qte with xformula returns valid IPW estimates", {
  set.seed(42)
  cq2 <- quiet_discrete(unc_qte(
    yname   = "re78", dname = "treat",
    data    = lalonde.psid,
    xformula = ~ age + education + black + hispanic + married + nodegree,
    probs   = seq(0.25, 0.75, 0.25),
    biters  = 20
  ))
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
                             xformula = ~ z + z2, est_method = em,
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
  for (em in c("ipw", "reg", "aipw")) {
    for (tg in c("qte", "qtt")) {
      est <- function(dat, wn) {
        r <- suppressWarnings(unc_qte(yname = "Y", dname = "D", data = dat,
                                      xformula = ~z, weightsname = wn,
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

test_that("unc_qte cband = TRUE gives a uniform band qte +/- c * se", {
  probs <- c(0.25, 0.5, 0.75)
  fit <- function(cband) {
    set.seed(7)
    suppressWarnings(unc_qte(yname = "re78", dname = "treat",
                             data = lalonde.psid, probs = probs,
                             biters = 20, cband = cband))
  }
  ub <- fit(TRUE)
  pw <- fit(FALSE)
  z <- qnorm(0.975)

  # same bootstrap draws, so the same estimates, SEs and critical value
  expect_equal(ub$qte.se, pw$qte.se)
  expect_equal(ub$c, pw$c)
  expect_true(isTRUE(ub$cband))
  expect_false(isTRUE(pw$cband))

  # uniform band uses c, pointwise uses z
  expect_equal(ub$qte.upper, ub$qte + ub$c * ub$qte.se)
  expect_equal(ub$qte.lower, ub$qte - ub$c * ub$qte.se)
  expect_equal(pw$qte.upper, pw$qte + z * pw$qte.se)

  # summary labels the band correctly
  expect_output(print(summary(ub)), "Simult.")
  expect_output(print(summary(pw)), "Pointwise")
})

test_that("unc_qte drops rows with missing values, with a warning", {
  d <- lalonde.psid
  d$re78[1:3] <- NA
  set.seed(1)
  expect_warning(
    r_na <- quiet_discrete(unc_qte(yname = "re78", dname = "treat",
                                   data = d, probs = 0.5, biters = 2)),
    "dropping 3 row"
  )
  set.seed(1)
  r_cc <- quiet_discrete(unc_qte(yname = "re78", dname = "treat",
                                 data = d[-(1:3), ], probs = 0.5,
                                 biters = 2))
  expect_equal(r_na$qte, r_cc$qte)
  expect_equal(r_na$qte.se, r_cc$qte.se)
})

test_that("unc_qte clustervars: row-level clusters match, cluster shocks raise SEs", {
  probs <- c(0.25, 0.5, 0.75)
  fit <- function(data, ...) {
    set.seed(5)
    suppressWarnings(unc_qte(yname = "y", dname = "d", data = data,
                             probs = probs, biters = 40, ...))
  }

  # 40 clusters with a common outcome shock and treatment assigned by
  # cluster (with treatment varying within clusters, the shocks would mostly
  # cancel between treated and untreated and clustering would matter little)
  set.seed(1)
  n_cl <- 40
  m <- 25
  cl_id <- rep(seq_len(n_cl), each = m)
  d <- data.frame(cl = cl_id, row = seq_len(n_cl * m),
                  d = rep(rep(0:1, n_cl / 2), each = m))
  d$y <- 2 * rnorm(n_cl)[cl_id] + rnorm(nrow(d)) + d$d

  # one cluster per row, in row order: exactly the unclustered bootstrap
  unclustered <- fit(d)
  by_row <- fit(d, clustervars = "row")
  expect_equal(by_row$qte.se, unclustered$qte.se)
  expect_equal(by_row$ate.se, unclustered$ate.se)

  # resampling whole clusters picks up the common shocks
  clustered <- fit(d, clustervars = "cl")
  expect_equal(clustered$qte, unclustered$qte)
  expect_gt(mean(clustered$qte.se), 1.5 * mean(unclustered$qte.se))

  expect_error(fit(d, clustervars = c("cl", "row")), "one column")
  expect_error(fit(d, clustervars = "nope"), "not found")
})

test_that("unc_qte's deprecated xformla still works, with a warning", {
  args <- list(yname = "re78", dname = "treat", data = lalonde.psid,
               probs = c(0.25, 0.5, 0.75), biters = 4)
  set.seed(3)
  new <- suppressWarnings(do.call(unc_qte, c(args, xformula = ~ age)))
  set.seed(3)
  expect_warning(old <- do.call(unc_qte, c(args, xformla = ~ age)),
                 "deprecated")
  expect_equal(old$qte, new$qte)
  expect_equal(old$qte.se, new$qte.se)
  expect_error(do.call(unc_qte, c(args, xformula = ~ age, xformla = ~ age)),
               "only")
})

test_that("unc_qte's pscore.reg is slim but still a usable glm", {
  set.seed(2)
  xf <- ~ age + education + black + nodegree
  res <- suppressWarnings(unc_qte(yname = "re78", dname = "treat",
                                  data = lalonde.psid, xformula = xf,
                                  probs = c(0.25, 0.5, 0.75), biters = 2))
  ps <- res$pscore.reg
  full <- glm(treat ~ age + education + black + nodegree, data = lalonde.psid,
              family = quasibinomial())

  expect_s3_class(ps, "glm")
  expect_equal(coef(ps), coef(full))
  expect_equal(summary(ps)$coefficients, summary(full)$coefficients)
  expect_equal(unname(fitted(ps)), unname(fitted(full)))
  expect_equal(predict(ps, newdata = lalonde.psid[1:5, ], type = "response"),
               predict(full, newdata = lalonde.psid[1:5, ], type = "response"))
  # no copies of the data, and the formula no longer carries the
  # estimator's working frame
  expect_null(ps$data)
  expect_null(ps$model)
  expect_lt(length(serialize(res, NULL)), 1e6)
})

test_that("unc_qte's deprecated est_method = \"or\" still works, as \"reg\"", {
  args <- list(yname = "re78", dname = "treat", data = lalonde.psid,
               xformula = ~ age + education, probs = c(0.25, 0.5, 0.75),
               biters = 2)
  set.seed(3)
  new <- quiet_discrete(do.call(unc_qte, c(args, est_method = "reg")))
  set.seed(3)
  expect_warning(
    old <- quiet_discrete(do.call(unc_qte, c(args, est_method = "or"))),
    "deprecated"
  )
  expect_equal(old$qte, new$qte)
})

test_that("bootstrap scale: IQR-based, SD inside a mass point, 0 if constant", {
  set.seed(1)
  B <- 100
  normal <- rnorm(B, sd = 2)
  # inside a mass point: nearly constant, a few draws jump
  mass <- c(rep(6e-5, 92), rnorm(B - 92, -170, 10)) + rnorm(B, sd = 1e-6)
  constant <- rep(0, B)
  draws <- cbind(normal, mass, constant)
  sc <- unname(boot_curve_scale(draws))

  iqr <- unname(diff(quantile(normal, c(0.25, 0.75), type = 1))) /
    (qnorm(0.75) - qnorm(0.25))
  expect_equal(sc[1], iqr)
  expect_equal(sc[2], sd(mass))
  expect_equal(sc[3], 0)

  # the sup-t critical value stays moderate: with the IQR scale at the
  # mass-point quantile it was about 1e8
  est <- c(0, 6e-5)
  cb <- apply(draws[, 1:2], 1, function(q) max(abs((q - est) / sc[1:2])))
  expect_lt(quantile(cb, 0.95, type = 1), sqrt(B))
})
