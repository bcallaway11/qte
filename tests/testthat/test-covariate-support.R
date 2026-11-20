# =============================================================================
# Title: Covariate support of the conditional outcome models
# Description: Conditional models fit on comparison units are evaluated at
#   treated units' covariates. Columns that are collinear in the fitting
#   sample are dropped; when treated units don't satisfy the same relation
#   (e.g. a dummy for a region with treated units only), the estimates
#   extrapolate, and every estimator warns once per call, listing the
#   cells (qte#16, via ptetools::unsupported_covariates() and
#   ptetools::pte_cell_warning()). Before, cond_dist = "qr" stopped with
#   "Singular design matrix" and cond_dist = "dr" extrapolated silently.
# Author: Brant Callaway
# Last update: 2026-11-20
# Date created: 2026-11-20
# =============================================================================

# Two periods, treated cohort 2. Region "c" (dummy reg_c) has treated units
# only; x2 equals x1 for comparison units and differs for some treated ones.
support_data <- function(seed = 1, n = 1500) {
  set.seed(seed)
  # G: group variable, named as in the estimators' examples
  G <- ifelse(runif(n) < 0.3, 2, 0) # nolint: object_name_linter
  region <- ifelse(G == 2 & runif(n) < 0.3, "c",
                   sample(c("a", "b"), n, TRUE))
  x1 <- rnorm(n)
  x2 <- ifelse(G == 2 & runif(n) < 0.2, x1 + 1, x1)
  u <- rnorm(n) + (region == "b") + 2 * (region == "c")
  d <- data.frame(id = seq_len(n), G = G, region = region,
                  reg_b = as.numeric(region == "b"),
                  reg_c = as.numeric(region == "c"),
                  zero = 0, x1 = x1, x2 = x2, w = runif(n, 0.5, 2))
  rbind(cbind(d, period = 1, Y = u),
        cbind(d, period = 2, Y = u + 0.5 + rnorm(n, sd = 0.1)))
}

# run an estimator and collect its warnings (messages)
fit_support <- function(f, xformula, data = support_data(), ...) {
  msgs <- character(0)
  set.seed(1)
  res <- withCallingHandlers(
    suppressMessages(f(yname = "Y", gname = "G", tname = "period",
                       idname = "id", data = data, xformula = xformula,
                       gt_type = "qtt", probs = c(0.25, 0.5, 0.75),
                       biters = 3, ...)),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(res = res, warnings = grep("Covariates dropped", msgs, value = TRUE))
}

test_that("a covariate with no comparison variation warns once, QR and DR", {
  ests <- list(
    cic_qr = function(...) cic(..., cond_dist = "qr"),
    cic_dr = function(...) cic(..., cond_dist = "dr"),
    qdid_dr = function(...) qdid(..., cond_dist = "dr"),
    ddid_qr = function(...) ddid(..., cond_dist = "qr"),
    mdid = mdid,
    cdf_did = function(...) cdf_did(..., est_method = "reg"),
    lou_dr = function(...) lou_qtt(..., est_method = "reg", cond_dist = "dr")
  )
  share <- with(subset(support_data(), period == 1 & G == 2), mean(reg_c))
  for (nm in names(ests)) {
    out <- fit_support(ests[[nm]], ~ reg_b + reg_c)
    expect_true(is.finite(out$res$att$overall$att), label = nm)
    expect_length(out$warnings, 1)
    expect_match(out$warnings, "group 2, period 2: reg_c", fixed = TRUE,
                 label = nm)
    expect_match(out$warnings,
                 paste0(round(100 * share, 1), "% of treated units"),
                 fixed = TRUE, label = nm)
  }
})

test_that("factor and character covariates are handled the same way", {
  for (xf in list(~ factor(region), ~ region)) {
    out <- fit_support(cic, xf, cond_dist = "qr")
    expect_length(out$warnings, 1)
    # column factor(region)c or regionc
    expect_match(out$warnings, "region\\)?c \\(")
  }
})

test_that("columns collinear only among comparison units are flagged", {
  out <- fit_support(cic, ~ x1 + x2, cond_dist = "dr")
  expect_length(out$warnings, 1)
  expect_match(out$warnings, "x2 (", fixed = TRUE)
})

test_that("harmless dropped columns give no warning and change nothing", {
  # `zero` is 0 for every unit: dropped, but treated units satisfy that too
  base <- fit_support(cic, ~ reg_b, cond_dist = "qr")
  with_zero <- fit_support(cic, ~ reg_b + zero, cond_dist = "qr")
  expect_length(with_zero$warnings, 0)
  expect_equal(with_zero$res$overall$qtt, base$res$overall$qtt)
  base_dr <- fit_support(cic, ~ reg_b, cond_dist = "dr")
  with_zero_dr <- fit_support(cic, ~ reg_b + zero, cond_dist = "dr")
  expect_length(with_zero_dr$warnings, 0)
  expect_equal(with_zero_dr$res$overall$qtt, base_dr$res$overall$qtt)
})

test_that("shares are weighted by the sampling weights", {
  dat <- support_data()
  out <- fit_support(cic, ~ reg_b + reg_c, data = dat, weightsname = "w")
  tr <- subset(dat, period == 1 & G == 2)
  share_w <- sum(tr$w * tr$reg_c) / sum(tr$w)
  expect_match(out$warnings, paste0(round(100 * share_w, 1), "%"),
               fixed = TRUE)
})

test_that("unc_qte() warns once, not again in bootstrap draws", {
  set.seed(2)
  n <- 1000
  d <- data.frame(treat = rbinom(n, 1, 0.4), x = rnorm(n))
  d$r <- ifelse(d$treat == 1 & runif(n) < 0.25, 1, 0)
  d$y <- d$x + 2 * d$r + d$treat + rnorm(n)
  for (cd in c("qr", "dr")) {
    msgs <- character(0)
    withCallingHandlers(
      unc_qte(yname = "y", dname = "treat", data = d, xformula = ~ x + r,
              est_method = "reg", target = "qtt", cond_dist = cd,
              probs = c(0.25, 0.5, 0.75), biters = 10),
      warning = function(w) {
        msgs <<- c(msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    msgs <- grep("Covariates dropped", msgs, value = TRUE)
    expect_length(msgs, 1)
    expect_match(msgs, "r (", fixed = TRUE)
  }
})
