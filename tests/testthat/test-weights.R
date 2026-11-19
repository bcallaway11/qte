# =============================================================================
# Title: Tests for sampling weights in the panel estimators
# Description: Sampling weights in the panel estimators (#15). Integer
#   weights must give the same estimates as duplicating units (row
#   replication); cohort-dependent weights make this also check that cohorts
#   are combined by weighted shares (#13, ptetools #23). Rescaling weights
#   changes nothing; invalid weights error (ptetools #24); and covariate
#   cic_gt()/qdid_gt() must not depend on the row order of gt_data (#14).
# Author: Brant Callaway
# Last update: 2026-11-19
# Date created: 2026-11-19
# =============================================================================

# Two-period panel cell with covariate x and unequal weights, as the
# ptetools subset functions hand it to the *_gt() functions.
weighted_gt_cell <- function(n = 400, seed = 11) {
  set.seed(seed)
  x <- rnorm(n)
  D <- rbinom(n, 1, plogis(0.5 * x))
  u <- rnorm(n, 0.4 * D + 0.6 * x)
  w <- sample(1:4, n, replace = TRUE)
  rbind(
    data.frame(id = 1:n, D = D, name = "pre",  Y = u + x, x = x, .w = w),
    data.frame(id = 1:n, D = D, name = "post", Y = 1 + 1.5 * u + x + D,
               x = x, .w = w)
  )
}

test_that("covariate cic_gt()/qdid_gt() ATT ignores gt_data row order (#14)", {
  gt <- weighted_gt_cell()
  # same cell, rows shuffled within each period (ids unchanged)
  set.seed(2)
  shuffled <- gt[c(sample(which(gt$name == "pre")),
                   sample(which(gt$name == "post"))), ]
  rownames(shuffled) <- NULL

  op <- options(qte.disable_qr_cache = TRUE)
  on.exit(options(op))
  for (f in list(cic_gt, qdid_gt)) {
    sorted_att   <- suppressWarnings(f(gt, ~x))$attgt
    shuffled_att <- suppressWarnings(f(shuffled, ~x))$attgt
    expect_equal(shuffled_att, sorted_att, tolerance = 1e-8)
  }
})


# --- Row replication -----------------------------------------------------

# Balanced panel: cohorts 3 and 4 and never-treated units, 4 periods, an
# integer-valued outcome, and integer weights that depend on the cohort, so
# weighted and unweighted cohort shares differ. The integer outcome puts
# every support point on the grid ptetools evaluates QTTs on (quantiles of
# the pooled outcome), so QTTs can be compared exactly too.
weighted_panel <- function() {
  set.seed(4)
  n <- 300
  g <- sample(c(0, 3, 4), n, replace = TRUE, prob = c(0.4, 0.3, 0.3))
  w <- ifelse(g == 3, sample(2:3, n, TRUE), sample(1:2, n, TRUE))
  d <- expand.grid(id = seq_len(n), period = 1:4)
  d$g <- g[d$id]
  d$w <- w[d$id]
  d$x <- rnorm(n)[d$id]
  d$y <- round(2 * (d$x + rnorm(n)[d$id] + 0.3 * d$period +
                      (d$g > 0 & d$period >= d$g) * (1 + d$g / 3) +
                      rnorm(nrow(d))))
  d
}

# each unit duplicated w times, every copy with its own id
replicate_units <- function(d) {
  first <- d[d$period == 1, ]
  ids <- rep(first$id, first$w)
  out <- d[unlist(lapply(ids, function(i) which(d$id == i))), ]
  out$id <- rep(seq_along(ids), each = 4)
  out$w <- NULL
  out
}

fit_panel <- function(estimator, data, ...) {
  set.seed(1)
  suppressMessages(suppressWarnings(estimator(
    yname = "y", gname = "g", tname = "period", idname = "id", data = data,
    biters = 2, ...
  )))
}

# cond_dist = "dr" with covariates: weighted distribution regression equals
# its fit on replicated rows. (Weighted quantile regression does not match
# exactly: the linear program's solution is not unique.)
estimators <- list(
  cic = list(f = cic, dr = TRUE), qdid = list(f = qdid, dr = TRUE),
  mdid = list(f = mdid, dr = FALSE), ddid = list(f = ddid, dr = TRUE),
  panel_qtt = list(f = panel_qtt, dr = TRUE),
  lou_qtt = list(f = lou_qtt, dr = TRUE),
  # cdf_did has no cond_dist; with covariates its default (aipw) uses
  # distribution regression, which matches replication exactly
  cdf_did = list(f = cdf_did, dr = FALSE)
)

test_that("integer weights equal row replication: ATT and event study", {
  d <- weighted_panel()
  dr <- replicate_units(d)
  op <- options(qte.disable_qr_cache = TRUE)
  on.exit(options(op))
  for (nm in names(estimators)) {
    e <- estimators[[nm]]
    for (xf in list(~1, ~x)) {
      extra <- if (e$dr && length(all.vars(xf)) > 0) list(cond_dist = "dr")
      a <- do.call(fit_panel, c(list(e$f, d, weightsname = "w",
                                     xformula = xf), extra))
      b <- do.call(fit_panel, c(list(e$f, dr, xformula = xf), extra))
      label <- paste(nm, deparse(xf))
      expect_equal(a$overall_results$att, b$overall_results$att,
                   tolerance = 1e-6, label = label)
      expect_equal(a$dyn_results$att.e, b$dyn_results$att.e,
                   tolerance = 1e-6, label = label)
    }
  }
})

test_that("integer weights equal row replication: QTT", {
  d <- weighted_panel()
  dr <- replicate_units(d)
  op <- options(qte.disable_qr_cache = TRUE)
  on.exit(options(op))
  # mdid() is left out: its counterfactual outcomes (pre-period outcome plus a
  # mean shift) are off the outcome grid, so its QTT moves with the grid,
  # which differs between the two data sets
  for (nm in setdiff(names(estimators), "mdid")) {
    e <- estimators[[nm]]
    for (xf in list(~1, ~x)) {
      extra <- if (e$dr && length(all.vars(xf)) > 0) list(cond_dist = "dr")
      args <- c(list(xformula = xf, gt_type = "qtt",
                     probs = c(0.25, 0.5, 0.75)), extra)
      a <- do.call(fit_panel, c(list(e$f, d, weightsname = "w"), args))
      b <- do.call(fit_panel, c(list(e$f, dr), args))
      expect_equal(a$overall$qtt, b$overall$qtt, tolerance = 1e-6,
                   label = paste(nm, deparse(xf)))
    }
  }
})

test_that("panel_qtt long copula weights cohorts by weighted shares (#13)", {
  d <- weighted_panel()
  dr <- replicate_units(d)
  args <- list(gt_type = "qtt", pre_copula = "long",
               probs = c(0.25, 0.5, 0.75))
  a <- do.call(fit_panel, c(list(panel_qtt, d, weightsname = "w"), args))
  b <- do.call(fit_panel, c(list(panel_qtt, dr), args))
  expect_equal(a$overall$qtt, b$overall$qtt, tolerance = 1e-6)
  expect_equal(a$dynamic$qtt, b$dynamic$qtt, tolerance = 1e-6)
})

# --- Scale invariance and validation -------------------------------------

test_that("rescaling the weights changes neither estimates nor SEs", {
  d <- weighted_panel()
  d7 <- d
  d7$w <- 7 * d7$w
  for (e in list(cic, ddid)) {
    for (type in c("att", "qtt")) {
      args <- list(gt_type = type, probs = c(0.25, 0.5, 0.75))
      a <- do.call(fit_panel, c(list(e, d, weightsname = "w"), args))
      b <- do.call(fit_panel, c(list(e, d7, weightsname = "w"), args))
      if (type == "att") {
        expect_equal(a$overall_results, b$overall_results)
      } else {
        expect_equal(a$overall, b$overall)
      }
    }
  }
})

test_that("invalid sampling weights error (ptetools #24)", {
  d <- weighted_panel()
  varying <- d
  varying$w <- varying$w + (varying$period == 2)
  expect_error(fit_panel(cic, varying, weightsname = "w"), "constant")
  negative <- d
  negative$w[1] <- -1
  expect_error(fit_panel(cic, negative, weightsname = "w"), "non-negative")
})
