# =============================================================================
# Title: Tests for sampling weights in the panel estimators
# Description: Checks that sampling weights stay attached to the right units.
#   cic_gt()/qdid_gt() with covariates must not depend on the row order of
#   gt_data (#14). More weight tests are planned in #15.
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
