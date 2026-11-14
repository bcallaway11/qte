# =============================================================================
# Title: Tests for conditional quantile helpers (issue #9)
# Description: qr_cdf()/qr_quantile() must give monotone, correctly indexed
#   conditional CDFs and quantiles (including for a single row of newdata),
#   and covariate-adjusted cic_gt() must run with one treated unit.
# Author: Brant Callaway
# Last update: 2026-11-14
# Date created: 2026-11-14
# =============================================================================

test_that("qr helpers: correct indexing, monotone under crossing, one row ok", {
  set.seed(1)
  d <- data.frame(x = rnorm(200))
  d$y <- 1 + d$x + rnorm(200) * exp(0.5 * d$x)
  taus <- seq(0.05, 0.95, 0.05)
  fit <- suppressWarnings(quantreg::rq(y ~ x, data = d, tau = taus))
  nd <- data.frame(x = c(-3, 0, 3))
  q <- qr_sorted_quantiles(fit, nd)

  # F evaluated at the k-th predicted quantile is tau_k (quantreg's
  # multi-row Fhat returns tau_{k-1})
  expect_equal(qr_cdf(q[c(2, 2), ], taus, q[2, c(5, 10)]), taus[c(5, 10)])
  expect_equal(qr_quantile(q[c(2, 2), ], taus, taus[c(5, 10)]), q[2, c(5, 10)])

  # monotone in y and in p for every row, crossing or not
  y_grid <- seq(-10, 10, length.out = 50)
  for (i in seq_len(nrow(q))) {
    expect_false(is.unsorted(qr_cdf(q[rep(i, 50), ], taus, y_grid)))
    expect_false(is.unsorted(qr_quantile(q[rep(i, 19), ], taus, taus)))
  }

  # a single row of newdata works
  expect_equal(dim(qr_sorted_quantiles(fit, nd[1, , drop = FALSE])),
               c(1L, length(taus)))
})

test_that("covariate cic_gt runs with a single treated unit (#9)", {
  set.seed(1)
  n <- 301
  X <- matrix(rnorm(n * 6), n, 6, dimnames = list(NULL, paste0("x", 1:6)))
  X[n, ] <- 3 # treated unit far out in covariate space -> quantiles cross
  D <- rep(c(0L, 1L), c(n - 1, 1))
  y_pre <- drop(X %*% rep(0.3, 6)) + rnorm(n) * exp(0.4 * X[, 1])
  gt <- data.frame(id = rep(seq_len(n), 2),
                   name = rep(c("pre", "post"), each = n), D = rep(D, 2),
                   Y = c(y_pre, y_pre + 0.5 + rnorm(n) + D), .w = 1,
                   rbind(X, X))
  res <- suppressWarnings(cic_gt(gt, xformula = ~ x1 + x2 + x3 + x4 + x5 + x6))
  expect_true(is.finite(res$attgt))
})
