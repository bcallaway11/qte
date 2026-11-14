# =============================================================================
# Title: Tests for the fast rank/quantile helpers and quantile-regression reuse
# Description: first_geq(), wcdf_at() and wquant() must match the loops they
#   replace, and reusing quantile regressions across (g,t) cells must not
#   change estimates.
# Author: Brant Callaway
# Last update: 2026-11-14
# Date created: 2026-11-14
# =============================================================================

test_that("fast rank and quantile helpers match the loop definitions", {
  set.seed(1)
  y <- round(rnorm(300), 1)                     # ties
  w <- runif(300)
  w <- w / sum(w)
  x <- c(min(y) - 1, y[1:50], max(y) + 1)       # includes out-of-range values
  expect_equal(wcdf_at(y, w, x), vapply(x, function(v) sum(w * (y <= v)), 1.0))

  cdf <- cumsum(w[order(y)])
  p <- c(0, runif(40), cdf[10], 1, 1.5, NA)     # exact hit, above 1, NA
  loop <- vapply(p, function(pp) {
    idx <- which(cdf >= pp)
    if (length(idx) == 0L) length(cdf) else idx[1L]
  }, integer(1))
  expect_identical(first_geq(cdf, p), loop)
})

test_that("reusing quantile regressions across cells leaves cic_gt unchanged", {
  data(mpdta, package = "did")
  ptep <- ptetools::setup_pte(yname = "lemp", gname = "first.treat", tname = "year",
                              idname = "countyreal", data = mpdta)
  run <- function() {
    res <- suppressWarnings(ptetools::compute.pte(ptep, ptetools::two_by_two_subset,
                                                  cic_gt, xformula = ~lpop))
    vapply(res$attgt.list, function(r) r$att, 1.0)
  }
  old <- options(qte.disable_qr_cache = TRUE)
  on.exit(options(old))
  no_cache <- run()
  options(qte.disable_qr_cache = FALSE)
  expect_identical(run(), no_cache)
})
