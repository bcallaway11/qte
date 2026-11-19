# =============================================================================
# Title: Tests for clustered standard errors
# Description: clustervars is passed on to ptetools (ptetools #6): clustering
#   on the unit id reproduces unclustered standard errors, and with
#   cluster-by-period shocks clustering raises them.
# Author: Brant Callaway
# Last update: 2026-11-19
# Date created: 2026-11-19
# =============================================================================

cluster_data <- function() {
  set.seed(10)
  n_cl <- 30
  per_cl <- 10
  n <- n_cl * per_cl
  n_periods <- 4
  cl <- rep(seq_len(n_cl), each = per_cl)
  g_cl <- sample(c(0, 3), n_cl, replace = TRUE)
  d <- expand.grid(id = seq_len(n), period = seq_len(n_periods))
  d$cl <- cl[d$id]
  d$grp <- g_cl[d$cl]
  shock <- matrix(rnorm(n_cl * n_periods), n_cl)
  d$y <- rnorm(n)[d$id] + shock[cbind(d$cl, d$period)] +
    (d$grp > 0 & d$period >= d$grp) + rnorm(nrow(d))
  d
}

run_qtt <- function(estimator, ...) {
  set.seed(1)
  suppressMessages(suppressWarnings(estimator(
    yname = "y", gname = "grp", tname = "period", idname = "id",
    data = cluster_data(), gt_type = "qtt", probs = c(0.25, 0.5, 0.75),
    biters = 40, ...
  )))$overall$se
}

test_that("clustervars clusters the QTT standard errors", {
  for (estimator in list(ddid, cic)) {
    unclustered <- run_qtt(estimator)
    expect_equal(run_qtt(estimator, clustervars = "id"), unclustered)
    expect_gt(mean(run_qtt(estimator, clustervars = "cl")),
              1.3 * mean(unclustered))
  }
})
