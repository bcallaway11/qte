# =============================================================================
# Title: Benchmark — CI and auxiliary estimators
# Description: Times ci.qte, ci.qtet, spatt with and without SE.
#   Run manually: Rscript bench/bench-ci.R
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

devtools::load_all()
data(lalonde.psid, package = "qte")
data(lalonde.psid.panel, package = "qte")

cat("--- ci.qte (se=FALSE) ---\n")
print(system.time(
  ci.qte(re78 ~ treat,
    data = lalonde.psid, se = FALSE,
    probs = seq(0.05, 0.95, 0.05)
  )
))

cat("--- ci.qte (se=TRUE, iters=100) ---\n")
print(system.time(
  ci.qte(re78 ~ treat,
    data = lalonde.psid, se = TRUE, iters = 100,
    probs = seq(0.05, 0.95, 0.05)
  )
))

cat("--- ci.qtet (se=FALSE) ---\n")
print(system.time(
  ci.qtet(re78 ~ treat,
    data = lalonde.psid, se = FALSE,
    probs = seq(0.05, 0.95, 0.05)
  )
))

cat("--- ci.qtet (se=TRUE, iters=100) ---\n")
print(system.time(
  ci.qtet(re78 ~ treat,
    data = lalonde.psid, se = TRUE, iters = 100,
    probs = seq(0.05, 0.95, 0.05)
  )
))

cat("--- spatt (se=FALSE) ---\n")
print(system.time(
  spatt(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id", se = FALSE
  )
))

cat("--- spatt (se=TRUE, analytic) ---\n")
print(system.time(
  spatt(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id", se = TRUE
  )
))
