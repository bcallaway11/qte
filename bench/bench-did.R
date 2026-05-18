# =============================================================================
# Title: Benchmark — DiD estimators
# Description: Times CiC, QDiD, MDiD, ddid2, panel.qtet with and without
#   bootstrap SE. Run manually: Rscript bench/bench-did.R
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

devtools::load_all()
data(lalonde.psid.panel, package = "qte")

cat("--- CiC (se=FALSE) ---\n")
print(system.time(
  CiC(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    se = FALSE, probs = seq(0.05, 0.95, 0.05)
  )
))

cat("--- CiC (se=TRUE, iters=100) ---\n")
print(system.time(
  CiC(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    se = TRUE, iters = 100, probs = seq(0.05, 0.95, 0.05)
  )
))

cat("--- QDiD (se=FALSE) ---\n")
print(system.time(
  QDiD(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    se = FALSE, probs = seq(0.05, 0.95, 0.05)
  )
))

cat("--- QDiD (se=TRUE, iters=100) ---\n")
print(system.time(
  QDiD(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    se = TRUE, iters = 100, probs = seq(0.05, 0.95, 0.05)
  )
))

cat("--- MDiD (se=FALSE) ---\n")
print(system.time(
  MDiD(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    se = FALSE, probs = seq(0.05, 0.95, 0.05)
  )
))

cat("--- MDiD (se=TRUE, iters=100) ---\n")
print(system.time(
  MDiD(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    se = TRUE, iters = 100, probs = seq(0.05, 0.95, 0.05)
  )
))

cat("--- ddid2 (se=FALSE) ---\n")
print(system.time(
  ddid2(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    se = FALSE, probs = seq(0.05, 0.95, 0.05)
  )
))

cat("--- ddid2 (se=TRUE, iters=100) ---\n")
print(system.time(
  ddid2(re ~ treat,
    t = 1978, tmin1 = 1975, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    se = TRUE, iters = 100, probs = seq(0.05, 0.95, 0.05)
  )
))

cat("--- panel.qtet (se=FALSE) ---\n")
print(system.time(
  suppressWarnings(panel.qtet(re ~ treat,
    t = 1978, tmin1 = 1975, tmin2 = 1974, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    se = FALSE, probs = seq(0.05, 0.95, 0.05)
  ))
))

cat("--- panel.qtet (se=TRUE, iters=100) ---\n")
print(system.time(
  suppressWarnings(panel.qtet(re ~ treat,
    t = 1978, tmin1 = 1975, tmin2 = 1974, tname = "year",
    data = lalonde.psid.panel, idname = "id",
    se = TRUE, iters = 100, probs = seq(0.05, 0.95, 0.05)
  ))
))
