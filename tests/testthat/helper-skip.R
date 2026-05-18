# =============================================================================
# Title: Test skip helpers
# Description: Provides skip_slow() to gate slow tests behind an env var.
#   Run slow tests with: R_SLOW_TESTS=true Rscript -e "devtools::test()"
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

skip_slow <- function() {
  if (!identical(Sys.getenv("R_SLOW_TESTS"), "true")) {
    testthat::skip("slow test — set R_SLOW_TESTS=true to run")
  }
}
