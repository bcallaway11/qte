# =============================================================================
# Title: Test helper for the discrete-outcome warning
# Description: quiet_discrete() muffles only the qte_discrete_outcome warning
#   (#17), for tests on lalonde.psid (12% zero earnings) that check other
#   things; other warnings still reach expect_warning() etc.
# Author: Brant Callaway
# Last update: 2026-11-19
# Date created: 2026-11-19
# =============================================================================

quiet_discrete <- function(expr) {
  withCallingHandlers(
    expr,
    qte_discrete_outcome = function(w) invokeRestart("muffleWarning")
  )
}
