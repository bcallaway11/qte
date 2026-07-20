# =============================================================================
# Title: Run revdeplite
# Description: Run lightweight reverse dependency checks for qte.
# Author: Brant Callaway
# Last update: 2026-07-20
# Date created: 2026-07-20
# =============================================================================

revdeplite::revdeplite(
  github_deps = c(
    "bcallaway11/csabounds"
  ),
  check_dir = ".revdeplite"
)
