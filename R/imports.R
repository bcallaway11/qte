# =============================================================================
# Title: Package-level imports and documentation
# Description: Declares package-level @import directives and the package
#   documentation entry point. All broad namespace imports live here.
# Author: Brant Callaway
# Last update: 2026-11-19
# Date created: 2026-05-18
# =============================================================================

#' qte: A package for computing quantile treatment effects
#'
#' @import BMisc
#' @import graphics
#' @import stats
#' @import pbapply
#' @importFrom utils write.table tail
#' @importFrom ggplot2 autoplot
#' @importFrom rlang .data
"_PACKAGE"

# Suppress the R CMD check note for the one NSE column reference left:
#   .w — the sampling-weights column, passed as `weights = .w` to rq() and
#   glm() so it is looked up in the data. (The names injected by the legacy
#   setupData(), and the G/period ptetools columns, went away in v2.1.0.)
utils::globalVariables(".w")
