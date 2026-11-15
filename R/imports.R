# =============================================================================
# Title: Package-level imports and documentation
# Description: Declares package-level @import directives and the package
#   documentation entry point. All broad namespace imports live here.
# Author: Brant Callaway
# Last update: 2026-11-15
# Date created: 2026-05-18
# =============================================================================

#' qte: A package for computing quantile treatment effects
#'
#' @import BMisc
#' @import graphics
#' @import stats
#' @import pbapply
#' @importFrom utils write.table tail
#' @import data.table
#' @importFrom ggplot2 autoplot
#' @importFrom rlang .data
"_PACKAGE"

# Suppress R CMD check notes for NSE variable references.
#   .w  — column in gt_data_frame used as weights in rq() calls
#   G, period — columns in ptetools data accessed via subset() NSE
#   Remaining names — injected into compute.unc_qte()'s environment by
#   setupData() via parent.frame(), plus eachIter in unc_qte(). Trimmed in
#   v2.1.0 after the deprecated two-period wrappers were removed.
utils::globalVariables(c(
  # ptetools / data.table column references
  ".w", "G", "period",
  # variables injected by setupData() (used by compute.unc_qte / unc_qte)
  "yname", "treat", "treated", "untreated", "wname", "probs", "method",
  "x", "xformla", "eachIter"
))
