# =============================================================================
# Title: Package-level imports and documentation
# Description: Declares package-level @import directives and the package
#   documentation entry point. All broad namespace imports live here.
# Author: Brant Callaway
# Last update: 2026-05-22
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
# TODO (v2.1): trim this list after deprecated wrappers and ggqte() are removed.
#   .w  — column in gt_data_frame used as weights in rq() calls
#   G, period — columns in ptetools data accessed via subset() NSE
#   aes, element_rect, element_text — ggplot2 NSE in deprecated ggqte()
#   Remaining names — injected into caller envs by setupData() via parent.frame()
utils::globalVariables(c(
  # ptetools / data.table column references
  ".w", "G", "period",
  # ggplot2 references in deprecated ggqte()
  "aes", "element_rect", "element_text",
  # variables injected by setupData() via assign(..., envir = parent.frame())
  "yname", "treat", "panel", "treated.t", "treated.tmin1",
  "untreated.t", "untreated.tmin1",
  "F.treated.t", "F.treated.tmin1", "F.untreated.t", "F.untreated.tmin1",
  "xformla", "data", "x", "wname", "probs", "method", "eachIter",
  "F.untreated.tmin2", "F.treated.tmin2", "untreated.tmin2", "treated.tmin2",
  "untreated.change.t", "F.untreated.change.t",
  "untreated.change.tmin1", "F.untreated.change.tmin1",
  "treated.change.tmin1", "F.treated.change.tmin1",
  "tmin1", "tmin2", "idname"
))
