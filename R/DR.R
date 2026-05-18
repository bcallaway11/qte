# =============================================================================
# Title: Distribution Regression (DR) class
# Description: Constructor for the DR S3 class, which holds a fitted
#   distribution regression object (a list of models, one per y-value).
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

#' @title DR
#'
#' @description A distribution regression object
#'
#' @param yvals A vector of values that y can take
#' @param drlist A list where for each value of y, a distribution regression
#'
#' @export
DR <- function(yvals, drlist) {
  out <- list(yvals = yvals, drlist = drlist)
  class(out) <- "DR"
  return(out)
}
