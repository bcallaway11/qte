# =============================================================================
# Title: Helper functions for qte package
# Description: Utility functions for plotting QTE objects. Provides
#   autoplot.QTE for ggplot2-based plotting of QTE results (the deprecated
#   ggqte() was removed in v2.1.0). Also keep_used_columns(),
#   used by the staggered estimators before calling ptetools::pte().
# Author: Brant Callaway
# Last update: 2026-11-20
# Date created: 2026-05-18
# =============================================================================

# --- keep_used_columns -------------------------------------------------------
# Drop the columns of `data` that the estimators never read. ptetools copies
# the data for every (g,t) cell and every bootstrap draw, so unused columns
# only cost time. The *_gt functions use the outcome, group, time, id,
# weights, the cluster variable and the covariates in xformula, plus columns
# ptetools adds itself.
keep_used_columns <- function(data, yname, gname, tname, idname = NULL,
                              weightsname = NULL, xformula = ~1,
                              clustervars = NULL) {
  data <- as.data.frame(data)
  used <- unique(c(yname, gname, tname, idname, weightsname, clustervars,
                   all.vars(xformula)))
  data[, intersect(used, names(data)), drop = FALSE]
}


# --- base_period_covariates --------------------------------------------------
# Give every row of a (g,t) cell its unit's covariate values from one
# reference period: the base period ("pre") for the two-period cells, the
# earliest period used ("pre2") for panel_qtt's three-period cells.
# Conditioning on covariates from a single pre-treatment period is the
# convention throughout (as in ptetools::pte_attgt() and did): the
# comparison group's post-period distribution is conditioned on the same
# covariate values that are plugged in for the treated, and the treated are
# never conditioned on covariates measured after treatment. Each *_gt
# function calls this itself, so a custom attgt_fun can choose differently.
#
# Repeated cross sections are untouched: a row matches a reference-period
# row only if the same id appears there, which happens only for panel data.
# Panel units with no reference-period row keep their own values.
#
# Args: gt_data (gt_data_frame for one cell); xformula; ref (value of
#   gt_data$name for the reference period).
# Returns: gt_data with the covariate columns replaced.
base_period_covariates <- function(gt_data, xformula, ref = "pre") {
  xvars <- intersect(all.vars(xformula), names(gt_data))
  if (length(xvars) == 0L) return(gt_data)
  ref_rows <- gt_data[gt_data$name == ref, , drop = FALSE]
  idx <- match(gt_data$id, ref_rows$id)
  has_ref <- !is.na(idx)
  for (v in xvars) gt_data[[v]][has_ref] <- ref_rows[[v]][idx[has_ref]]
  gt_data
}


# --- autoplot.QTE ------------------------------------------------------------

#' @title autoplot.QTE
#'
#' @description Plot a \code{QTE} object as a quantile treatment effect curve
#'   with optional confidence bands.
#'
#' @param object a \code{QTE} object, as returned by \code{\link{unc_qte}}.
#' @param cband logical; if \code{TRUE} (default), show the uniform confidence
#'   band, \code{qte} \eqn{\pm} \code{c} \eqn{\times} \code{qte.se} (as
#'   stored by \code{\link{unc_qte}} with \code{cband = TRUE}). If
#'   \code{FALSE}, show pointwise intervals computed from \code{object$qte.se}.
#' @param ylab label for the y-axis. Default \code{"QTE"}.
#' @param ... unused.
#'
#' @return a \code{ggplot} object.
#' @method autoplot QTE
#' @export
autoplot.QTE <- function(object, cband = TRUE, ylab = "QTE", ...) {
  probs <- object$probs
  qte   <- object$qte
  alp   <- if (is.null(object$alp)) 0.05 else object$alp

  df <- data.frame(probs = probs, qte = qte)

  has_se   <- !is.null(object$qte.se) && !anyNA(object$qte.se)
  has_band <- !is.null(object$qte.upper) && !anyNA(object$qte.upper)
  has_cval <- !is.null(object$c) && !anyNA(object$c)

  if (cband && isTRUE(object$cband) && has_band) {
    # stored bounds are already the uniform band (qte +/- c * se)
    df$lower <- object$qte.lower
    df$upper <- object$qte.upper
  } else if (cband && has_se && has_cval) {
    # objects from qte < 2.1.0 stored pointwise bounds but kept the uniform
    # critical value c, so the band can be rebuilt the same way
    cval     <- as.numeric(object$c)
    df$lower <- qte - cval * object$qte.se
    df$upper <- qte + cval * object$qte.se
  } else if (has_se) {
    cval     <- qnorm(1 - alp / 2)
    df$lower <- qte - cval * object$qte.se
    df$upper <- qte + cval * object$qte.se
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = probs, y = qte)) +
    ggplot2::geom_hline(yintercept = 0, color = "gray50") +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous("Quantile", limits = c(0, 1)) +
    ggplot2::scale_y_continuous(ylab) +
    ggplot2::theme_bw()

  if ("lower" %in% names(df)) {
    p <- p +
      ggplot2::geom_line(ggplot2::aes(y = .data$lower), linetype = "dashed") +
      ggplot2::geom_line(ggplot2::aes(y = .data$upper), linetype = "dashed")
  }

  p
}
