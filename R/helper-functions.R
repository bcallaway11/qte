# =============================================================================
# Title: Helper functions for qte package
# Description: Utility functions for plotting QTE objects. Provides
#   autoplot.QTE for ggplot2-based plotting of QTE results (the deprecated
#   ggqte() was removed in v2.1.0). Also keep_used_columns(),
#   used by the staggered estimators before calling ptetools::pte().
# Author: Brant Callaway
# Last update: 2026-11-19
# Date created: 2026-05-18
# =============================================================================

# --- keep_used_columns -------------------------------------------------------
# Drop the columns of `data` that the estimators never read. ptetools copies
# the data for every (g,t) cell and every bootstrap draw, so unused columns
# only cost time. The *_gt functions use the outcome, group, time, id,
# weights and the covariates in xformula, plus columns ptetools adds itself.
keep_used_columns <- function(data, yname, gname, tname, idname = NULL,
                              weightsname = NULL, xformula = ~1) {
  data <- as.data.frame(data)
  used <- unique(c(yname, gname, tname, idname, weightsname, all.vars(xformula)))
  data[, intersect(used, names(data)), drop = FALSE]
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
