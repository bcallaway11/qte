# =============================================================================
# Title: Helper functions for qte package
# Description: Utility functions for plotting QTE objects. Provides
#   autoplot.QTE and plot.QTE for ggplot2-based plotting of QTE results,
#   and retains ggqte() as a deprecated wrapper.
# Author: Brant Callaway
# Last update: 2026-05-22
# Date created: 2026-05-18
# =============================================================================


# --- autoplot.QTE ------------------------------------------------------------

#' @title autoplot.QTE
#'
#' @description Plot a \code{QTE} object as a quantile treatment effect curve
#'   with optional confidence bands.
#'
#' @param object a \code{QTE} object, as returned by \code{\link{unc_qte}}.
#' @param cband logical; if \code{TRUE} (default), show the uniform confidence
#'   band stored in \code{object$qte.upper} / \code{object$qte.lower}.
#'   If \code{FALSE}, show pointwise intervals computed from \code{object$qte.se}.
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

  if (has_band && cband) {
    df$lower <- object$qte.lower
    df$upper <- object$qte.upper
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


# --- ggqte (deprecated) ------------------------------------------------------

#' @title ggqte
#'
#' @description \strong{Deprecated.} Use \code{autoplot()} on a \code{QTE}
#'   object instead.
#'
#' @param qteobj a QTE object
#' @param main optional title
#' @param ylab optional y axis label
#' @param ylim optional limits of y axis
#' @param ybreaks optional breaks in y axis
#' @param xbreaks optional breaks in x axis
#' @param setype options are "pointwise", "uniform" or both
#' @param alp significance level override
#'
#' @return a ggplot object
#' @export
ggqte <- function(qteobj, main = "", ylab = "QTE", ylim = NULL, ybreaks = NULL, xbreaks = c(.1, .3, .5, .7, .9), setype = "pointwise", alp = qteobj$alp) {
  .Deprecated(msg = paste0(
    "ggqte() is deprecated. Use autoplot() on a QTE object instead:\n",
    "  autoplot(result)              # uniform band (default)\n",
    "  autoplot(result, cband=FALSE) # pointwise intervals"
  ))
  tau <- qteobj$probs
  qte <- qteobj$qte
  qte.se <- qteobj$qte.se
  c <- qteobj$c
  if (is.null(alp)) alp <- 0.05
  if (!is.null(qte.se)) {
    cmat <- data.frame(tau, qte = qteobj$qte, qte.se = qteobj$qte.se)
  } else {
    cmat <- data.frame(tau, qte = qteobj$qte)
  }
  qp <- ggplot2::ggplot(data = cmat, aes(tau, qte)) +
    ggplot2::geom_line(aes(tau, qte)) +
    ## geom_errorbar(size=.3, width=.02) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_point(aes(tau, qte)) +
    ggplot2::ggtitle(main) +
    ggplot2::scale_x_continuous("tau", limits = c(0, 1), breaks = xbreaks) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.border = element_rect(
        colour = "black", size = 1,
        fill = NA,
        linetype = "solid"
      ),
      plot.title = element_text(hjust = 0.5)
    )
  if (is.null(ylim) && is.null(ybreaks)) {
    qp <- qp + ggplot2::scale_y_continuous(ylab)
  } else if (is.null(ylim) && !is.null(ybreaks)) {
    qp <- qp + ggplot2::scale_y_continuous(ylab, breaks = ybreaks)
  } else if (!is.null(ylim) && is.null(ybreaks)) {
    qp <- qp + ggplot2::scale_y_continuous(ylab, limits = ylim)
  } else {
    qp <- qp + ggplot2::scale_y_continuous(ylab, limits = ylim, breaks = ybreaks)
  }

  if (!is.null(qte.se)) {
    if (setype == "both" || setype == "pointwise") {
      qp <- qp + ggplot2::geom_line(aes(tau, qte + qnorm(1 - alp / 2) * qte.se), linetype = "dashed")
      qp <- qp + ggplot2::geom_line(aes(tau, qte - qnorm(1 - alp / 2) * qte.se), linetype = "dashed")
    }
    if (setype == "both" || setype == "uniform") {
      qp <- qp + ggplot2::geom_line(aes(tau, qte + c * qte.se), linetype = "dashed") + ggplot2::geom_line(aes(tau, qte - c * qte.se), linetype = "dashed")
    }
  }
  qp
}

