# =============================================================================
# Title: Helper functions for qte package
# Description: Utility functions for plotting QTE objects.
# Author: Brant Callaway
# Last update: 2026-05-21
# Date created: 2026-05-18
# =============================================================================

#' @title ggqte
#'
#' @description Makes somewhat nicer plots of quantile treatment effects
#'  by using ggplot
#'
#' @import ggplot2
#'
#' @param qteobj a QTE object
#' @param main optional title
#' @param ylab optional y axis label
#' @param ylim optional limits of y axis
#' @param ybreaks optional breaks in y axis
#' @param xbreaks optional breaks in x axis
#' @param setype options are "pointwise", "uniform" or both; pointwise confidence
#'  intervals cover the QTE at each particular point with a fixed probability,
#'  uniform confidence bands cover the entire curve with a fixed
#'  probability.  Uniform confidence bands will tend to be wider.  The option
#'  "both" will plot both types of confidence intervals
#' @param alp gives a way to override the significance level in the case where
#'  `setype="pointwise"`.
#'
#' @return a ggplot object
#' @export
ggqte <- function(qteobj, main = "", ylab = "QTE", ylim = NULL, ybreaks = NULL, xbreaks = c(.1, .3, .5, .7, .9), setype = "pointwise", alp = qteobj$alp) {
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

