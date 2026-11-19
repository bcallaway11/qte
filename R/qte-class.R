# =============================================================================
# Title: QTE class and its methods
# Description: Defines the QTE S3 class (returned by unc_qte()) and its
#   print/summary/plot methods. The legacy computeSE()/SE() bootstrap
#   helpers were removed in v2.1.0 (unc_qte() has its own bootstrap), as was
#   the deprecated panel.qtet(). Dataset documentation is in data.R.
# Author: Brant Callaway
# Last update: 2026-11-19
# Date created: 2026-05-18
# =============================================================================

## summary function for QTE objects
#' @title Summary
#'
#' @description \code{summary.QTE} summarizes a \code{QTE} object, returning
#'   formatted data frames for the overall ATE and the QTE table suitable for
#'   printing.
#'
#' @param object A \code{QTE} object, as returned by \code{\link{unc_qte}}.
#' @param ... unused.
#'
#' @return A \code{summary.QTE} object (a list with \code{overall_ate} and
#'   \code{qte_table} data frames).
#' @export
summary.QTE <- function(object, ...) {
  alp  <- if (!is.null(object$alp)) object$alp else 0.05
  z    <- qnorm(1 - alp / 2)
  on_treated <- isTRUE(grepl("Treated", object$type, ignore.case = TRUE))

  # --- overall ATE/ATT ---
  ate    <- object$ate
  ate.se <- object$ate.se
  if (!is.null(ate.se) && !is.na(ate.se)) {
    ate.lower <- if (!is.null(object$ate.lower)) object$ate.lower else ate - z * ate.se
    ate.upper <- if (!is.null(object$ate.upper)) object$ate.upper else ate + z * ate.se
    ate_sig   <- !is.na(ate.upper) & ((ate.upper < 0) | (ate.lower > 0))
    out_ate   <- cbind.data.frame(
      round(ate,       4),
      round(ate.se,    4),
      round(ate.lower, 4),
      round(ate.upper, 4),
      ifelse(ate_sig, "*", "")
    )
  } else {
    out_ate <- cbind.data.frame(round(ate, 4))
  }

  # --- QTE table ---
  if (!is.null(object$qte)) {
    probs  <- object$probs
    qte    <- object$qte
    qte.se <- object$qte.se
    if (!is.null(qte.se)) {
      # prefer stored uniform band; fall back to pointwise from SE
      qte.lower <- if (!is.null(object$qte.lower)) object$qte.lower else qte - z * qte.se
      qte.upper <- if (!is.null(object$qte.upper)) object$qte.upper else qte + z * qte.se
      qte_sig   <- !is.na(qte.upper) & ((qte.upper < 0) | (qte.lower > 0))
      out_qte   <- cbind.data.frame(
        probs,
        round(qte,       4),
        round(qte.se,    4),
        round(qte.lower, 4),
        round(qte.upper, 4),
        ifelse(qte_sig, "*", "")
      )
    } else {
      out_qte <- cbind.data.frame(probs, round(qte, 4))
    }
  } else {
    out_qte <- NULL
  }

  out <- list(
    overall_ate = out_ate,
    qte_table   = out_qte,
    alp         = alp,
    on_treated  = on_treated,
    has_se      = !is.null(object$qte.se),
    # label the band "Simult." only when the object says it is one; objects
    # from qte < 2.1.0 have no `cband` field and stored pointwise intervals
    has_cband   = isTRUE(object$cband)
  )
  class(out) <- "summary.QTE"
  out
}

#' @title Print summary.QTE
#'
#' @description Prints a \code{summary.QTE} object.
#'
#' @param x A \code{summary.QTE} object.
#' @param ... unused.
#'
#' @return None. Called for its side effect of printing.
#' @export
print.summary.QTE <- function(x, ...) {
  alp        <- x$alp
  ci_header  <- paste0("[ ", 100 * (1 - alp), "% ")
  band_label <- if (x$has_cband) "Simult. " else "Pointwise "
  ate_label  <- if (x$on_treated) "Overall ATT" else "Overall ATE"
  qte_label  <- if (x$on_treated) "QTT" else "QTE"

  # --- overall ATE/ATT ---
  cat("\n")
  cat(paste0(ate_label, ":  \n"))
  out1 <- x$overall_ate
  if (x$has_se) {
    colnames(out1) <- c(
      if (x$on_treated) "ATT" else "ATE",
      "   Std. Error",
      paste0("    ", ci_header),
      "Conf. Int.]",
      ""
    )
  }
  print(out1, row.names = FALSE)
  cat("\n\n")

  # --- QTE table ---
  if (!is.null(x$qte_table)) {
    cat(paste0(qte_label, ":\n"))
    out2 <- x$qte_table
    if (x$has_se) {
      colnames(out2) <- c(
        "Tau",
        qte_label,
        "Std. Error",
        paste0(ci_header, band_label),
        "Conf. Band]",
        ""
      )
    } else {
      colnames(out2) <- c("Tau", qte_label)
    }
    print(out2, row.names = FALSE, justify = "centre")
    cat("---\n")
    cat("Signif. codes: `*' confidence band does not cover 0\n\n")
  }
}

#' @keywords internal
#' @noRd
qte_print_mat <- function(m, probs = NULL, header = NULL, digits = 2, nsmall = 2) {
  write.table(
    cbind(
      probs,
      format(m,
        justify = "right",
        digits = digits, nsmall = nsmall
      )
    ),
    row.names = FALSE, col.names = header, quote = FALSE, sep = "\t"
  )
  ## print(m, print.gap=3, right=T)
}

##
#' @title plot.QTE
#'
#' @description Plots a \code{QTE} object using \code{\link{autoplot.QTE}}.
#'
#' @param x a \code{QTE} object, as returned by \code{\link{unc_qte}}.
#' @param cband logical; if \code{TRUE} (default), show the uniform confidence
#'   band. If \code{FALSE}, show pointwise intervals.
#' @param ylab label for the y-axis. Default \code{"QTE"}.
#' @param ... passed to \code{\link{autoplot.QTE}}.
#'
#' @return invisibly returns the \code{ggplot} object.
#' @export
plot.QTE <- function(x, cband = TRUE, ylab = "QTE", ...) {
  p <- ggplot2::autoplot(x, cband = cband, ylab = ylab, ...)
  print(p)
  invisible(p)
}


##### SETUP CLASSES################
#' @title QTE
#'
#' @description Main class of objects.  A \code{QTE} object is returned by
#'  all of the methods that compute the QTE or QTET.
#'
#' @param qte The Quantile Treatment Effect at each value of probs
#' @param qte.se A vector of standard errors for each qte
#' @param qte.upper A vector of upper confidence intervals for each qte (it is
#'  based on the bootstrap confidence interval -- not the se -- so it may not
#'  be symmetric about the qte
#' @param qte.lower A vector of lower confidence intervals for each qte (it is
#'  based on the bootstrap confidence interval -- not the se -- so it may not
#'  be symmyetric about the qte
#' @param ate The Average Treatment Effect (or Average Treatment Effect on
#'  the Treated)
#' @param ate.se The standard error for the ATE
#' @param ate.lower Lower confidence interval for the ATE (it is based on the
#'  bootstrap confidence intervall -- not the se -- so it may not be symmetric
#'  about the ATE
#' @param ate.upper Upper confidence interval for the ATE (it is based on the
#'  bootstrap confidence interval -- not the se -- so it may not be symmetric
#'  about the ATE
#' @param c The critical value from a KS-type statistic used for creating
#'  uniform confidence bands
#' @param alp The significance level
#' @param cband Whether \code{qte.lower}/\code{qte.upper} are a uniform
#'  confidence band (\code{TRUE}) or pointwise intervals (\code{FALSE}).
#'  Default \code{NULL} (unknown, e.g. objects from qte < 2.1.0), treated as
#'  pointwise.
#' @param pscore.reg The results of propensity score regression, if specified
#' @param probs The values for which the qte is computed
#' @param type Takes the values "On the Treated" or "Population" to indicate
#'  whether the estimated QTE is for the treated group or for the entire
#'  population
#' @param F.treated.t Distribution of treated outcomes for the treated group at
#'  period t
#' @param F.untreated.t Distribution of untreated potential outcomes for the
#'  untreated group at period t
#' @param F.treated.t.cf Counterfactual distribution of untreated potential
#'  outcomes for the treated group at period t
#' @param F.treated.tmin1 Distribution of treated outcomes for the
#'  treated group at period tmin1
#' @param F.treated.tmin2 Distribution of treated outcomes for the
#'  treated group at period tmin2
#' @param F.treated.change.tmin1 Distribution of the change in outcomes for
#'  the treated group between periods tmin1 and tmin2
#' @param F.untreated.change.t Distribution of the change in outcomes for the
#'  untreated group between periods t and tmin1
#' @param F.untreated.change.tmin1 Distribution of the change in outcomes for
#'  the untreated group between periods tmin1 and tmin2
#' @param F.untreated.tmin1 Distribution of outcomes for the untreated group
#'  in period tmin1
#' @param F.untreated.tmin2 Distribution of outcomes for the untreated group
#'  in period tmin2
#' @param condQ.treated.t Conditional quantiles for the treated group in
#'  period t
#' @param condQ.treated.t.cf Counterfactual conditional quantiles for the treated
#'  group in period t
#' @param eachIterList An optional list of the outcome of each bootstrap
#'  iteration
#' @param inffunct The influence function for the treated group;
#'  used for inference when there are multiple
#'  periods and in the case with panel data.  It is needed for computing covariance
#'  terms in the variance-covariance matrix.
#' @param inffuncu The influence function for the untreated group
#'
#' @keywords internal
#' @export
QTE <- function(qte, ate = NULL, qte.se = NULL, qte.lower = NULL,
                qte.upper = NULL, ate.se = NULL, ate.lower = NULL, ate.upper = NULL,
                c = NULL, alp = 0.05, pscore.reg = NULL, probs, type = "On the Treated",
                F.treated.t = NULL, F.untreated.t = NULL, F.treated.t.cf = NULL,
                F.treated.tmin1 = NULL, F.treated.tmin2 = NULL,
                F.treated.change.tmin1 = NULL,
                F.untreated.change.t = NULL,
                F.untreated.change.tmin1 = NULL,
                F.untreated.tmin1 = NULL,
                F.untreated.tmin2 = NULL,
                condQ.treated.t = NULL,
                condQ.treated.t.cf = NULL,
                eachIterList = NULL, inffunct = NULL, inffuncu = NULL,
                cband = NULL) {
  out <- list(
    qte = qte, ate = ate, qte.se = qte.se, qte.lower = qte.lower,
    qte.upper = qte.upper, ate.se = ate.se, ate.lower = ate.lower,
    ate.upper = ate.upper, c = c, alp = alp, cband = cband,
    pscore.reg = pscore.reg, probs = probs,
    type = type, F.treated.t = F.treated.t, F.untreated.t = F.untreated.t,
    F.treated.t.cf = F.treated.t.cf,
    F.treated.tmin1 = F.treated.tmin1,
    F.treated.tmin2 = F.treated.tmin2,
    F.treated.change.tmin1 = F.treated.change.tmin1,
    F.untreated.change.t = F.untreated.change.t,
    F.untreated.change.tmin1 = F.untreated.change.tmin1,
    F.untreated.tmin1 = F.untreated.tmin1,
    F.untreated.tmin2 = F.untreated.tmin2,
    condQ.treated.t = condQ.treated.t,
    condQ.treated.t.cf = condQ.treated.t.cf,
    eachIterList = eachIterList,
    inffunct = inffunct,
    inffuncu = inffuncu
  )
  class(out) <- "QTE"
  return(out)
}
