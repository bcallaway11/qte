# =============================================================================
# Title: QTE class and its methods
# Description: Defines the QTE and SE S3 classes and their print/summary/plot
#   methods, and computeSE() for the legacy empirical bootstrap used by
#   unc_qte(). (The deprecated panel.qtet() was removed in v2.1.0; use
#   panel_qtt().) Dataset documentation is in data.R.
# Author: Brant Callaway
# Last update: 2026-11-15
# Date created: 2026-05-18
# =============================================================================

#' @title computeSE
#'
#' @description Computes standard errors from bootstrap results.  This function
#'  is called by several functions in the qte package
#'
#' @param bootIters List of bootstrap iterations
#' @param qteobj the \code{QTE} object holding the full-sample estimates
#' @param alp significance level for confidence intervals
#'
#' @keywords internal
#'
#' @return SEObj
computeSE <- function(bootIters, qteobj, alp = 0.05) {
  ## For now, just plot the qte and att with standard errors
  ## helper function to get the first element out of a list
  qte <- qteobj$qte
  ate <- qteobj$ate
  iters <- length(bootIters)

  getElement <- function(Lst, elemNum) {
    return(as.numeric(unlist((Lst[elemNum])))) # as.numeric is a trick to
    ## get numerical value of qte
  }
  all.ate <- unlist(sapply(bootIters, FUN = getElement, elemNum = 2))
  ## get se
  ate.se <- sd(all.ate)
  ate.upper <- ate + qnorm(1 - alp / 2) * ate.se
  ate.lower <- ate - qnorm(1 - alp / 2) * ate.se
  ## reorder asc
  ## all.ate = all.ate[order(all.ate)]
  ## ate.upper = all.ate[min(iters,round((1-alp/2)*iters))]
  ## ate.lower = all.ate[max(1,round((alp/2)*iters))]

  ## now get CI for qte:
  all.qte <- lapply(bootIters, FUN = getElement, elemNum = 1)
  qte.mat <- do.call(rbind, lapply(all.qte, FUN = as.numeric, ncol = length(all.qte[[1]]), byrow = TRUE))
  ## standard error
  qte.se <- apply(qte.mat, FUN = sd, MARGIN = 2)

  sigmahalf <- (apply(qte.mat, 2, function(b) quantile(b, .75, type = 1)) -
    apply(qte.mat, 2, function(b) quantile(b, .25, type = 1))) / (qnorm(.75) - qnorm(.25))

  ## this seems to work a bit better in practice when some QTEs are 0
  if (any(sigmahalf == 0)) {
    sigmahalf <- apply(qte.mat, 2, sd)
    sigmahalf <- sapply(seq_along(qte), function(i) max(sigmahalf[i], .000000001))
  }
  cb <- apply(qte.mat, 1, function(q) max(abs((q - qte) / sigmahalf)))
  c <- quantile(cb, (1 - alp), type = 1)
  ## qte se by quantiles
  ## sorted.qtemat = apply(qte.mat, 2, sort)
  ## qte.upper = sorted.qtemat[round((1-alp/2)*iters),]
  ## qte.lower = sorted.qtemat[max(1,round((alp/2)*iters)),]
  ## qte se by sd
  qte.upper <- qte + qnorm(1 - alp / 2) * qte.se
  qte.lower <- qte - qnorm(1 - alp / 2) * qte.se

  out <- SE(
    ate.se = ate.se, ate.upper = ate.upper, ate.lower = ate.lower,
    qte.se = qte.se, qte.upper = qte.upper, qte.lower = qte.lower,
    c = c
  )
  return(out)
}


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
    has_cband   = !is.null(object$qte.upper)
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
                eachIterList = NULL, inffunct = NULL, inffuncu = NULL) {
  out <- list(
    qte = qte, ate = ate, qte.se = qte.se, qte.lower = qte.lower,
    qte.upper = qte.upper, ate.se = ate.se, ate.lower = ate.lower,
    ate.upper = ate.upper, c = c, alp = alp,
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

#' @title SE
#'
#' @description Class for Standard Error Objects
#'
#' @param qte.se The QTE Standard Error
#' @param ate.se The ATE Standard Error
#' @param qte.upper The QTE upper CI
#' @param qte.lower The QTE lower CI
#' @param ate.upper The ATE upper CI
#' @param ate.lower The ATE lower CI
#' @param c The critical value from a KS-type statistic used for creating
#'  uniform confidence bands
#' @param probs The values at which the QTE is computed
#'
#' @keywords internal
SE <- function(qte.se = NULL, ate.se = NULL, qte.upper = NULL, qte.lower = NULL,
               ate.upper = NULL, ate.lower = NULL, c = NULL, probs = NULL) {
  out <- list(
    qte.se = qte.se, qte.upper = qte.upper, qte.lower = qte.lower,
    ate.se = ate.se, ate.upper = ate.upper, ate.lower = ate.lower,
    c = c, probs = probs
  )
  class(out) <- "SE"
  return(out)
}
