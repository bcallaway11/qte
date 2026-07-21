# =============================================================================
# Title: QTE class and panel QTET estimator
# Description: Defines the QTE and SE S3 classes and their print/summary/plot
#   methods. Also contains the panel.qtet estimator (Callaway-Li 2019) and its
#   internal compute function. Dataset documentation is in data.R.
# Author: Brant Callaway
# Last update: 2026-07-20
# Date created: 2026-05-18
# =============================================================================

##### MAIN FUNCTIONS#####

##### Panel QTET#####
## Idea here is that we can use information from a third period
## to point identify counterfactual distribution of outcomes
## for the treated group
## call plot function, summary function, formula function, etc. later
## add functionality to pass in pscore
#' @title compute.panel.qtet
#'
#' @description
#' \code{compute.panel.qtet} uses third period of data,
#' combined with Distributional
#' Difference in Differences assumption (Fan and Yu, 2012)
#' to point identify QTET.
#'
#' @inheritParams compute.unc_qte
#'
#' @return QTE object
#'
#' @keywords internal
compute.panel.qtet <- function(qp) {
  setupData(qp)


  ## calculate the distribution of the change for the treated group;
  ## this will be changed if there are covariates
  F.treated.change.t <- F.untreated.change.t


  ## calculate the att; this will be changed if there are covariates
  att <- mean(treated.t[, yname]) - mean(treated.tmin1[, yname]) -
    (mean(untreated.t[, yname]) - mean(untreated.tmin1[, yname]))


  ## build counterfactual distribution depending on which case we are in

  ## a.1) If there are covariates need to satisfy the Distributional D-i-D
  ## then we will need to modify the distribution of the changes in outcomes
  ## using the method presented in the paper.

  pscore.reg <- NULL # do this in case no covariates as we return this value
  qr.reg <- NULL

  ## setup the data some
  if (!(is.null(x))) {
    treated.t$dy <- treated.t[, yname] - treated.tmin1[, yname]
    treated.tmin1$dy <- treated.tmin1[, yname] - treated.tmin2[, yname]
    untreated.t$dy <- untreated.t[, yname] - untreated.tmin1[, yname]
    untreated.tmin1$dy <- untreated.tmin1[, yname] -
      untreated.tmin2[, yname]
  }

  if (is.null(x)) method <- "pscore" ## just trick to reuse some code in pscore and no covariates case

  if (method == "pscore") {
    if (!(is.null(x))) {
      ## estimate pscore
      this.formla <- y ~ x
      formula.tools::rhs(this.formla) <- formula.tools::rhs(xformla)
      formula.tools::lhs(this.formla) <- as.name(treat)
      pscore.data <- rbind(treated.tmin2, untreated.tmin2)
      pscore.reg <- glm(this.formla,
        data = pscore.data,
        family = binomial(link = "logit")
      )
      pscore <- fitted(pscore.reg)


      ## calculate ATT and QTT
      dtat <- rbind.data.frame(treated.t, untreated.t)
      dtat$pscore <- pscore
      pD1 <- nrow(treated.t) / nrow(untreated.t)

      ## this contains the support of the change in y
      p.dy.seq <- unique(dtat$dy)
      p.dy.seq <- p.dy.seq[order(p.dy.seq)]
      distvals <- rep(0, length(p.dy.seq))
      for (i in seq_along(p.dy.seq)) {
        distvals[i] <- getWeightedMean(1 * (dtat$dy <= p.dy.seq[i]), (1 - dtat[, treat]) * pscore / ((1 - pscore) * pD1), norm = TRUE)
        ## mean(1*(pscore.data$changey <= p.dy.seq[i])*
        ##             (1-pscore.data[,treat])*pscore/((1-pscore)*pD1)) /
        ## mean( (1-pscore.data[,treat])*pscore/((1-pscore)*pD1) )
      }

      F.untreated.change.t <- make_dist(p.dy.seq, distvals)

      ## after we have the propensity score
      ## use it to estimate the att using abadie's method.
      att <- getWeightedMean(
        y = dtat$dy,
        weights = (dtat[, treat] - pscore) / ((1 - pscore) * pD1), norm = TRUE
      ) ## mean(((pscore.data$changey)/pval)*(pscore.data[,treat] - pscore) /
      ##      (1-pscore))

      ## update the lag of the untreated change so that we can
      ## do pre-testing if desired
      pscore.data.tmin1 <- rbind(treated.tmin1, untreated.tmin1)
      posvals.seq <- pscore.data.tmin1$dy
      distvals.tmin1 <- rep(0, length(posvals.seq))
      for (dy in posvals.seq) {
        distvals.tmin1[which(dy == posvals.seq)] <-
          mean(1 * (pscore.data.tmin1$dy <= dy) *
            (1 - pscore.data.tmin1[, treat]) * pscore / ((1 - pscore) * pD1))
      }
      pscore.data.tmin1$distvals <- distvals.tmin1
      pscore.data1.tmin1 <- pscore.data.tmin1[order(pscore.data.tmin1$dy), ]
      F.untreated.change.tmin1 <- make_dist(
        pscore.data1.tmin1$dy,
        pscore.data1.tmin1$distvals
      )
    }


    ## compute counterfactual distribution
    quantys1 <- quantile(F.treated.tmin1,
      probs = F.treated.tmin2(treated.tmin2[, yname])
    )

    quantys2 <- quantile(F.untreated.change.t,
      probs = F.treated.change.tmin1(treated.tmin1[, yname] -
        treated.tmin2[, yname])
    )

    y.seq <- (quantys1 + quantys2)[order(quantys1 + quantys2)]

    F.treated.t.cf.val <- vapply(y.seq,
      FUN = function(y) {
        mean(1 * (quantys2 <=
          y - quantys1))
      }, FUN.VALUE = 1
    )

    F.treated.t.cf <- make_dist(y.seq, F.treated.t.cf.val)
  } else if (method == "qr") {
    u <- seq(.01, .99, .01) ## hard-coded for now
    yformla <- BMisc::toformula("y", BMisc::rhs_vars(xformla))
    dyformla <- BMisc::toformula("dy", BMisc::rhs_vars(xformla))

    dQRt <- quantreg::rq(dyformla, data = untreated.t, tau = u) ## holds by conditional did assumption
    dQRtmin1 <- quantreg::rq(dyformla, data = treated.tmin1, tau = u)
    QRtmin1 <- quantreg::rq(yformla, data = treated.tmin1, tau = u)
    QRtmin2 <- quantreg::rq(yformla, data = treated.tmin2, tau = u)

    ## use csa-type result; exploit that we average over X_i
    n1 <- nrow(treated.t)
    n0 <- nrow(untreated.t)
    QRtmin2F <- predict(QRtmin2, type = "Fhat", stepfun = TRUE)
    Ftmin2 <- sapply(1:n1, function(i) QRtmin2F[[i]](treated.tmin2$y[i]))
    QRtmin1Q <- predict(QRtmin1, type = "Qhat", stepfun = TRUE)
    Qtmin1 <- sapply(1:n1, function(i) QRtmin1Q[[i]](Ftmin2[i]))

    dQRtmin1F <- predict(dQRtmin1, type = "Fhat", stepfun = TRUE)
    dFtmin1 <- sapply(1:n1, function(i) dQRtmin1F[[i]](treated.tmin1$dy[i]))
    dQRtQ <- predict(dQRt, newdata = treated.t, type = "Qhat", stepfun = TRUE) ## predict for the treated guys even though estimate with the untreated groups

    dQt <- sapply(1:n1, function(i) dQRtQ[[i]](dFtmin1[i]))

    yvals <- unique(rbind.data.frame(treated.t, untreated.t)$y)
    yvals <- sort(yvals)
    F.treated.t.cf.val <- sapply(yvals, function(yy) mean(1 * (dQt + Qtmin1 <= yy)))

    F.treated.t.cf <- make_dist(yvals, F.treated.t.cf.val)

    att <- mean(treated.t$y) - sum(quantile(F.treated.t.cf, probs = u, type = 1)) / length(u)
  } else {
    stop("invalid method supplied")
  }


  ## QTE
  F.treated.t <- ecdf(treated.t[, yname])

  qte <- quantile(F.treated.t, probs = probs) -
    quantile(F.treated.t.cf, probs = probs)


  out <- QTE(
    F.treated.t = F.treated.t,
    F.treated.tmin1 = F.treated.tmin1,
    F.treated.tmin2 = F.treated.tmin2,
    F.treated.change.tmin1 = F.treated.change.tmin1,
    F.untreated.t = F.untreated.t,
    F.untreated.tmin1 = F.untreated.tmin1,
    F.untreated.tmin2 = F.untreated.tmin2,
    F.untreated.change.t = F.untreated.change.t,
    F.untreated.change.tmin1 = F.untreated.change.tmin1,
    F.treated.t.cf = F.treated.t.cf,
    qte = qte, pscore.reg = pscore.reg, ate = att, probs = probs
  )
  class(out) <- "QTE"
  return(out)
}


#' @title panel.qtet
#'
#' @description \strong{Deprecated.} Use \code{\link{panel_qtt}} instead.
#'   \code{panel_qtt} supports staggered treatment adoption and exposes the
#'   \code{pre_copula} option for controlling how the copula is transferred
#'   across pre-treatment periods.
#'
#'   \code{panel.qtet} computes the Quantile Treatment Effect
#' on the Treated (QTET) using the method of Callaway and Li (2015).  This
#' method should be used when the researcher wants to invoke a Difference
#' in Differences assumption to identify the QTET.  Relative to the other
#' Difference in Differences methods available in the \code{qte} package,
#' this method's assumptions are more intuitively similar to the identifying
#' assumptions used in identifying the Average Treatment Effect on the Treated
#' (ATT).
#'
#' Additionally, this method can accommodate covariates in a more
#' flexible way than the other Difference in Differences methods available.
#' In order to accommodate covariates, the user should specify a vector \code{x}
#' of covariate names.  The user also may specify a method for estimating
#' the propensity score.  The default is logit.
#'
#' \code{panel.qtet} can only be used in some situations, however.  The
#' method requires three periods of panel data where individuals
#' are not treated until the last period.  The data should be formatted
#' as a panel; the names of columns containing time periods and ids
#' for each cross sectional unit need to be passed to the method.
#'
#' @param formla The formula y ~ d where y is the outcome and d is the
#'  treatment indicator (d should be binary), d should be equal to one
#'  in all time periods for individuals that are eventually treated
#' @param xformla A optional one sided formula for additional covariates that
#'  will be adjusted for.  E.g ~ age + education.  Additional covariates can
#'  also be passed by name using the x paramater.
#' @param t The 3rd time period in the sample.  Treated individuals should
#'  be treated in this time period and untreated individuals should not be
#'  treated.  The code attempts to enforce this condition, but it is good
#'  try to handle this outside the panel.qtet method.
#' @param tmin1 The 2nd time period in the sample.  This should be a
#'  pre-treatment period for all individuals in the sample.
#' @param tmin2 The 1st time period in the sample.  This should be a
#'  pre-treatment period for all individuals in the sample.
#' @param tname The name of the column containing the time periods
#' @param data A data.frame containing all the variables used
#' @param idname The individual (cross-sectional unit) id name
#' @param probs A vector of values between 0 and 1 to compute the QTET at
#' @param iters The number of iterations to compute bootstrap standard errors.
#'  This is only used if se=TRUE
#' @param alp The significance level used for constructing bootstrap
#'  confidence intervals
#' @param method The method for including covariates, should either be "QR"
#'  for quantile regression or "pscore" for propensity score
#' @param se Boolean whether or not to compute standard errors
#' @param retEachIter Boolean whether or not to return list of results
#'  from each iteration of the bootstrap procedure (default is FALSE).
#'  This is potentially useful for debugging but can cause errors due
#'  to running out of memory.
#' @param pl Whether or not to compute standard errors in parallel
#' @param cores Number of cores to use if computing in parallel
#' @param biters Number of bootstrap iterations; alias for \code{iters}
#'   matching the \code{did}/\code{ptetools} naming convention. If both are
#'   supplied, \code{biters} takes precedence.
#' @param cl Number of cores for parallel bootstrap; alias for
#'   \code{pl}/\code{cores}. \code{cl = 1} (default) runs sequentially;
#'   \code{cl > 1} enables parallelism.
#'
#' @examples
#' # See ?panel_qtt for the modern replacement.
#' \dontrun{
#' data(lalonde)
#'
#' ## Run the panel.qtet method on the experimental data with no covariates
#' pq1 <- panel.qtet(re ~ treat,
#'   t = 1978, tmin1 = 1975, tmin2 = 1974, tname = "year",
#'   data = lalonde.exp.panel, idname = "id", se = FALSE,
#'   probs = seq(0.05, 0.95, 0.05)
#' )
#' summary(pq1)
#'
#' ## Run the panel.qtet method on the observational data with no covariates
#' pq2 <- panel.qtet(re ~ treat,
#'   t = 1978, tmin1 = 1975, tmin2 = 1974, tname = "year",
#'   data = lalonde.psid.panel, idname = "id", se = FALSE,
#'   probs = seq(0.05, 0.95, 0.05)
#' )
#' summary(pq2)
#'
#' ## Run the panel.qtet method on the observational data conditioning on
#' ## age, education, black, hispanic, married, and nodegree.
#' ## The propensity score will be estimated using the default logit method.
#' pq3 <- panel.qtet(re ~ treat,
#'   t = 1978, tmin1 = 1975, tmin2 = 1974, tname = "year",
#'   xformla = ~ age + I(age^2) + education + black + hispanic + married + nodegree,
#'   data = lalonde.psid.panel, idname = "id", se = FALSE, method = "pscore",
#'   probs = seq(0.05, 0.95, 0.05)
#' )
#' summary(pq3)
#'
#' pq4 <- panel.qtet(re ~ treat,
#'   t = 1978, tmin1 = 1975, tmin2 = 1974, tname = "year",
#'   xformla = ~ age + I(age^2) + education + black + hispanic + married + nodegree,
#'   data = lalonde.psid.panel, idname = "id", se = FALSE, method = "qr",
#'   probs = seq(0.05, 0.95, 0.05)
#' )
#' summary(pq4)
#' }
#'
#' @references
#' Callaway, Brantly and Tong Li.  ``Quantile Treatment Effects in Difference
#'  in Differences Models with Panel Data.'' Working Paper, 2019.
#'
#' @seealso \code{\link{panel_qtt}} for the modern replacement supporting
#'   staggered adoption.
#'
#' @return \code{QTE} object
#'
#' @keywords internal
#' @export
panel.qtet <- function(formla, xformla = NULL, t, tmin1, tmin2,
                       tname, data,
                       idname, probs = seq(0.05, 0.95, 0.05),
                       iters = 100, alp = 0.05, method = c("qr", "pscore"), se = TRUE,
                       retEachIter = FALSE, pl = FALSE, cores = NULL,
                       biters = NULL, cl = NULL) {
  .Deprecated(msg = paste0(
    "panel.qtet() is deprecated. Use panel_qtt() instead.\n",
    "panel_qtt() supports staggered treatment adoption and uses a ",
    "yname/gname/tname interface. See ?panel_qtt."
  ))
  # biters and cl are aliases matching the did/ptetools API
  if (!is.null(biters)) iters <- biters
  if (!is.null(cl)) {
    pl    <- cl > 1
    cores <- cl
  }
  method <- method[1]

  data <- panelize.data(data, idname, tname, t, tmin1, tmin2)

  qp <- QTEparams(
    formla = formla, xformla = xformla,
    t = t, tmin1 = tmin1, tmin2 = tmin2,
    tname = tname, data = data,
    idname = idname, probs = probs,
    iters = iters, alp = alp, method = method,
    se = se, retEachIter = retEachIter,
    pl = pl, cores = cores, panel = TRUE, bootstrapiter = FALSE
  )


  ## Do some error handling for common cases...

  ## setup the data as build on this a bit...
  setupData(qp)

  ## do some checking that format of data ok
  panel.checks(qp)

  ## first calculate the actual estimate
  pqte <- compute.panel.qtet(qp)


  ## compute standard errors if so desired
  if (se) {
    qp$bootstrapiter <- TRUE

    ## bootstrap the standard errors
    ## the bootstrap method does this generically
    SEobj <- bootstrap(qp, pqte, compute.panel.qtet)

    ## set the results
    out <- QTE(
      qte = pqte$qte, qte.upper = SEobj$qte.upper,
      qte.lower = SEobj$qte.lower, ate = pqte$ate,
      ate.upper = SEobj$ate.upper, ate.lower = SEobj$ate.lower,
      qte.se = SEobj$qte.se, ate.se = SEobj$ate.se,
      c = SEobj$c,
      alp = alp,
      F.treated.t = pqte$F.treated.t,
      F.untreated.t = pqte$F.untreated.t,
      F.treated.t.cf = pqte$F.treated.t.cf,
      F.treated.tmin1 = pqte$F.treated.tmin1,
      F.treated.tmin2 = pqte$F.treated.tmin2,
      F.treated.change.tmin1 = pqte$F.treated.change.tmin1,
      F.untreated.change.t = pqte$F.untreated.change.t,
      F.untreated.change.tmin1 = pqte$F.untreated.change.tmin1,
      F.untreated.tmin1 = pqte$F.untreated.tmin1,
      F.untreated.tmin2 = pqte$F.untreated.tmin2,
      pscore.reg = pqte$pscore.reg,
      eachIterList = eachIter,
      probs = probs
    )
    return(out)
  } else {
    return(pqte)
  }
}


###### GENERAL HELPER FUNCTIONS#######

## return an SE object
## bootIters should contain ATT as first object in list
#' @title computeSE
#'
#' @description Computes standard errors from bootstrap results.  This function
#'  is called by several functions in the qte package
#'
#' @param bootIters List of bootstrap iterations
#' @inheritParams panel.qtet
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
