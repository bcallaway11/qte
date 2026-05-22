# =============================================================================
# Title: Cross-sectional QTE under unconfoundedness (unc_qte)
# Description: Implements Firpo (2007) propensity-score reweighting estimator
#   for the Quantile Treatment Effect (QTE) under unconfoundedness (selection
#   on observables). Treated and untreated quantile functions are reweighted
#   by IPW weights derived from an estimated propensity score.
# Author: Brant Callaway
# Last update: 2026-05-22
# Date created: 2026-05-18
# =============================================================================

utils::globalVariables(c("yname", "treat", "treated", "x", "xformla", "data",
                         "wname", "probs", "method", "untreated", "eachIter"))

# --- Internal compute function -----------------------------------------------

#' @title compute.unc_qte
#'
#' @description Computes the QTE under unconfoundedness using Firpo (2007)
#'   propensity score reweighting. Called internally by \code{unc_qte} and
#'   during each bootstrap iteration.
#'
#' @param qp QTEparams object
#'
#' @import Hmisc
#'
#' @keywords internal
#'
#' @return QTE object
compute.unc_qte <- function(qp) {
  setupData(qp)
  bootstrapiter <- qp$bootstrapiter

  ate <- weighted_mean(treated[, yname], treated[, wname]) -
    weighted_mean(untreated[, yname], untreated[, wname])

  treated.firpo.quantiles <- weighted_quantile(
    probs, treated[, yname], treated[, wname]
  )
  untreated.firpo.quantiles <- weighted_quantile(
    probs, untreated[, yname], untreated[, wname]
  )
  qte <- treated.firpo.quantiles - untreated.firpo.quantiles

  n <- nrow(data)
  pscore.reg <- NULL

  if (!is.null(x)) {
    D <- data[, treat]
    y <- data[, yname]
    w <- data[, wname]

    this.formla <- y ~ x
    formula.tools::lhs(this.formla) <- as.name(treat)
    formula.tools::rhs(this.formla) <- formula.tools::rhs(xformla)
    pscore.reg <- glm(this.formla, data = data,
                      family = binomial(link = method))
    pscore <- fitted(pscore.reg)

    treated.weights   <- w * D / pscore
    untreated.weights <- w * (1 - D) / (1 - pscore)

    treated.firpo.quantiles   <- weighted_quantile(probs, y,
                                                   treated.weights,
                                                   norm = TRUE)
    untreated.firpo.quantiles <- weighted_quantile(probs, y,
                                                   untreated.weights,
                                                   norm = TRUE)
    qte <- treated.firpo.quantiles - untreated.firpo.quantiles
    ate <- weighted_mean(y, treated.weights) -
      weighted_mean(y, untreated.weights)
  }

  if (is.null(untreated.firpo.quantiles)) {
    F.treated.t.cf <- ecdf(untreated[, yname])
    F.treated.t    <- ecdf(treated[, yname])
  } else {
    F.treated.t.cf <- ecdf(untreated.firpo.quantiles)
    F.treated.t    <- ecdf(treated.firpo.quantiles)
  }

  if (bootstrapiter) {
    out <- QTE(qte = qte, ate = ate, probs = probs)
  } else {
    out <- QTE(F.treated.t = F.treated.t, F.treated.t.cf = F.treated.t.cf,
               qte = qte, pscore.reg = pscore.reg, ate = ate, probs = probs)
  }
  return(out)
}


# --- User-facing function -----------------------------------------------------

#' @title unc_qte
#'
#' @description Estimates the Quantile Treatment Effect (QTE) under
#'   unconfoundedness (selection on observables) using the method of
#'   Firpo (2007). Propensity score reweighting constructs the counterfactual
#'   quantile function. Standard errors and confidence bands are computed via
#'   the empirical bootstrap.
#'
#' @param yname character; name of the outcome variable in \code{data}.
#' @param dname character; name of the binary treatment indicator in
#'   \code{data} (1 = treated, 0 = untreated).
#' @param data data.frame containing the analysis data.
#' @param xformla one-sided formula for covariates used to estimate the
#'   propensity score, e.g. \code{~ age + education}. Default \code{~1}
#'   uses no covariates (simple quantile differences).
#' @param weightsname character; name of a column in \code{data} containing
#'   sampling weights. Default \code{NULL} applies equal weights.
#' @param probs numeric vector of quantile levels at which to evaluate the
#'   QTE. Default \code{seq(0.05, 0.95, 0.05)}.
#' @param alp significance level for confidence intervals. Default 0.05.
#' @param biters number of bootstrap iterations. Default 100.
#' @param cband logical; if \code{TRUE} (default), compute a uniform
#'   confidence band in addition to pointwise intervals.
#' @param boot_type bootstrap type. Currently only \code{"empirical"} is
#'   supported.
#' @param method propensity score model. One of \code{"logit"} (default) or
#'   \code{"probit"}.
#' @param cl number of cores for parallel bootstrap. Default 1 (sequential).
#'
#' @references
#' Firpo, Sergio. ``Efficient Semiparametric Estimation of Quantile Treatment
#'   Effects.'' Econometrica 75.1, pp. 259-276, 2007.
#'
#' @examples
#' data(lalonde)
#'
#' ## No covariates
#' q1 <- unc_qte(yname = "re78", dname = "treat", data = lalonde.psid,
#'               biters = 50, probs = seq(0.05, 0.95, 0.05))
#' summary(q1)
#'
#' ## With covariates
#' q2 <- unc_qte(yname = "re78", dname = "treat", data = lalonde.psid,
#'               xformla = ~ age + I(age^2) + education + black +
#'                           hispanic + married + nodegree,
#'               biters = 50, probs = seq(0.05, 0.95, 0.05))
#' summary(q2)
#'
#' @return QTE object
#' @export
unc_qte <- function(yname, dname, data,
                    xformla = ~1,
                    weightsname = NULL,
                    probs = seq(0.05, 0.95, 0.05),
                    alp = 0.05,
                    biters = 100,
                    cband = TRUE,
                    boot_type = "empirical",
                    method = c("logit", "probit"),
                    cl = 1) {
  method    <- match.arg(method)
  boot_type <- match.arg(boot_type)

  formla <- as.formula(paste(yname, "~", dname))
  w      <- if (!is.null(weightsname)) data[[weightsname]] else NULL
  pl     <- cl > 1
  cores  <- cl

  qp <- QTEparams(formla, xformla, t = NULL, tmin1 = NULL, tmin2 = NULL,
                  tname = NULL, data = data, w = w, idname = NULL,
                  probs = probs, iters = biters, alp = alp, method = method,
                  se = TRUE, retEachIter = FALSE, bootstrapiter = FALSE,
                  pl = pl, cores = cores)

  firpo.qte <- compute.unc_qte(qp)

  qp$bootstrapiter <- TRUE
  SEobj <- bootstrap(qp, firpo.qte, compute.unc_qte)

  QTE(qte       = firpo.qte$qte,
      qte.upper = SEobj$qte.upper,  qte.lower = SEobj$qte.lower,
      ate       = firpo.qte$ate,
      ate.upper = SEobj$ate.upper,  ate.lower = SEobj$ate.lower,
      qte.se    = SEobj$qte.se,     ate.se    = SEobj$ate.se,
      c = SEobj$c, alp = alp,
      pscore.reg    = firpo.qte$pscore.reg,
      F.treated.t   = firpo.qte$F.treated.t,
      F.treated.t.cf = firpo.qte$F.treated.t.cf,
      eachIterList  = eachIter,
      probs = probs)
}


# --- Deprecated wrapper -------------------------------------------------------

#' @title ci.qte
#'
#' @description Deprecated. Use \code{\link{unc_qte}} instead.
#'
#' @param formla formula of the form \code{y ~ treatment}.
#' @param xformla one-sided covariate formula. Default \code{NULL}.
#' @param x ignored (deprecated argument).
#' @param data data.frame.
#' @param w sampling weight vector. Cannot be forwarded automatically;
#'   use \code{weightsname} in \code{unc_qte} instead.
#' @param probs quantile grid.
#' @param se ignored; bootstrap is always run in \code{unc_qte}.
#' @param iters number of bootstrap iterations.
#' @param alp significance level.
#' @param method propensity score model.
#' @param retEachIter ignored.
#' @param printIter ignored.
#' @param pl logical; whether to run bootstrap in parallel.
#' @param cores number of cores.
#' @param biters alias for \code{iters}; takes precedence if supplied.
#' @param cl alias for cores; takes precedence if supplied.
#'
#' @keywords internal
#' @export
ci.qte <- function(formla, xformla = NULL, x = NULL, data, w = NULL,
                   probs = seq(0.05, 0.95, 0.05), se = TRUE,
                   iters = 100, alp = 0.05, method = "logit",
                   retEachIter = FALSE, printIter = FALSE,
                   pl = FALSE, cores = 2, biters = NULL, cl = NULL) {
  .Deprecated("unc_qte")
  if (!is.null(w)) {
    warning("ci.qte: 'w' cannot be forwarded to unc_qte; ",
            "use 'weightsname' in unc_qte instead.", call. = FALSE)
  }
  if (!is.null(biters)) iters <- biters
  if (!is.null(cl)) { pl <- cl > 1; cores <- cl }

  form  <- as.formula(formla)
  yname <- as.character(form[[2]])
  dname <- as.character(form[[3]])

  unc_qte(yname = yname, dname = dname, data = data,
          xformla = if (is.null(xformla)) ~1 else xformla,
          probs = probs, alp = alp, biters = iters,
          method = method, cl = if (pl) cores else 1L)
}
