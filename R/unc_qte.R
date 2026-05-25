# =============================================================================
# Title: Cross-sectional QTE/QTT under unconfoundedness (unc_qte / ci.qte / ci.qtet)
# Description: Estimates the Quantile Treatment Effect (QTE) or Quantile
#   Treatment Effect on the Treated (QTT) under unconfoundedness via three
#   methods: IPW propensity-score reweighting (Firpo 2007), outcome regression
#   via quantile regression inversion (OR), and a doubly-robust AIPW
#   combination of IPW and OR. Standard errors via empirical bootstrap.
#   Also contains deprecated wrappers ci.qte and ci.qtet.
# Author: Brant Callaway
# Last update: 2026-05-24
# Date created: 2026-05-18
# =============================================================================

utils::globalVariables(c("yname", "treat", "treated", "x", "xformla", "data",
                         "wname", "probs", "method", "untreated", "eachIter"))

# --- Internal compute function -----------------------------------------------

#' @title compute.unc_qte
#'
#' @description Computes the QTE or QTT under unconfoundedness. Dispatches to
#'   IPW, OR (outcome regression), or AIPW based on \code{qp$est_method} and
#'   \code{qp$target}. Called internally by \code{unc_qte} and during each
#'   bootstrap iteration.
#'
#' @param qp QTEparams object
#'
#' @keywords internal
#'
#' @return QTE object
compute.unc_qte <- function(qp) {
  setupData(qp)
  bootstrapiter <- qp$bootstrapiter
  est_method    <- if (is.null(qp$est_method)) "ipw" else qp$est_method
  target        <- if (is.null(qp$target))     "qte" else qp$target

  D    <- data[, treat]
  y    <- data[, yname]
  w    <- data[, wname]

  pscore.reg <- NULL

  # Formulas built from the injected xformla (RHS extracted via [[2]])
  xrhs      <- deparse(xformla[[2]])
  ps_formla <- as.formula(paste(treat, "~", xrhs))
  or_formla <- as.formula(paste(yname, "~", xrhs))

  # Dense tau grid for the OR / AIPW outcome models
  tau_grid <- seq(0.01, 0.99, 0.01)
  K        <- length(tau_grid)

  # ---------------------------------------------------------------------------
  # IPW branch
  # ---------------------------------------------------------------------------
  if (est_method == "ipw") {
    if (is.null(x)) {
      # No covariates: naive weighted quantile differences
      q1  <- weighted_quantile(probs, treated[, yname], treated[, wname])
      q0  <- weighted_quantile(probs, untreated[, yname], untreated[, wname])
      ate <- weighted_mean(treated[, yname], treated[, wname]) -
             weighted_mean(untreated[, yname], untreated[, wname])
    } else {
      pscore.reg <- glm(ps_formla, data = data, family = binomial(link = method))
      pscore     <- pmax(1e-6, pmin(1 - 1e-6, fitted(pscore.reg)))

      if (target == "qte") {
        t_wts <- w * D       / pscore
        u_wts <- w * (1 - D) / (1 - pscore)
        q1  <- weighted_quantile(probs, y, t_wts, norm = TRUE)
        q0  <- weighted_quantile(probs, y, u_wts, norm = TRUE)
        ate <- weighted_mean(y, t_wts) - weighted_mean(y, u_wts)
      } else {
        # QTT: Q_Y(1)|D=1 from treated data; Q_Y(0)|D=1 by reweighting untreated
        # by p(x)/(1-p(x)), tilting toward the treated covariate distribution
        q1     <- weighted_quantile(probs, treated[, yname], treated[, wname])
        att1   <- weighted_mean(treated[, yname], treated[, wname])
        ps_u   <- pscore[D == 0]
        w_u    <- w[D == 0]
        u_qtt  <- w_u * ps_u / (1 - ps_u)
        q0     <- weighted_quantile(probs, untreated[, yname], u_qtt, norm = TRUE)
        att0   <- weighted_mean(untreated[, yname], u_qtt)
        ate    <- att1 - att0
      }
    }

  # ---------------------------------------------------------------------------
  # OR branch: invert quantile regression to approximate marginal distribution
  # ---------------------------------------------------------------------------
  } else if (est_method == "or") {
    # Fit conditional quantile model on untreated (counterfactual arm).
    # Suppress "Solution may be nonunique" from quantreg LP solver — benign.
    M0 <- suppressWarnings(
      quantreg::rq(or_formla, tau = tau_grid, data = untreated)
    )

    if (target == "qte") {
      # Predict at all n units; integrate out X to get marginal distribution
      pred0 <- as.matrix(predict(M0, newdata = data))
      # Rearrange each unit's conditional quantile function to be monotone
      pred0 <- t(apply(pred0, 1, sort))
      w_rep <- rep(w, each = K)
      q0    <- weighted_quantile(probs, c(t(pred0)), w_rep, norm = TRUE)

      M1    <- suppressWarnings(
        quantreg::rq(or_formla, tau = tau_grid, data = treated)
      )
      pred1 <- as.matrix(predict(M1, newdata = data))
      pred1 <- t(apply(pred1, 1, sort))
      q1    <- weighted_quantile(probs, c(t(pred1)), w_rep, norm = TRUE)

      ate <- weighted_mean(c(t(pred1)), w_rep) -
             weighted_mean(c(t(pred0)), w_rep)
    } else {
      # QTT: counterfactual predicted at treated units only
      pred0  <- as.matrix(predict(M0, newdata = treated))
      pred0  <- t(apply(pred0, 1, sort))
      wt_rep <- rep(treated[, wname], each = K)
      q0     <- weighted_quantile(probs, c(t(pred0)), wt_rep, norm = TRUE)

      q1   <- weighted_quantile(probs, treated[, yname], treated[, wname])
      att1 <- weighted_mean(treated[, yname], treated[, wname])
      att0 <- weighted_mean(c(t(pred0)), wt_rep)
      ate  <- att1 - att0
    }

  # ---------------------------------------------------------------------------
  # AIPW branch: doubly-robust CDF estimator, then invert
  # ---------------------------------------------------------------------------
  } else if (est_method == "aipw") {
    pscore.reg <- glm(ps_formla, data = data, family = binomial(link = method))
    pscore     <- pmax(1e-6, pmin(1 - 1e-6, fitted(pscore.reg)))

    # Outcome model for counterfactual arm.
    # Suppress "Solution may be nonunique" from quantreg LP solver — benign.
    M0        <- suppressWarnings(
      quantreg::rq(or_formla, tau = tau_grid, data = untreated)
    )
    pred0_all <- as.matrix(predict(M0, newdata = data))
    pred0_all <- t(apply(pred0_all, 1, sort))

    # Evaluate CDFs at each observed Y value; approx() inverts
    y_grid <- sort(unique(y))

    if (target == "qte") {
      # AIPW CDF for Y(0): F_c(y) from OR + IPW augmentation for bias correction
      F0 <- vapply(y_grid, function(yy) {
        F_cond <- rowMeans(pred0_all <= yy)
        mean(F_cond + (1 - D) / (1 - pscore) * (as.numeric(y <= yy) - F_cond))
      }, numeric(1))
      F0 <- cummax(pmax(0, pmin(1, F0)))
      q0 <- approx(F0, y_grid, xout = probs, rule = 2, ties = min)$y

      # AIPW CDF for Y(1)
      M1        <- suppressWarnings(
        quantreg::rq(or_formla, tau = tau_grid, data = treated)
      )
      pred1_all <- as.matrix(predict(M1, newdata = data))
      pred1_all <- t(apply(pred1_all, 1, sort))
      F1 <- vapply(y_grid, function(yy) {
        F_cond <- rowMeans(pred1_all <= yy)
        mean(F_cond + D / pscore * (as.numeric(y <= yy) - F_cond))
      }, numeric(1))
      F1 <- cummax(pmax(0, pmin(1, F1)))
      q1 <- approx(F1, y_grid, xout = probs, rule = 2, ties = min)$y

      t_wts <- w * D       / pscore
      u_wts <- w * (1 - D) / (1 - pscore)
      ate   <- weighted_mean(y, t_wts) - weighted_mean(y, u_wts)

    } else {
      # AIPW CDF for Y(0)|D=1: weighted by p(x)/pbar to target treated distribution
      pbar   <- mean(D)
      F_qtt  <- vapply(y_grid, function(yy) {
        F_cond    <- rowMeans(pred0_all <= yy)
        aipw_term <- D * F_cond +
                     (1 - D) * pscore / (1 - pscore) * (as.numeric(y <= yy) - F_cond)
        mean(aipw_term) / pbar
      }, numeric(1))
      F_qtt <- cummax(pmax(0, pmin(1, F_qtt)))
      q0    <- approx(F_qtt, y_grid, xout = probs, rule = 2, ties = min)$y

      q1   <- weighted_quantile(probs, treated[, yname], treated[, wname])
      att1 <- weighted_mean(treated[, yname], treated[, wname])
      ps_u <- pscore[D == 0]
      w_u  <- w[D == 0]
      att0 <- weighted_mean(untreated[, yname], w_u * ps_u / (1 - ps_u))
      ate  <- att1 - att0
    }
  }

  qte <- q1 - q0
  F.treated.t    <- ecdf(q1)
  F.treated.t.cf <- ecdf(q0)

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
#' @description Estimates the Quantile Treatment Effect (QTE) or Quantile
#'   Treatment Effect on the Treated (QTT) under unconfoundedness, also known
#'   as selection on observables. The key identifying assumption is
#'   \eqn{(Y(0), Y(1)) \perp D \mid X}, i.e., potential outcomes are
#'   independent of treatment conditional on covariates \eqn{X}.
#'
#'   Three estimation methods are available via \code{est_method}:
#'
#'   \describe{
#'     \item{\code{"ipw"}}{Propensity-score reweighting (Firpo 2007). The
#'       propensity score \eqn{p(X) = P(D=1|X)} is estimated by logit or
#'       probit and used to reweight the sample so that the covariate
#'       distributions of the treated and untreated groups match. For the QTE,
#'       both groups are reweighted toward the population distribution; for the
#'       QTT, only the untreated group is reweighted toward the treated
#'       covariate distribution.}
#'     \item{\code{"or"}}{Outcome regression via quantile regression inversion.
#'       \code{\link[quantreg]{rq}} is fit on a dense internal grid of
#'       \eqn{\tau} values, yielding a conditional quantile function
#'       \eqn{Q_{Y|X,D=j}(\tau)} for each arm \eqn{j}. The marginal
#'       distribution is recovered by averaging conditional quantile functions
#'       over the empirical covariate distribution (Melly 2006; Chernozhukov,
#'       Fernandez-Val, and Melly 2013).}
#'     \item{\code{"aipw"}}{Doubly-robust augmented IPW. Combines the
#'       propensity score model with the conditional quantile outcome model.
#'       The CDF estimator is consistent if either the propensity score or the
#'       outcome model is correctly specified (semiparametric efficiency when
#'       both are correct).}
#'   }
#'
#'   When \code{xformla = ~1} (the default), all three methods reduce to
#'   simple quantile differences between the treated and untreated groups, which
#'   is consistent under unconditional unconfoundedness. Covariates are required
#'   for covariate-adjusted estimation.
#'
#'   Standard errors and uniform confidence bands are computed via the
#'   empirical bootstrap. Setting \code{biters = 0} is planned for a future
#'   version to skip inference and return point estimates only.
#'
#' @param yname character; name of the outcome variable in \code{data}.
#' @param dname character; name of the binary treatment indicator in
#'   \code{data} (1 = treated, 0 = untreated).
#' @param data data.frame containing the analysis data.
#' @param xformla one-sided formula for covariates used in the propensity
#'   score and/or outcome model, e.g. \code{~ age + I(age^2) + education}.
#'   Default \code{~1} uses no covariates (simple quantile differences).
#' @param weightsname character; name of a column in \code{data} containing
#'   sampling weights. Default \code{NULL} applies equal weights.
#' @param probs numeric vector of quantile levels at which to evaluate the
#'   QTE or QTT. Default \code{seq(0.05, 0.95, 0.05)}.
#' @param alp significance level for confidence intervals. Default 0.05.
#' @param biters number of bootstrap iterations for standard errors and
#'   confidence bands. Default 100.
#' @param cband logical; if \code{TRUE} (default), report a uniform confidence
#'   band (simultaneous over all quantiles in \code{probs}) in addition to
#'   pointwise intervals.
#' @param boot_type bootstrap variant. Currently only \code{"empirical"} is
#'   supported.
#' @param method propensity score link function: \code{"logit"} (default) or
#'   \code{"probit"}. Used by \code{est_method = "ipw"} and \code{"aipw"}.
#' @param est_method estimation method; one of:
#'   \itemize{
#'     \item \code{"ipw"} (default) — inverse propensity weighting
#'     \item \code{"or"} — outcome regression via quantile regression inversion
#'     \item \code{"aipw"} — doubly-robust augmented IPW
#'   }
#' @param target target parameter; one of:
#'   \itemize{
#'     \item \code{"qte"} (default) — population QTE:
#'       \eqn{Q_{Y(1)}(\tau) - Q_{Y(0)}(\tau)}
#'     \item \code{"qtt"} — QTT (effect on the treated):
#'       \eqn{Q_{Y(1)|D=1}(\tau) - Q_{Y(0)|D=1}(\tau)}
#'   }
#' @param cl number of cores for parallel bootstrap. Default 1 (sequential).
#'
#' @references
#' Firpo, Sergio. ``Efficient Semiparametric Estimation of Quantile Treatment
#'   Effects.'' \emph{Econometrica} 75(1), pp. 259--276, 2007.
#'
#' Melly, Blaise. ``Estimation of Counterfactual Distributions Using Quantile
#'   Regression.'' Working paper, University of St. Gallen, 2006.
#'
#' Chernozhukov, Victor, Ivan Fernandez-Val, and Blaise Melly.
#'   ``Inference on Counterfactual Distributions.''
#'   \emph{Econometrica} 81(6), pp. 2205--2268, 2013.
#'
#' @return An object of class \code{QTE} containing:
#'   \describe{
#'     \item{\code{qte}}{numeric vector of estimated QTE (or QTT) at each
#'       element of \code{probs}.}
#'     \item{\code{ate}}{estimated ATE (or ATT when \code{target = "qtt"}).}
#'     \item{\code{qte.se}}{bootstrap standard errors for \code{qte}.}
#'     \item{\code{qte.lower}, \code{qte.upper}}{confidence interval bounds
#'       for \code{qte}; uniform when \code{cband = TRUE}, pointwise otherwise.}
#'     \item{\code{ate.se}, \code{ate.lower}, \code{ate.upper}}{SE and CI
#'       for \code{ate}.}
#'     \item{\code{probs}}{the \code{probs} vector passed by the user.}
#'     \item{\code{pscore.reg}}{fitted propensity score \code{glm} object, or
#'       \code{NULL} when \code{est_method = "or"} or \code{xformla = ~1}.}
#'   }
#'   Use \code{summary()} to print a formatted table and \code{plot()} to
#'   display the QTE curve with confidence bands.
#'
#' @examples
#' data(lalonde)
#'
#' ## IPW, no covariates (simple quantile differences)
#' \donttest{
#' q1 <- unc_qte(yname = "re78", dname = "treat", data = lalonde.psid,
#'               biters = 20, probs = seq(0.05, 0.95, 0.05))
#' summary(q1)
#' plot(q1)
#' }
#'
#' ## OR with covariates, QTE
#' \donttest{
#' xf <- ~ age + I(age^2) + education + black + hispanic + married + nodegree
#' q2 <- unc_qte(yname = "re78", dname = "treat", data = lalonde.psid,
#'               xformla = xf, est_method = "or",
#'               biters = 20, probs = seq(0.05, 0.95, 0.05))
#' summary(q2)
#' }
#'
#' ## AIPW with covariates, QTT
#' \donttest{
#' q3 <- unc_qte(yname = "re78", dname = "treat", data = lalonde.psid,
#'               xformla = xf, est_method = "aipw", target = "qtt",
#'               biters = 20, probs = seq(0.05, 0.95, 0.05))
#' summary(q3)
#' }
#'
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
                    est_method = c("ipw", "or", "aipw"),
                    target = c("qte", "qtt"),
                    cl = 1) {
  method     <- match.arg(method)
  est_method <- match.arg(est_method)
  target     <- match.arg(target)
  boot_type  <- match.arg(boot_type)

  formla <- as.formula(paste(yname, "~", dname))
  w      <- if (!is.null(weightsname)) data[[weightsname]] else NULL
  pl     <- cl > 1
  cores  <- cl

  qp <- QTEparams(formla, xformla, t = NULL, tmin1 = NULL, tmin2 = NULL,
                  tname = NULL, data = data, w = w, idname = NULL,
                  probs = probs, iters = biters, alp = alp, method = method,
                  se = TRUE, retEachIter = FALSE, bootstrapiter = FALSE,
                  pl = pl, cores = cores,
                  est_method = est_method, target = target)

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
      probs = probs,
      type = if (target == "qtt") "On the Treated" else "Overall")
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
#' @return A \code{QTE} object; same structure as \code{\link{unc_qte}}.
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


# --- ci.qtet (deprecated) ----------------------------------------------------

#' @title ci.qtet
#'
#' @description \strong{Deprecated.} Use \code{\link{unc_qte}}\code{(target = "qtt")}
#'   for a cross-sectional QTT under unconfoundedness, or \code{\link{lou_qtt}}
#'   for staggered treatment adoption with optional lagged-outcome conditioning.
#'
#'   \code{ci.qtet} estimates the Quantile Treatment Effect on the Treated
#'   (QTET) under a Conditional Independence Assumption (sometimes called
#'   Selection on Observables) following Firpo (2007). It uses propensity score
#'   re-weighting to estimate the counterfactual distribution. Standard errors
#'   are computed via the bootstrap.
#'
#' @inheritParams panel.qtet
#' @inheritParams ci.qte
#' @param method Method to compute propensity score.  Default is logit; other
#'  option is probit.
#' @param indsample Binary variable for whether to treat the samples as
#'  independent or dependent.  This affects bootstrap standard errors.  In
#'  the job training example, the samples are independent because they
#'  are two samples collected independently and then merged.  If the data is
#'  from the same source, usually should set this option to be FALSE.
#' @param printIter For debugging only; should leave at default FALSE unless
#'  you want to see a lot of output
#' @param biters Number of bootstrap iterations; alias for \code{iters}
#'   matching the \code{did}/\code{ptetools} naming convention. If both are
#'   supplied, \code{biters} takes precedence.
#' @param cl Number of cores for parallel bootstrap; alias for
#'   \code{pl}/\code{cores}. \code{cl = 1} (default) runs sequentially;
#'   \code{cl > 1} enables parallelism.
#'
#' @references
#' Firpo, Sergio. ``Efficient Semiparametric Estimation of Quantile Treatment
#'  Effects.'' Econometrica 75.1, pp. 259-276, 2007.
#'
#' @examples
#' # See ?unc_qte for the modern replacement.
#' \dontrun{
#' data(lalonde)
#'
#' ## Estimate the QTET of participating in the job training program;
#' ## This is the no covariate case.  Note: Because individuals that participate
#' ## in the job training program are likely to be much different than
#' ## individuals that do not (e.g. less experience and less education), this
#' ## method is likely to perform poorly at estimating the true QTET
#' q1 <- ci.qtet(re78 ~ treat,
#'   x = NULL, data = lalonde.psid, se = FALSE,
#'   probs = seq(0.05, 0.95, 0.05)
#' )
#' summary(q1)
#'
#' ## This estimation controls for all the available background characteristics.
#' q2 <- ci.qtet(re78 ~ treat,
#'   xformla = ~ age + I(age^2) + education + black + hispanic + married + nodegree,
#'   data = lalonde.psid, se = FALSE, probs = seq(0.05, 0.95, 0.05)
#' )
#' summary(q2)
#' }
#'
#' @return A \code{QTE} object; same structure as \code{\link{unc_qte}}.
#' @keywords internal
#' @export
ci.qtet <- function(formla, xformla = NULL, w = NULL, data,
                    probs = seq(0.05, 0.95, 0.05), se = TRUE,
                    iters = 100, alp = 0.05, method = "logit",
                    retEachIter = FALSE, indsample = TRUE,
                    printIter = FALSE, pl = FALSE, cores = 2,
                    biters = NULL, cl = NULL) {
  .Deprecated(
    msg = paste0(
      "'ci.qtet' is deprecated.\n",
      "For a cross-sectional QTT under unconfoundedness use:\n",
      "  unc_qte(..., target = 'qtt')\n",
      "For staggered treatment adoption with optional lagged-outcome conditioning use:\n",
      "  lou_qtt(...)"
    )
  )
  if (!is.null(w)) {
    warning("ci.qtet: 'w' cannot be forwarded to unc_qte; ",
            "use 'weightsname' in unc_qte instead.", call. = FALSE)
  }
  if (!is.null(biters)) iters <- biters
  if (!is.null(cl)) { pl <- cl > 1; cores <- cl }

  form  <- as.formula(formla)
  yname <- as.character(form[[2]])
  dname <- as.character(form[[3]])

  unc_qte(yname = yname, dname = dname, data = data,
          xformla = if (is.null(xformla)) ~1 else xformla,
          target = "qtt",
          probs = probs, alp = alp, biters = iters,
          method = method, cl = if (pl) cores else 1L)
}
