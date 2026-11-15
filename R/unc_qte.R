# =============================================================================
# Title: Cross-sectional QTE/QTT under unconfoundedness (unc_qte)
# Description: Estimates the Quantile Treatment Effect (QTE) or Quantile
#   Treatment Effect on the Treated (QTT) under unconfoundedness via three
#   methods: IPW propensity-score reweighting (Firpo 2007), outcome regression
#   via quantile regression inversion (OR), and a doubly-robust AIPW
#   combination of IPW and OR. Standard errors via empirical bootstrap.
#   (The deprecated wrappers ci.qte() and ci.qtet() were removed in v2.1.0.)
# Author: Brant Callaway
# Last update: 2026-11-15
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
  cond_dist     <- if (is.null(qp$cond_dist))  "qr"  else qp$cond_dist
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
  # OR branch: average conditional distributions over X (via predicted
  # conditional quantiles) to approximate the marginal distribution
  # ---------------------------------------------------------------------------
  } else if (est_method == "or") {
    # Conditional outcome model on untreated (counterfactual arm): quantile or
    # distribution regression (cond_dist); see cond_fit().
    M0 <- cond_fit(or_formla, untreated, cond_dist, tau_grid, weighted = FALSE)

    if (target == "qte") {
      # Predict at all n units; integrate out X to get marginal distribution
      pred0 <- cond_quantile_matrix(M0, data, tau_grid)
      w_rep <- rep(w, each = K)
      q0    <- weighted_quantile(probs, c(t(pred0)), w_rep, norm = TRUE)

      M1    <- cond_fit(or_formla, treated, cond_dist, tau_grid, weighted = FALSE)
      pred1 <- cond_quantile_matrix(M1, data, tau_grid)
      q1    <- weighted_quantile(probs, c(t(pred1)), w_rep, norm = TRUE)

      ate <- weighted_mean(c(t(pred1)), w_rep) -
             weighted_mean(c(t(pred0)), w_rep)
    } else {
      # QTT: counterfactual predicted at treated units only
      pred0  <- cond_quantile_matrix(M0, treated, tau_grid)
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

    # Outcome model for counterfactual arm (quantile or distribution
    # regression, cond_dist; see cond_fit()).
    M0        <- cond_fit(or_formla, untreated, cond_dist, tau_grid, weighted = FALSE)
    pred0_all <- cond_quantile_matrix(M0, data, tau_grid)

    # Evaluate CDFs at each observed Y value; approx() inverts
    y_grid <- sort(unique(y))

    if (target == "qte") {
      # AIPW CDF for Y(0): F_c(y) from OR + IPW augmentation for bias correction
      F0 <- vapply(y_grid, function(yy) {
        F_cond <- rowMeans(pred0_all <= yy)
        mean(F_cond + (1 - D) / (1 - pscore) * (as.numeric(y <= yy) - F_cond))
      }, numeric(1))
      F0 <- cummax(pmax(0, pmin(1, F0)))
      q0 <- approx(F0[!duplicated(F0)], y_grid[!duplicated(F0)], xout = probs, rule = 2)$y

      # AIPW CDF for Y(1)
      M1        <- cond_fit(or_formla, treated, cond_dist, tau_grid, weighted = FALSE)
      pred1_all <- cond_quantile_matrix(M1, data, tau_grid)
      F1 <- vapply(y_grid, function(yy) {
        F_cond <- rowMeans(pred1_all <= yy)
        mean(F_cond + D / pscore * (as.numeric(y <= yy) - F_cond))
      }, numeric(1))
      F1 <- cummax(pmax(0, pmin(1, F1)))
      q1 <- approx(F1[!duplicated(F1)], y_grid[!duplicated(F1)], xout = probs, rule = 2)$y

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
      q0    <- approx(F_qtt[!duplicated(F_qtt)], y_grid[!duplicated(F_qtt)], xout = probs, rule = 2)$y

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
#' @param cond_dist For \code{est_method = "or"} or \code{"aipw"}, how the
#'   conditional outcome distributions are estimated: \code{"qr"} (default)
#'   for quantile regression or \code{"dr"} for distribution regression (a
#'   logit of \eqn{1\{Y \le y\}} on the covariates at a grid of thresholds;
#'   Chernozhukov, Fernandez-Val and Melly 2013). Ignored for \code{"ipw"}.
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
                    cond_dist = c("qr", "dr"),
                    cl = 1) {
  method     <- match.arg(method)
  est_method <- match.arg(est_method)
  target     <- match.arg(target)
  cond_dist  <- match.arg(cond_dist)
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
  qp$cond_dist <- cond_dist

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
