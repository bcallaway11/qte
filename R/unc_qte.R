# =============================================================================
# Title: Cross-sectional QTE/QTT under unconfoundedness (unc_qte)
# Description: Estimates the Quantile Treatment Effect (QTE) or Quantile
#   Treatment Effect on the Treated (QTT) under unconfoundedness via three
#   methods: IPW propensity-score reweighting (Firpo 2007), outcome regression
#   via quantile regression inversion (OR), and a doubly-robust AIPW
#   combination of IPW and OR. Standard errors via empirical bootstrap.
#   Self-contained: unc_qte() is the cross-sectional estimator and does not go
#   through ptetools (the staggered analogue is lou_qtt()). As of v2.1.0 it no
#   longer uses the legacy QTEparams/setupData/bootstrap machinery.
# Author: Brant Callaway
# Last update: 2026-11-19
# Date created: 2026-05-18
# =============================================================================

# --- Internal: data preparation ----------------------------------------------

# unc_qte_data
# Builds the working data frame used by the estimator and every bootstrap
# draw: the outcome, treatment and sampling weights under fixed internal
# column names (.y, .d, .w, and .cluster with clustervars), plus the
# covariates named in xformula under their own names. Rows with a missing
# value in any of these are dropped (with a warning), so the outcome,
# treatment, weights, clusters and covariates stay aligned.
#
# Arguments: as in unc_qte(). Returns a data.frame.
unc_qte_data <- function(data, yname, dname, xformula, weightsname,
                         clustervars = NULL) {
  data <- as.data.frame(data)
  if (!is.null(clustervars) &&
        (!is.character(clustervars) || length(clustervars) != 1L)) {
    stop("`clustervars` must be the name of one column of `data`.")
  }
  needed <- c(yname, dname, weightsname, clustervars, all.vars(xformula))
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0) {
    stop("column(s) not found in `data`: ",
         paste(missing_cols, collapse = ", "), ".")
  }

  dat <- data[, unique(all.vars(xformula)), drop = FALSE]
  dat$.y <- data[[yname]]
  dat$.d <- as.numeric(data[[dname]])
  dat$.w <- if (is.null(weightsname)) 1 else data[[weightsname]]
  if (!is.null(clustervars)) dat$.cluster <- data[[clustervars]]

  keep <- stats::complete.cases(dat)
  if (!all(keep)) {
    warning("dropping ", sum(!keep), " row(s) with missing values in the ",
            "outcome, treatment, weights, clusters or covariates.")
    dat <- dat[keep, , drop = FALSE]
  }
  rownames(dat) <- NULL
  dat
}

# --- Internal: point estimates -----------------------------------------------

# unc_qte_est
# Computes the QTE or QTT (and the ATE or ATT) under unconfoundedness for one
# data set; called once on the full sample and once per bootstrap draw.
#
# Arguments:
#   dat         data.frame from unc_qte_data() (columns .y, .d, .w + covariates)
#   xformula     one-sided covariate formula
#   probs       quantile levels
#   method      propensity-score link ("logit" or "probit")
#   est_method  "ipw", "reg" or "aipw"
#   target      "qte" (population) or "qtt" (on the treated)
#   cond_dist   "qr" or "dr": how conditional outcome distributions are
#               estimated in the OR/AIPW branches (see cond_fit())
#   full        if TRUE also return pscore.reg and the ECDFs of the estimated
#               quantiles; bootstrap draws only need qte and ate
#
# Returns a list with qte, ate and, if full, pscore.reg, F.treated.t,
# F.treated.t.cf.
unc_qte_est <- function(dat, xformula, probs, method, est_method, target,
                        cond_dist, full = TRUE) {
  D <- dat$.d
  y <- dat$.y
  w <- dat$.w
  treated   <- dat[D == 1, , drop = FALSE]
  untreated <- dat[D == 0, , drop = FALSE]

  pscore.reg <- NULL

  # Formulas from the RHS of xformula (deparse of xformula[[2]]; "1" for ~1).
  # Note that ~1 still goes through the model-based branches below (with an
  # intercept-only model), which is what earlier versions did too.
  xrhs      <- paste(deparse(xformula[[2]]), collapse = " ")
  ps_formla <- stats::as.formula(paste(".d ~", xrhs))
  or_formla <- stats::as.formula(paste(".y ~", xrhs))

  # Dense tau grid for the OR / AIPW outcome models
  tau_grid <- seq(0.01, 0.99, 0.01)
  K        <- length(tau_grid)

  # Propensity score p(x) = P(D = 1 | X = x), fit by weighted (quasi)binomial
  # glm so sampling weights enter the score; trimmed away from 0 and 1.
  # weights = .w is evaluated in dat (the weights column), so a user
  # covariate that happens to be called `w` can't be picked up instead.
  fit_pscore <- function() {
    stats::glm(ps_formla, data = dat,
               family = stats::quasibinomial(link = method), weights = .w)
  }

  # ---------------------------------------------------------------------------
  # IPW branch (Firpo 2007)
  # ---------------------------------------------------------------------------
  if (est_method == "ipw") {
    pscore.reg <- fit_pscore()
    pscore     <- pmax(1e-6, pmin(1 - 1e-6, stats::fitted(pscore.reg)))

    if (target == "qte") {
      # reweight each arm to the population covariate distribution
      t_wts <- w * D       / pscore
      u_wts <- w * (1 - D) / (1 - pscore)
      q1  <- weighted_quantile(probs, y, t_wts, norm = TRUE)
      q0  <- weighted_quantile(probs, y, u_wts, norm = TRUE)
      ate <- weighted_mean(y, t_wts) - weighted_mean(y, u_wts)
    } else {
      # QTT: Q_Y(1)|D=1 from treated data; Q_Y(0)|D=1 by reweighting untreated
      # by p(x)/(1-p(x)), tilting toward the treated covariate distribution
      q1    <- weighted_quantile(probs, treated$.y, treated$.w)
      att1  <- weighted_mean(treated$.y, treated$.w)
      ps_u  <- pscore[D == 0]
      u_qtt <- w[D == 0] * ps_u / (1 - ps_u)
      q0    <- weighted_quantile(probs, untreated$.y, u_qtt, norm = TRUE)
      att0  <- weighted_mean(untreated$.y, u_qtt)
      ate   <- att1 - att0
    }

  # ---------------------------------------------------------------------------
  # OR branch: average conditional distributions over X (via predicted
  # conditional quantiles) to get the marginal distribution (Melly 2006;
  # Chernozhukov, Fernandez-Val and Melly 2013)
  # ---------------------------------------------------------------------------
  } else if (est_method == "reg") {
    # conditional outcome model on the untreated (counterfactual) arm
    M0 <- cond_fit(or_formla, untreated, cond_dist, tau_grid)

    if (target == "qte") {
      # predict at all n units; integrating out X gives the marginal
      pred0 <- cond_quantile_matrix(M0, dat, tau_grid)
      w_rep <- rep(w, each = K)
      q0    <- weighted_quantile(probs, c(t(pred0)), w_rep, norm = TRUE)

      M1    <- cond_fit(or_formla, treated, cond_dist, tau_grid)
      pred1 <- cond_quantile_matrix(M1, dat, tau_grid)
      q1    <- weighted_quantile(probs, c(t(pred1)), w_rep, norm = TRUE)

      ate <- weighted_mean(c(t(pred1)), w_rep) -
             weighted_mean(c(t(pred0)), w_rep)
    } else {
      # QTT: counterfactual predicted at treated units only
      pred0  <- cond_quantile_matrix(M0, treated, tau_grid)
      wt_rep <- rep(treated$.w, each = K)
      q0     <- weighted_quantile(probs, c(t(pred0)), wt_rep, norm = TRUE)

      q1   <- weighted_quantile(probs, treated$.y, treated$.w)
      att1 <- weighted_mean(treated$.y, treated$.w)
      att0 <- weighted_mean(c(t(pred0)), wt_rep)
      ate  <- att1 - att0
    }

  # ---------------------------------------------------------------------------
  # AIPW branch: doubly-robust CDF estimator, then invert
  # ---------------------------------------------------------------------------
  } else if (est_method == "aipw") {
    pscore.reg <- fit_pscore()
    pscore     <- pmax(1e-6, pmin(1 - 1e-6, stats::fitted(pscore.reg)))

    # outcome model for the counterfactual arm, predicted at all units
    M0        <- cond_fit(or_formla, untreated, cond_dist, tau_grid)
    pred0_all <- cond_quantile_matrix(M0, dat, tau_grid)

    # evaluate the CDFs at each observed outcome value; approx() inverts
    y_grid <- sort(unique(y))

    # invert a (monotonized) CDF evaluated on y_grid at the levels probs
    invert_cdf <- function(F_vals) {
      F_vals <- cummax(pmax(0, pmin(1, F_vals)))
      keep   <- !duplicated(F_vals)
      stats::approx(F_vals[keep], y_grid[keep], xout = probs, rule = 2)$y
    }

    if (target == "qte") {
      # AIPW CDF for Y(0): F_c(y) from OR + IPW augmentation for bias correction
      F0 <- vapply(y_grid, function(yy) {
        F_cond <- rowMeans(pred0_all <= yy)
        stats::weighted.mean(
          F_cond + (1 - D) / (1 - pscore) * (as.numeric(y <= yy) - F_cond), w
        )
      }, numeric(1))
      q0 <- invert_cdf(F0)

      # AIPW CDF for Y(1)
      M1        <- cond_fit(or_formla, treated, cond_dist, tau_grid)
      pred1_all <- cond_quantile_matrix(M1, dat, tau_grid)
      F1 <- vapply(y_grid, function(yy) {
        F_cond <- rowMeans(pred1_all <= yy)
        stats::weighted.mean(
          F_cond + D / pscore * (as.numeric(y <= yy) - F_cond), w
        )
      }, numeric(1))
      q1 <- invert_cdf(F1)

      # Doubly-robust ATE: AIPW means with the outcome models' conditional
      # means m_d(x) (row means of predicted quantiles), consistent with the
      # AIPW CDFs above (issue #11).
      m0  <- rowMeans(pred0_all)
      m1  <- rowMeans(pred1_all)
      ate <- weighted_mean(m1 + D / pscore * (y - m1), w) -
             weighted_mean(m0 + (1 - D) / (1 - pscore) * (y - m0), w)

    } else {
      # AIPW CDF for Y(0)|D=1: weighted by p(x)/pbar to target the treated
      pbar  <- stats::weighted.mean(D, w)
      F_qtt <- vapply(y_grid, function(yy) {
        F_cond    <- rowMeans(pred0_all <= yy)
        aipw_term <- D * F_cond +
          (1 - D) * pscore / (1 - pscore) * (as.numeric(y <= yy) - F_cond)
        stats::weighted.mean(aipw_term, w) / pbar
      }, numeric(1))
      q0 <- invert_cdf(F_qtt)

      q1   <- weighted_quantile(probs, treated$.y, treated$.w)
      att1 <- weighted_mean(treated$.y, treated$.w)
      # Doubly-robust ATT: E[Y(0)|D=1] = E[D*m0 + (1-D)*(p/(1-p))*(Y - m0)]/E[D]
      # (the mean of the AIPW counterfactual CDF above; issue #11)
      m0   <- rowMeans(pred0_all)
      att0 <- sum(w * (D * m0 + (1 - D) * pscore / (1 - pscore) * (y - m0))) /
              sum(w * D)
      ate  <- att1 - att0
    }
  }

  out <- list(qte = q1 - q0, ate = ate)
  if (full) {
    out$pscore.reg     <- trim_glm(pscore.reg)
    out$F.treated.t    <- stats::ecdf(q1)
    out$F.treated.t.cf <- stats::ecdf(q0)
  }
  out
}

# --- Internal: est_method "or" -> "reg" -------------------------------------

# "or" was renamed "reg" in 2.1.0 (the name did and DRDID use for outcome
# regression); "or" still works, with a deprecation warning.
resolve_est_method_or <- function(est_method) {
  if (identical(est_method, "or")) {
    warning('est_method = "or" is deprecated; use "reg" instead.',
            call. = FALSE)
    return("reg")
  }
  est_method
}

# --- Internal: slimmer glm for the returned object ----------------------------

# trim_glm
# Drops the parts of a fitted glm that only duplicate the training data
# ($data, $model, $y, $effects) and points its formula and terms at baseenv().
# The formula was built inside unc_qte_est(), so its environment held that
# function's whole frame (working data, treated/untreated subsets, ...), which
# saveRDS() wrote out with the fit: about 3.6MB of a 3.6MB QTE object for
# lalonde.psid, 0.4MB after trimming. coef(), vcov(), summary(), fitted(),
# residuals(), nobs() and predict(newdata = ) give identical results;
# anova() and update(), which refit on the data, no longer work.
#
# Args: fit a glm, or NULL. Returns the trimmed glm (or NULL).
trim_glm <- function(fit) {
  if (is.null(fit)) return(NULL)
  fit[c("data", "model", "y", "effects")] <- NULL
  environment(fit$formula) <- baseenv()
  attr(fit$terms, ".Environment") <- baseenv()
  fit
}

# --- Internal: bootstrap inference -------------------------------------------

# unc_qte_boot
# Empirical (nonparametric) bootstrap: resample rows of dat with replacement
# (or whole clusters, when dat has a .cluster column), re-estimate, and
# summarize the draws into standard errors and intervals.
#
# Arguments:
#   dat, est_args  working data and the remaining unc_qte_est() arguments
#                  (as a list)
#   est            full-sample list from unc_qte_est() (qte, ate)
#   biters         number of bootstrap draws
#   alp            significance level
#   cband          if TRUE, qte.lower/qte.upper are a uniform band over probs
#                  (qte +/- c * se); if FALSE, pointwise (qte +/- z * se)
#   cl             passed to pbapply::pblapply() (parallel workers)
#
# Returns a list with qte.se, qte.lower, qte.upper, ate.se, ate.lower,
# ate.upper and c (critical value for a uniform band over probs).
unc_qte_boot <- function(dat, est_args, est, biters, alp, cband, cl) {
  n <- nrow(dat)
  if (is.null(dat$.cluster)) {
    draw_rows <- function() sample(seq_len(n), n, replace = TRUE)
  } else {
    # cluster bootstrap: draw clusters with replacement and keep all rows of
    # each drawn cluster, so within-cluster dependence is preserved (as the
    # ptetools empirical bootstrap does with clustervars). With one row per
    # cluster, in row order, this gives exactly the same draws as resampling
    # rows, since sample(x, n, TRUE) is x[sample.int(length(x), n, TRUE)].
    clusters <- unique(dat$.cluster)
    rows <- split(seq_len(n), factor(dat$.cluster, levels = clusters))
    draw_rows <- function() {
      pos <- sample.int(length(clusters), length(clusters), replace = TRUE)
      unlist(rows[pos], use.names = FALSE)
    }
  }
  draws <- pbapply::pblapply(seq_len(biters), function(b) {
    idx <- draw_rows()
    do.call(unc_qte_est, c(list(dat = dat[idx, , drop = FALSE], full = FALSE),
                           est_args))
  }, cl = cl)

  qte_mat <- do.call(rbind, lapply(draws, function(r) as.numeric(r$qte)))
  ate_b   <- vapply(draws, function(r) as.numeric(r$ate), numeric(1))
  z       <- stats::qnorm(1 - alp / 2)

  # ATE/ATT: bootstrap SE, normal-approximation interval
  ate.se <- stats::sd(ate_b)

  # QTE: bootstrap SE per quantile, normal-approximation pointwise interval
  qte.se <- apply(qte_mat, 2, stats::sd)

  # Critical value for a uniform band: the (1 - alp) quantile of the sup over
  # probs of |boot - est| / sigma, with sigma a robust (IQR-based) scale; fall
  # back to the bootstrap SE when the IQR is 0 (e.g. QTEs that are exactly 0)
  sigmahalf <- (apply(qte_mat, 2, function(b) stats::quantile(b, .75, type = 1)) -
                  apply(qte_mat, 2, function(b) stats::quantile(b, .25, type = 1))) /
    (stats::qnorm(.75) - stats::qnorm(.25))
  if (any(sigmahalf == 0)) {
    sigmahalf <- pmax(apply(qte_mat, 2, stats::sd), 1e-9)
  }
  cb <- apply(qte_mat, 1, function(q) max(abs((q - est$qte) / sigmahalf)))
  c_val <- as.numeric(stats::quantile(cb, 1 - alp, type = 1))

  # Uniform band qte +/- c * se, exactly as ptetools builds lower_ub/upper_ub
  # for the panel estimators (same c, scaled by the bootstrap SE); pointwise
  # band qte +/- z * se otherwise
  band_cval <- if (cband) c_val else z

  list(
    qte.se = qte.se,
    qte.lower = est$qte - band_cval * qte.se,
    qte.upper = est$qte + band_cval * qte.se,
    ate.se = ate.se,
    ate.lower = est$ate - z * ate.se, ate.upper = est$ate + z * ate.se,
    c = c_val
  )
}


# --- User-facing function -----------------------------------------------------

#' @title unc_qte
#'
#' @description Estimates the Quantile Treatment Effect (QTE) or Quantile
#'   Treatment Effect on the Treated (QTT) under unconfoundedness, also known
#'   as selection on observables, with cross-sectional data. The key
#'   identifying assumption is \eqn{(Y(0), Y(1)) \perp D \mid X}, i.e.,
#'   potential outcomes are independent of treatment conditional on
#'   covariates \eqn{X}. For staggered treatment adoption with panel or
#'   repeated cross-section data, see \code{\link{lou_qtt}}.
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
#'     \item{\code{"reg"}}{Outcome regression via quantile regression inversion.
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
#'   When \code{xformula = ~1} (the default), all three methods reduce to
#'   simple quantile differences between the treated and untreated groups, which
#'   is consistent under unconditional unconfoundedness. Covariates are required
#'   for covariate-adjusted estimation.
#'
#'   Rows with a missing value in the outcome, treatment, weights or any
#'   covariate in \code{xformula} are dropped, with a warning.
#'
#'   Standard errors and confidence intervals are computed via the empirical
#'   bootstrap, resampling rows of \code{data}, or whole clusters when
#'   \code{clustervars} is given.
#'
#' @param yname character; name of the outcome variable in \code{data}.
#' @param dname character; name of the binary treatment indicator in
#'   \code{data} (1 = treated, 0 = untreated).
#' @param data data.frame containing the analysis data.
#' @param xformula one-sided formula for covariates used in the propensity
#'   score and/or outcome model, e.g. \code{~ age + I(age^2) + education}.
#'   Default \code{~1} uses no covariates (simple quantile differences).
#' @param weightsname character; name of a column in \code{data} containing
#'   sampling weights. Default \code{NULL} applies equal weights.
#' @param clustervars character; name of a column in \code{data} to cluster
#'   standard errors on (one variable). The bootstrap then resamples whole
#'   clusters instead of rows. Default \code{NULL} treats rows as
#'   independent.
#' @param probs numeric vector of quantile levels at which to evaluate the
#'   QTE or QTT. Default \code{seq(0.05, 0.95, 0.05)}.
#' @param alp significance level for confidence intervals. Default 0.05.
#' @param biters number of bootstrap iterations for standard errors and
#'   confidence bands. Default 100.
#' @param cband logical; if \code{TRUE} (default), \code{qte.lower} and
#'   \code{qte.upper} form a uniform confidence band, simultaneous over all
#'   quantiles in \code{probs}: \code{qte} \eqn{\pm} \code{c} \eqn{\times}
#'   \code{qte.se}, with the bootstrap critical value \code{c} (the same
#'   construction as the panel estimators' uniform bands). If \code{FALSE},
#'   they are pointwise intervals, \code{qte} \eqn{\pm z_{1-\alpha/2}}
#'   \code{qte.se}.
#' @param boot_type bootstrap variant. Currently only \code{"empirical"} is
#'   supported.
#' @param method propensity score link function: \code{"logit"} (default) or
#'   \code{"probit"}. Used by \code{est_method = "ipw"} and \code{"aipw"}.
#' @param est_method estimation method; one of:
#'   \itemize{
#'     \item \code{"ipw"} (default) — inverse propensity weighting
#'     \item \code{"reg"} — outcome regression via quantile regression
#'       inversion (\code{"or"} is a deprecated alias)
#'     \item \code{"aipw"} — doubly-robust augmented IPW
#'   }
#' @param cond_dist For \code{est_method = "reg"} or \code{"aipw"}, how the
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
#' @param xformla Deprecated; use \code{xformula}. Still accepted, with a
#'   warning.
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
#'     \item{\code{qte.lower}, \code{qte.upper}}{confidence bounds for
#'       \code{qte}: a uniform band if \code{cband = TRUE}, pointwise
#'       intervals otherwise (see \code{cband}).}
#'     \item{\code{cband}}{whether \code{qte.lower}/\code{qte.upper} are a
#'       uniform band.}
#'     \item{\code{ate.se}, \code{ate.lower}, \code{ate.upper}}{SE and CI
#'       for \code{ate}.}
#'     \item{\code{c}}{bootstrap critical value for a uniform band over
#'       \code{probs}.}
#'     \item{\code{probs}}{the \code{probs} vector passed by the user.}
#'     \item{\code{pscore.reg}}{fitted propensity score \code{glm} object
#'       (without copies of the data: \code{coef()}, \code{vcov()},
#'       \code{summary()}, \code{fitted()} and \code{predict(newdata = )}
#'       work; \code{anova()} and \code{update()} do not), or
#'       \code{NULL} when \code{est_method = "reg"}.}
#'     \item{\code{F.treated.t}, \code{F.treated.t.cf}}{ECDFs of the
#'       estimated quantiles of \eqn{Y(1)} and \eqn{Y(0)} at \code{probs}.}
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
#'               xformula = xf, est_method = "reg",
#'               biters = 20, probs = seq(0.05, 0.95, 0.05))
#' summary(q2)
#' }
#'
#' ## AIPW with covariates, QTT
#' \donttest{
#' q3 <- unc_qte(yname = "re78", dname = "treat", data = lalonde.psid,
#'               xformula = xf, est_method = "aipw", target = "qtt",
#'               biters = 20, probs = seq(0.05, 0.95, 0.05))
#' summary(q3)
#' }
#'
#' @inheritSection cic Discrete outcomes
#' @export
unc_qte <- function(yname, dname, data,
                    xformula = ~1,
                    weightsname = NULL,
                    clustervars = NULL,
                    probs = seq(0.05, 0.95, 0.05),
                    alp = 0.05,
                    biters = 100,
                    cband = TRUE,
                    boot_type = "empirical",
                    method = c("logit", "probit"),
                    est_method = c("ipw", "reg", "aipw"),
                    target = c("qte", "qtt"),
                    cond_dist = c("qr", "dr"),
                    cl = 1,
                    xformla = NULL) {
  # xformla was renamed xformula in 2.1.0, matching the panel estimators
  if (!is.null(xformla)) {
    if (!missing(xformula)) {
      stop("Use `xformula` only; `xformla` is its deprecated old name.")
    }
    warning("`xformla` is deprecated; use `xformula` instead.", call. = FALSE)
    xformula <- xformla
  }
  method     <- match.arg(method)
  est_method <- match.arg(resolve_est_method_or(est_method),
                          c("ipw", "reg", "aipw"))
  target     <- match.arg(target)
  cond_dist  <- match.arg(cond_dist)
  boot_type  <- match.arg(boot_type)

  # The treatment must be a 0/1 indicator: it is the response of the
  # propensity-score glm and is used to split treated/untreated units. Other
  # codings (e.g. 1/2) otherwise fail deep inside glm() with an opaque
  # "y values must be 0 <= y <= 1" error (GitHub issue #1).
  if (!dname %in% names(data)) {
    stop("`dname` (\"", dname, "\") is not a column of `data`.")
  }
  d_vals <- unique(data[[dname]])
  d_vals <- d_vals[!is.na(d_vals)]
  if (!is.numeric(d_vals) && !is.logical(d_vals) ||
        !all(d_vals %in% c(0, 1)) || length(d_vals) != 2) {
    stop("`dname` (\"", dname, "\") must be a binary treatment indicator ",
         "coded 0 (untreated) / 1 (treated); found values: ",
         paste(utils::head(sort(d_vals), 5), collapse = ", "), ".")
  }

  # warn once if the outcome looks discrete (#17): QTEs are quantile
  # differences, which degenerate at mass points
  disc <- check_discrete_outcome(as.data.frame(data), yname, weightsname)

  dat <- unc_qte_data(data, yname, dname, xformula, weightsname, clustervars)

  est_args <- list(xformula = xformula, probs = probs, method = method,
                   est_method = est_method, target = target,
                   cond_dist = cond_dist)
  est  <- do.call(unc_qte_est, c(list(dat = dat, full = TRUE), est_args))
  boot <- unc_qte_boot(dat, est_args, est, biters = biters, alp = alp,
                       cband = cband, cl = cl)

  res <- QTE(qte       = est$qte,
             qte.upper = boot$qte.upper, qte.lower = boot$qte.lower,
             ate       = est$ate,
             ate.upper = boot$ate.upper, ate.lower = boot$ate.lower,
             qte.se    = boot$qte.se,    ate.se    = boot$ate.se,
             c = boot$c, alp = alp, cband = cband,
             pscore.reg     = est$pscore.reg,
             F.treated.t    = est$F.treated.t,
             F.treated.t.cf = est$F.treated.t.cf,
             probs = probs,
             type = if (target == "qtt") "On the Treated" else "Overall")
  res$outcome_discreteness <- disc
  res
}
