# =============================================================================
# Title: Distributional Difference-in-Differences (ddid / ddid2)
# Description: Implements the Callaway-Li-Oka (2018) distributional DiD
#   estimator. ddid_gt is the 2x2 workhorse used directly as a ptetools
#   attgt_fun. ddid is the unified user-facing function for staggered adoption.
#   ddid2 is a deprecated alias retained for backwards compatibility.
#   compute.ddid2 is an unexported helper kept only to support ddid2.
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

# --- 2x2 workhorse -----------------------------------------------------------

#' @title Distributional DiD: group-time estimator
#'
#' @description Computes the distributional DiD ATT and counterfactual outcome
#'   distribution for a single 2x2 (pre/post x treated/control) data subset.
#'   Serves directly as the \code{attgt_fun} argument to \code{ptetools::pte}.
#'
#'   \strong{Identification.} Under distributional parallel trends and a copula
#'   restriction (Callaway, Li, and Oka 2018), the counterfactual outcome for
#'   each control unit \eqn{j} is
#'   \deqn{kcf_j = \Delta Y_{\text{ctrl},j} + Q_{1,\text{pre}}(u_j)}
#'   where \eqn{\Delta Y_{\text{ctrl},j} = Y_{\text{post},j} - Y_{\text{pre},j}}
#'   is the observed change for control unit \eqn{j}, \eqn{u_j =
#'   F_{0,\text{pre}}(Y_{\text{pre},j})} is that unit's rank in the control
#'   pre-period distribution, and \eqn{Q_{1,\text{pre}}} is the quantile
#'   function of the treated pre-period distribution. The unconditional
#'   counterfactual distribution \eqn{F_{Y(0),\text{post}|D=1}} is then the
#'   (weighted) empirical CDF of \eqn{\{kcf_j\}}.
#'
#'   Unlike CiC, QDiD, and MDiD, the counterfactual is indexed over
#'   \emph{control} units, not treated units. Consequently \code{F0} and the
#'   ATT counterfactual term are weighted by \code{w_pre_ctrl}, and no
#'   individual treatment effect distribution (\code{Fte}) is returned.
#'
#'   \strong{Panel data required.} The estimator needs the actual change
#'   \eqn{\Delta Y_{\text{ctrl},j}} for each control unit, which requires
#'   observing the same units in both periods.
#'
#' @param gt_data A data frame (typically a \code{gt_data_frame} from
#'   \code{ptetools}) with columns \code{name} (\code{"pre"} or \code{"post"}),
#'   \code{D} (treatment dummy), \code{Y} (outcome), \code{id} (unit
#'   identifier), \code{.w} (sampling weights), and any covariate columns
#'   referenced by \code{xformula}. Control units must be observed in both
#'   periods.
#' @param xformula One-sided formula for covariates. Default \code{~1} uses no
#'   covariates. With covariates, the unconditional rank \eqn{u_j} is replaced
#'   by a conditional rank estimated via quantile regression on the control
#'   pre-period (\code{QR0tmin1}), and the treated pre-period quantile is
#'   replaced by a conditional quantile (\code{QR1tmin1}) evaluated at that
#'   rank and the control unit's own covariate values.
#' @param ... Additional arguments passed through by \code{ptetools}; not used
#'   directly.
#'
#' @return A \code{ptetools::attgt_noif} object with the ATT estimate and, in
#'   \code{extra_gt_returns}, \code{F0} (weighted ECDF of counterfactual
#'   outcomes indexed over control units), \code{F1} (weighted ECDF of observed
#'   treated post-period outcomes), and \code{Fte = NULL} (individual treatment
#'   effect distribution is not identified for this estimator).
#'
#' @references
#'   Callaway, Brantly, Tong Li, and Tatsushi Oka. ``Quantile Treatment Effects
#'   in Difference in Differences Models under Dependence Restrictions and with
#'   Only Two Time Periods.'' Journal of Econometrics 206(2), pp. 395-413, 2018.
#'
#' @importFrom quantreg rq
#'
#' @export
ddid_gt <- function(gt_data, xformula = ~1, ...) {
  gt_data <- droplevels(gt_data)

  # extract and sort the control group by id so pre[i] and post[i] are the
  # same unit — needed to compute within-unit changes delta_ctrl
  ctrl_pre  <- gt_data[gt_data$name == "pre"  & gt_data$D == 0, ]
  ctrl_post <- gt_data[gt_data$name == "post" & gt_data$D == 0, ]
  ctrl_pre  <- ctrl_pre[order(ctrl_pre$id), ]
  ctrl_post <- ctrl_post[order(ctrl_post$id), ]

  if (!identical(ctrl_pre$id, ctrl_post$id)) {
    stop("ddid_gt requires panel data: control units must appear in both pre and post periods")
  }

  # extract treated groups (pre for rank mapping, post for ATT)
  trt_pre  <- gt_data[gt_data$name == "pre"  & gt_data$D == 1, ]
  trt_post <- gt_data[gt_data$name == "post" & gt_data$D == 1, ]

  # guard: ptetools may call attgt_fun for (g,t) cells with no treated obs
  if (nrow(trt_pre) == 0 || nrow(trt_post) == 0) {
    return(ptetools::attgt_noif(attgt = NA))
  }

  # uppercase Y and distribution objects follow econometric notation # nolint: object_name_linter
  Y_pre_ctrl  <- ctrl_pre$Y  # nolint: object_name_linter
  Y_post_ctrl <- ctrl_post$Y # nolint: object_name_linter
  Y_pre_trt   <- trt_pre$Y   # nolint: object_name_linter
  Y_post_trt  <- trt_post$Y  # nolint: object_name_linter

  # normalise sampling weights within each (period x treatment) cell
  w_pre_ctrl  <- ctrl_pre$.w  / sum(ctrl_pre$.w)  # nolint: object_name_linter
  w_post_ctrl <- ctrl_post$.w / sum(ctrl_post$.w) # nolint: object_name_linter
  w_pre_trt   <- trt_pre$.w   / sum(trt_pre$.w)   # nolint: object_name_linter
  w_post_trt  <- trt_post$.w  / sum(trt_post$.w)  # nolint: object_name_linter

  # observed within-unit changes for each control unit (requires panel)
  delta_ctrl <- Y_post_ctrl - Y_pre_ctrl # nolint: object_name_linter

  # Helper: weighted type-1 quantile (smallest y with weighted CDF >= p).
  wquant <- function(y_vec, w_vec, probs) {
    ord <- order(y_vec)
    y_s <- y_vec[ord]
    cdf <- cumsum(w_vec[ord])
    y_s[vapply(probs, function(p) {
      idx <- which(cdf >= p)
      if (length(idx) == 0L) length(y_s) else idx[1L]
    }, integer(1))]
  }

  # DDiD transform (no covariates):
  # u_j   = weighted rank of control unit j in the control pre-period distribution
  # kcf_j = delta_ctrl_j + Q_{1,pre}^w(u_j)
  # NOTE: kcf is indexed over CONTROL units; F0 and ATT counterfactual
  # term are weighted by w_pre_ctrl, not w_pre_trt.
  u <- vapply(Y_pre_ctrl,
              function(y) sum(w_pre_ctrl * (Y_pre_ctrl <= y)),
              numeric(1))

  kcf <- delta_ctrl + wquant(Y_pre_trt, w_pre_trt, u)

  att <- weighted.mean(Y_post_trt, w_post_trt) - weighted.mean(kcf, w_pre_ctrl)
  F0  <- BMisc::weighted_ecdf(kcf,        weights = w_pre_ctrl) # nolint: object_name_linter
  F1  <- BMisc::weighted_ecdf(Y_post_trt, weights = w_post_trt) # nolint: object_name_linter

  # covariate adjustment: replace unconditional rank u_j with a conditional
  # rank from QR0tmin1 and unconditional Q_{1,pre}(u_j) with a conditional
  # quantile from QR1tmin1, both evaluated at the control unit's covariates.
  # Sampling weights are passed to rq().
  if (length(BMisc::rhs.vars(xformula)) > 0) {
    u_seq    <- seq(0.01, 0.99, 0.01)
    n0       <- nrow(ctrl_pre)
    yformla  <- BMisc::toformula("Y", BMisc::rhs.vars(xformula)) # nolint: object_name_linter

    QR0tmin1 <- rq(yformla, data = ctrl_pre, tau = u_seq, weights = ctrl_pre$.w)  # nolint: object_name_linter
    QR1tmin1 <- rq(yformla, data = trt_pre,  tau = u_seq, weights = trt_pre$.w)   # nolint: object_name_linter

    # conditional rank of each control unit in the control pre-period distribution
    QR0tmin1F  <- predict(QR0tmin1, type = "Fhat", stepfun = TRUE)                # nolint: object_name_linter
    rank0tmin1 <- sapply(seq_len(n0), function(i) QR0tmin1F[[i]](ctrl_pre$Y[i]))

    # conditional treated pre-period quantile at each control unit's rank and X
    ytmin1 <- predict(QR1tmin1, newdata = ctrl_pre, type = "Qhat", stepfun = TRUE) # nolint: object_name_linter
    ytmin1 <- sapply(seq_len(n0), function(i) ytmin1[[i]](rank0tmin1[i]))

    kcf <- delta_ctrl + ytmin1

    att <- weighted.mean(Y_post_trt, w_post_trt) - weighted.mean(kcf, w_pre_ctrl)
    F0  <- BMisc::weighted_ecdf(kcf, weights = w_pre_ctrl) # nolint: object_name_linter
  }

  # Fte (individual TE distribution) is not identified: kcf is indexed over
  # control units, so there is no unit-level pairing with Y_post_trt.
  ptetools::attgt_noif(attgt = att, extra_gt_returns = list(F0 = F0, F1 = F1, Fte = NULL))
}


# --- Unified user-facing function --------------------------------------------

#' @title Distributional Difference-in-Differences
#'
#' @description Computes Quantile Treatment effects on the Treated (QTT) and
#'   the Average Treatment Effect on the Treated (ATT) using the distributional
#'   DiD identification strategy of Callaway, Li, and Oka (2018). Handles
#'   two-period data and staggered treatment adoption uniformly via
#'   \code{ptetools}. Requires panel data.
#'
#'   \strong{Identification.} Under distributional parallel trends and a copula
#'   restriction (the rank correlation of untreated potential outcomes between
#'   the treated and control groups in the pre-period is preserved), the
#'   counterfactual distribution \eqn{F_{Y(0),\text{post}|D=1}} is recovered
#'   by adding each control unit's actual change \eqn{\Delta Y_{\text{ctrl},j}}
#'   to the treated pre-period quantile at that control unit's rank. This
#'   point-identifies the QTT with only two time periods, in contrast to
#'   methods (e.g. \code{panel_qtt}) that require three periods.
#'
#'   \strong{Panel data required.} Unlike \code{cic}, \code{qdid}, and
#'   \code{mdid}, this estimator cannot be applied to repeated cross sections
#'   because it requires observing the within-unit outcome change for each
#'   control unit.
#'
#' @param yname Name of the outcome variable in \code{data}.
#' @param gname Name of the treatment group variable (first treatment period;
#'   0 for never-treated units).
#' @param tname Name of the time period variable.
#' @param idname Name of the unit id variable (required).
#' @param data A data frame.
#' @param xformula One-sided formula for covariates. Default \code{~1} uses no
#'   covariates.
#' @param weightsname Name of the column in \code{data} containing sampling
#'   weights. Default \code{NULL} uses equal weights.
#' @param control_group Which units to use as the comparison group:
#'   \code{"notyettreated"} (default) or \code{"nevertreated"}.
#' @param anticipation Number of periods of anticipation. Default \code{0}.
#' @param alp Significance level for confidence bands. Default \code{0.05}.
#' @param cband Logical; if \code{TRUE} (default) compute a uniform confidence
#'   band rather than pointwise intervals.
#' @param biters Number of bootstrap iterations. Default \code{100}.
#' @param cl Number of clusters for parallel computation. Default \code{1}.
#' @param gt_type Type of group-time effect to compute. \code{"att"} (default)
#'   returns ATT(g,t). \code{"qtt"} returns the full QTT curve over \code{probs}
#'   using mixture-CDF aggregation.
#' @param probs For \code{gt_type = "qtt"}, the quantile grid. Default is
#'   \code{seq(0.05, 0.95, 0.05)}.
#'
#' @return For \code{gt_type = "att"}, a \code{pte_results} object from
#'   \code{ptetools}. For \code{gt_type = "qtt"}, a \code{pte_qtt} object.
#'
#' @references
#'   Callaway, Brantly, Tong Li, and Tatsushi Oka. ``Quantile Treatment Effects
#'   in Difference in Differences Models under Dependence Restrictions and with
#'   Only Two Time Periods.'' Journal of Econometrics 206(2), pp. 395-413, 2018.
#'
#' @export
ddid <- function(yname,
                 gname,
                 tname,
                 idname,
                 data,
                 xformula      = ~1,
                 weightsname   = NULL,
                 control_group = "notyettreated",
                 anticipation  = 0,
                 alp           = 0.05,
                 cband         = TRUE,
                 biters        = 100,
                 cl            = 1,
                 gt_type       = "att",
                 probs         = NULL) {

  ptetools::pte(
    yname         = yname,
    gname         = gname,
    tname         = tname,
    idname        = idname,
    data          = data,
    panel         = TRUE,
    setup_pte_fun = ptetools::setup_pte,
    subset_fun    = ptetools::two_by_two_subset,
    attgt_fun     = ddid_gt,
    xformula      = xformula,
    weightsname   = weightsname,
    control_group = control_group,
    anticipation  = anticipation,
    cband         = cband,
    alp           = alp,
    boot_type     = "empirical",
    biters        = biters,
    cl            = cl,
    gt_type       = gt_type,
    probs         = probs
  )
}


# --- Deprecated two-period workhorse (supports ddid2 only) ------------------

# Kept only to support the deprecated ddid2() wrapper below. Do not use in
# new code; use ddid_gt() instead.
#' @keywords internal
compute.ddid2 <- function(qp) { # nolint: object_name_linter
  setupData(qp)
  bootstrapiter <- qp$bootstrapiter

  att <- mean(treated.t[, yname]) - mean(treated.tmin1[, yname]) - # nolint: object_usage_linter
    (mean(untreated.t[, yname]) - mean(untreated.tmin1[, yname]))   # nolint: object_usage_linter

  quantys1 <- untreated.change.t # nolint: object_usage_linter

  quantys2 <- stats::quantile(F.treated.tmin1,  # nolint: object_usage_linter
    probs = F.untreated.tmin1(untreated.tmin1[, yname]), # nolint: object_usage_linter
    type = 1
  )

  # nolint start: object_name_linter, object_usage_linter
  qr0.reg <- NULL
  u <- seq(.01, .99, .01)
  if (!(is.null(x))) {
    n0 <- nrow(untreated.t)
    yformla <- toformula("y", BMisc::rhs.vars(xformla))
    QR1tmin1 <- rq(yformla, data = treated.tmin1, tau = u)
    QR0tmin1 <- rq(yformla, data = untreated.tmin1, tau = u)
    rank0tmin1 <- predict(QR0tmin1, type = "Fhat", stepfun = TRUE)
    rank0tmin1 <- sapply(1:n0, function(i) rank0tmin1[[i]](untreated.tmin1$y[i]))
    ytmin1 <- predict(QR1tmin1, newdata = untreated.tmin1, type = "Qhat", stepfun = TRUE)
    ytmin1 <- sapply(1:n0, function(i) ytmin1[[i]](rank0tmin1[i]))
    quantys2 <- ytmin1

    att <- mean(treated.t$y) -
      sum(quantile(F.treated.t.cf, probs = u, type = 1)) / length(u)
  }
  # nolint end: object_name_linter, object_usage_linter

  y.seq <- unique((quantys1 + quantys2)[order(quantys1 + quantys2)])
  F.treated.t.cf.val <- vapply(y.seq,   # nolint: object_name_linter
    FUN = function(y) mean(1 * (quantys1 + quantys2 <= y)),
    FUN.VALUE = 1
  )
  F.treated.t.cf <- makeDist(y.seq, F.treated.t.cf.val) # nolint: object_name_linter
  F.treated.t    <- stats::ecdf(treated.t[, yname])     # nolint: object_name_linter, object_usage_linter

  qte <- stats::quantile(F.treated.t, probs = probs) - # nolint: object_usage_linter
    stats::quantile(F.treated.t.cf, probs = probs)

  if (!is.null(x)) { # nolint: object_usage_linter
    att <- mean(treated.t$y) - # nolint: object_usage_linter
      sum(quantile(F.treated.t.cf, probs = u, type = 1)) / length(u)
  }

  out <- QTE(
    F.treated.t       = F.treated.t,
    F.treated.tmin1   = F.treated.tmin1,   # nolint: object_usage_linter
    F.untreated.change.t = F.untreated.change.t, # nolint: object_usage_linter
    F.untreated.t     = F.untreated.t,     # nolint: object_usage_linter
    F.untreated.tmin1 = F.untreated.tmin1, # nolint: object_usage_linter
    F.treated.t.cf    = F.treated.t.cf,
    qte = qte, pscore.reg = NULL, ate = att, probs = probs # nolint: object_usage_linter
  )
  class(out) <- "QTE"
  out
}


# --- Deprecated two-period wrapper -------------------------------------------

#' @rdname ddid
#' @export
ddid2 <- function(formla, xformla = NULL, t, tmin1,
                  tname, data, panel = TRUE,
                  dropalwaystreated = TRUE, idname = NULL,
                  probs = seq(0.05, 0.95, 0.05),
                  iters = 100, alp = 0.05, method = "logit", se = TRUE,
                  retEachIter = FALSE, seedvec = NULL, pl = FALSE, cores = NULL,
                  biters = NULL, cl = NULL) {
  .Deprecated(msg = paste0(
    "ddid2() is deprecated. The replacement ddid() uses a different data model: ",
    "specify outcome, group, and time variables by name (yname, gname, tname) ",
    "and provide gname as a column holding each unit's first treatment period ",
    "(0 = never-treated), rather than a formula and explicit t/tmin1 values. ",
    "See ?ddid."
  ))
  if (!is.null(biters)) iters <- biters
  if (!is.null(cl)) {
    pl    <- cl > 1
    cores <- cl
  }
  if (!panel) {
    stop("method not implemented with repeated cross sections data...\n  In this case, try change in changes method...")
  }

  data <- panelize.data(data, idname, tname, t, tmin1)

  qp <- QTEparams(
    formla = formla, xformla = xformla, t = t, tmin1 = tmin1,
    tname = tname, data = data, panel = panel,
    idname = idname, probs = probs,
    iters = iters, alp = alp, method = method,
    se = se, retEachIter = retEachIter, seedvec = seedvec,
    pl = pl, cores = cores
  )

  panel.checks(qp)

  pqte <- compute.ddid2(qp)

  if (se) {
    qp$bootstrapiter <- TRUE
    SEobj <- bootstrap(qp, pqte, compute.ddid2)
    out <- QTE(
      qte = pqte$qte, qte.upper = SEobj$qte.upper,
      qte.lower = SEobj$qte.lower, ate = pqte$ate,
      ate.upper = SEobj$ate.upper, ate.lower = SEobj$ate.lower,
      qte.se = SEobj$qte.se, ate.se = SEobj$ate.se,
      c = SEobj$c, alp = alp,
      F.treated.t          = pqte$F.treated.t,
      F.untreated.t        = pqte$F.untreated.t,
      F.treated.t.cf       = pqte$F.treated.t.cf,
      F.treated.tmin1      = pqte$F.treated.tmin1,
      F.untreated.change.t = pqte$F.untreated.change.t,
      F.untreated.tmin1    = pqte$F.untreated.tmin1,
      pscore.reg           = pqte$pscore.reg,
      eachIterList         = eachIter,
      probs                = probs
    )
    return(out)
  } else {
    return(pqte)
  }
}
