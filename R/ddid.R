# =============================================================================
# Title: Distributional Difference-in-Differences (ddid)
# Description: Implements the Callaway-Li-Oka (2018) distributional DiD
#   estimator. ddid_gt is the 2x2 workhorse used directly as a ptetools
#   attgt_fun. ddid is the unified user-facing function for staggered adoption.
#   (The deprecated two-period ddid2() was removed in v2.1.0.)
# Author: Brant Callaway
# Last update: 2026-11-19
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
#' @param cond_dist How the conditional distributions are estimated when
#'   \code{xformula} has covariates: \code{"qr"} (default) for quantile
#'   regression or \code{"dr"} for distribution regression. See
#'   \code{\link{cic_gt}}. Ignored without covariates.
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
ddid_gt <- function(gt_data, xformula = ~1, cond_dist = c("qr", "dr"), ...) {
  cond_dist <- match.arg(cond_dist)
  gt_data <- droplevels(gt_data)

  # extract and sort all four cells by id (panel structure enforced by ddid())
  ctrl_pre  <- gt_data[gt_data$name == "pre"  & gt_data$D == 0, ]
  ctrl_post <- gt_data[gt_data$name == "post" & gt_data$D == 0, ]
  trt_pre   <- gt_data[gt_data$name == "pre"  & gt_data$D == 1, ]
  trt_post  <- gt_data[gt_data$name == "post" & gt_data$D == 1, ]

  # guard: ptetools may call attgt_fun for (g,t) cells with no treated obs
  if (nrow(trt_pre) == 0 || nrow(trt_post) == 0) {
    return(ptetools::attgt_noif(attgt = NA))
  }

  ctrl_pre  <- ctrl_pre[order(ctrl_pre$id), ]
  ctrl_post <- ctrl_post[order(ctrl_post$id), ]
  trt_pre   <- trt_pre[order(trt_pre$id), ]
  trt_post  <- trt_post[order(trt_post$id), ]

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

  # DDiD transform (no covariates):
  # u_j   = weighted rank of control unit j in the control pre-period distribution
  # kcf_j = delta_ctrl_j + Q_{1,pre}^w(u_j)
  # NOTE: kcf is indexed over CONTROL units; F0 and ATT counterfactual
  # term are weighted by w_pre_ctrl, not w_pre_trt.
  u <- wcdf_at(Y_pre_ctrl, w_pre_ctrl, Y_pre_ctrl)

  kcf <- delta_ctrl + wquant(Y_pre_trt, w_pre_trt, u)

  att <- weighted.mean(Y_post_trt, w_post_trt) - weighted.mean(kcf, w_pre_ctrl)
  F0  <- BMisc::weighted_ecdf(kcf,        weights = w_pre_ctrl) # nolint: object_name_linter
  F1  <- BMisc::weighted_ecdf(Y_post_trt, weights = w_post_trt) # nolint: object_name_linter

  # covariate adjustment: replace unconditional rank u_j with a conditional
  # rank in the control pre-period distribution and unconditional
  # Q_{1,pre}(u_j) with a conditional quantile of the treated pre-period
  # distribution, both evaluated at the control unit's covariates. Estimated
  # by quantile or distribution regression (cond_dist), with sampling weights.
  #
  # kcf_j is a draw from the counterfactual distribution *given X_j*, so the
  # treated group's counterfactual integrates it over the treated covariate
  # distribution: control units are reweighted by the odds p(X)/(1 - p(X)),
  # p(x) = P(D = 1 | X = x) (issue #10; previously the control weights were
  # used, which integrates over the controls' X and is biased whenever
  # treatment depends on X).
  if (length(BMisc::rhs_vars(xformula)) > 0) {
    u_seq    <- seq(0.01, 0.99, 0.01)
    yformla  <- BMisc::toformula("Y", BMisc::rhs_vars(xformula)) # nolint: object_name_linter

    fit0tmin1 <- cond_fit(yformla, ctrl_pre, cond_dist, u_seq)
    fit1tmin1 <- cond_fit(yformla, trt_pre,  cond_dist, u_seq)

    # conditional rank of each control unit in the control pre-period
    # distribution (rearranged; see qr_sorted_quantiles(), dr_cdf_matrix())
    rank0tmin1 <- cond_cdf(fit0tmin1, ctrl_pre, ctrl_pre$Y)

    # conditional treated pre-period quantile at each control unit's rank and X
    ytmin1 <- cond_quantile(fit1tmin1, ctrl_pre, rank0tmin1)

    kcf <- delta_ctrl + ytmin1

    # propensity score from the pre-period covariates (weighted logit;
    # quasibinomial gives the binomial MLE without non-integer-weight warnings)
    pre_both   <- rbind(ctrl_pre, trt_pre)
    ps_formla  <- BMisc::toformula("D", BMisc::rhs_vars(xformula))
    ps_fit     <- stats::glm(ps_formla, data = pre_both,
                             family = stats::quasibinomial(), weights = .w)
    ps_ctrl    <- stats::predict(ps_fit, newdata = ctrl_pre, type = "response")
    ps_ctrl    <- pmin(pmax(ps_ctrl, 1e-6), 1 - 1e-6)
    w_ctrl_trt <- ctrl_pre$.w * ps_ctrl / (1 - ps_ctrl)
    w_ctrl_trt <- w_ctrl_trt / sum(w_ctrl_trt)

    att <- weighted.mean(Y_post_trt, w_post_trt) - weighted.mean(kcf, w_ctrl_trt)
    F0  <- BMisc::weighted_ecdf(kcf, weights = w_ctrl_trt) # nolint: object_name_linter
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
#' @param cond_dist How conditional distributions are estimated when
#'   \code{xformula} has covariates: \code{"qr"} (default) for quantile
#'   regression or \code{"dr"} for distribution regression. See
#'   \code{\link{cic_gt}}. Ignored without covariates.
#' @param weightsname Name of the column in \code{data} containing sampling
#'   weights. Default \code{NULL} uses equal weights.
#' @param clustervars Name of a column in \code{data} to cluster standard
#'   errors on, in addition to \code{idname} (as in \code{did::att_gt()}):
#'   one variable, constant over time within each unit. Default \code{NULL}
#'   clusters at the unit level. With the empirical bootstrap, whole clusters
#'   are resampled.
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
#' @examples
#' \donttest{
#' data(mpdta, package = "did")
#'
#' ## ATT aggregated across all groups and periods
#' res_att <- ddid(yname = "lemp", gname = "first.treat", tname = "year",
#'                 idname = "countyreal", data = mpdta,
#'                 gt_type = "att", biters = 20)
#' summary(res_att)
#'
#' ## Full QTT curve at selected quantiles
#' res_qtt <- ddid(yname = "lemp", gname = "first.treat", tname = "year",
#'                 idname = "countyreal", data = mpdta,
#'                 gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 20)
#' summary(res_qtt)
#' }
#'
#' @export
ddid <- function(yname,
                 gname,
                 tname,
                 idname,
                 data,
                 xformula      = ~1,
                 cond_dist     = c("qr", "dr"),
                 weightsname   = NULL,
                 clustervars   = NULL,
                 control_group = "notyettreated",
                 anticipation  = 0,
                 alp           = 0.05,
                 cband         = TRUE,
                 biters        = 100,
                 cl            = 1,
                 gt_type       = "att",
                 probs         = NULL) {
  cond_dist <- match.arg(cond_dist)

  aggte_fun <- if (gt_type == "qtt") {
    ptetools::qtt_pte_aggregations
  } else {
    function(al, p, eg) ptetools::attgt_pte_aggregations(al, p)
  }

  ptetools::pte(
    yname           = yname,
    gname           = gname,
    tname           = tname,
    idname          = idname,
    data            = keep_used_columns(data, yname, gname, tname, idname,
                                        weightsname, xformula,
                                        clustervars = clustervars),
    panel           = TRUE,
    setup_pte_fun   = ptetools::setup_pte,
    subset_fun      = ptetools::two_by_two_subset,
    attgt_fun       = ddid_gt,
    aggte_fun       = aggte_fun,
    xformula        = xformula,
    cond_dist       = cond_dist,
    weightsname     = weightsname,
    clustervars     = clustervars,
    control_group   = control_group,
    anticipation    = anticipation,
    cband           = cband,
    alp             = alp,
    boot_type       = "empirical",
    biters          = biters,
    cl              = cl,
    gt_type         = gt_type,
    probs           = probs
  )
}
