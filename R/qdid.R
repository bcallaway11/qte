# =============================================================================
# Title: Quantile Difference-in-Differences (QDiD)
# Description: Implements the Quantile DiD estimator (Athey-Imbens 2006).
#   qdid_gt is the 2x2 workhorse used directly as a ptetools attgt_fun.
#   qdid is the unified user-facing function covering two-period and staggered
#   treatment adoption for both panel and repeated cross sections data.
#   (The deprecated two-period QDiD() was removed in v2.1.0.)
# Author: Brant Callaway
# Last update: 2026-11-15
# Date created: 2026-05-18
# =============================================================================

# --- 2x2 workhorse -----------------------------------------------------------

#' @title Quantile Difference-in-Differences: group-time estimator
#'
#' @description Computes the QDiD ATT and counterfactual outcome distribution
#'   for a single 2x2 (pre/post x treated/control) data subset. Serves
#'   directly as the \code{attgt_fun} argument to \code{ptetools::pte}. Panel
#'   vs. repeated cross sections is detected automatically from whether the
#'   same unit ids appear in both periods.
#'
#' @param gt_data A data frame (typically a \code{gt_data_frame} from
#'   \code{ptetools}) with columns \code{name} (\code{"pre"} or
#'   \code{"post"}), \code{D} (treatment dummy), \code{Y} (outcome), \code{id}
#'   (unit identifier), \code{.w} (sampling weights), and any covariate columns
#'   referenced by \code{xformula}.
#' @param xformula One-sided formula for covariates. Default \code{~1} uses no
#'   covariates. With covariates, conditional quantile regressions are used
#'   following Athey and Imbens (2006).
#' @param ... Additional arguments passed through by \code{ptetools}; not used
#'   directly.
#'
#' @return A \code{ptetools::attgt_noif} object with the ATT estimate and, in
#'   \code{extra_gt_returns}, three objects: \code{F1} (weighted ECDF of
#'   observed treated outcomes in the post period), \code{F0} (weighted ECDF
#'   of counterfactual untreated outcomes for the treated group), and \code{Fte}
#'   (weighted ECDF of individual treatment effects; \code{NULL} for repeated
#'   cross sections).
#'
#' @references
#'   Athey, Susan and Guido Imbens. ``Identification and Inference in Nonlinear
#'   Difference-in-Differences Models.'' Econometrica 74(2), pp. 431-497, 2006.
#'
#' @importFrom quantreg rq
#'
#' @export
qdid_gt <- function(gt_data, xformula = ~1, ...) {
  gt_data <- droplevels(gt_data)

  # extract the four (period x treatment-status) cells
  # uppercase Y and distribution objects follow econometric notation # nolint: object_name_linter
  Y_pre_ctrl  <- gt_data$Y[gt_data$name == "pre"  & gt_data$D == 0] # nolint: object_name_linter
  Y_post_ctrl <- gt_data$Y[gt_data$name == "post" & gt_data$D == 0] # nolint: object_name_linter
  Y_pre_trt   <- gt_data$Y[gt_data$name == "pre"  & gt_data$D == 1] # nolint: object_name_linter
  Y_post_trt  <- gt_data$Y[gt_data$name == "post" & gt_data$D == 1] # nolint: object_name_linter

  # extract and normalise sampling weights within each (period x treatment) cell;
  # pte_params always creates .w (= rep(1, n) when weightsname is NULL)
  w_pre_ctrl  <- gt_data$.w[gt_data$name == "pre"  & gt_data$D == 0] # nolint: object_name_linter
  w_post_ctrl <- gt_data$.w[gt_data$name == "post" & gt_data$D == 0] # nolint: object_name_linter
  w_pre_trt   <- gt_data$.w[gt_data$name == "pre"  & gt_data$D == 1] # nolint: object_name_linter
  w_post_trt  <- gt_data$.w[gt_data$name == "post" & gt_data$D == 1] # nolint: object_name_linter
  w_pre_ctrl  <- w_pre_ctrl  / sum(w_pre_ctrl)
  w_post_ctrl <- w_post_ctrl / sum(w_post_ctrl)
  w_pre_trt   <- w_pre_trt   / sum(w_pre_trt)
  w_post_trt  <- w_post_trt  / sum(w_post_trt)

  # QDiD transform (A-I 2006): rank each treated pre-period unit in its own
  # pre-period distribution, then apply the control group's DiD shift at
  # that rank.  kcf[i] = Y_pre_trt[i] + Q_{0,post}(u[i]) - Q_{0,pre}(u[i]).
  # Weighted CDFs/quantiles are used throughout.

  # weighted rank in the treated pre-period distribution
  u <- wcdf_at(Y_pre_trt, w_pre_trt, Y_pre_trt)

  Q_pre_ctrl  <- wquant(Y_pre_ctrl,  w_pre_ctrl,  u) # nolint: object_name_linter
  Q_post_ctrl <- wquant(Y_post_ctrl, w_post_ctrl, u) # nolint: object_name_linter

  kcf <- Y_pre_trt + Q_post_ctrl - Q_pre_ctrl

  att <- weighted.mean(Y_post_trt, w_post_trt) - weighted.mean(kcf, w_pre_trt)
  F0  <- BMisc::weighted_ecdf(kcf,        weights = w_pre_trt)  # nolint: object_name_linter
  F1  <- BMisc::weighted_ecdf(Y_post_trt, weights = w_post_trt) # nolint: object_name_linter

  # detect panel: same unit ids appear in both periods for the treated group
  trt_pre_ids  <- gt_data$id[gt_data$D == 1 & gt_data$name == "pre"]
  trt_post_ids <- gt_data$id[gt_data$D == 1 & gt_data$name == "post"]
  panel <- length(intersect(trt_pre_ids, trt_post_ids)) > 0

  # distribution of individual treatment effects (panel only)
  Fte <- NULL # nolint: object_name_linter
  if (panel) {
    pre_s  <- gt_data[gt_data$name == "pre"  & gt_data$D == 1, ]
    post_s <- gt_data[gt_data$name == "post" & gt_data$D == 1, ]
    pre_s  <- pre_s[order(pre_s$id), ]
    post_s <- post_s[order(post_s$id), ]
    u_s       <- wcdf_at(Y_pre_trt, w_pre_trt, pre_s$Y)
    kcf_s     <- pre_s$Y + wquant(Y_post_ctrl, w_post_ctrl, u_s) -
                            wquant(Y_pre_ctrl,  w_pre_ctrl,  u_s)
    w_post_s  <- post_s$.w / sum(post_s$.w)
    Fte <- BMisc::weighted_ecdf(post_s$Y - kcf_s, weights = w_post_s) # nolint: object_name_linter
  }

  # covariate adjustment via conditional quantile regression (Athey-Imbens 2006).
  # Conditional rank is estimated from the treated pre-period QR (QR1tmin1).
  # Sampling weights are passed to rq() for consistency with the weighted ATT.
  if (length(BMisc::rhs_vars(xformula)) > 0) {
    u_seq <- seq(0.01, 0.99, 0.01)

    pre_ctrl  <- gt_data[gt_data$name == "pre"  & gt_data$D == 0, ]
    post_ctrl <- gt_data[gt_data$name == "post" & gt_data$D == 0, ]
    pre_trt   <- gt_data[gt_data$name == "pre"  & gt_data$D == 1, ]
    post_trt  <- gt_data[gt_data$name == "post" & gt_data$D == 1, ]

    if (panel) {
      pre_trt  <- pre_trt[order(pre_trt$id), ]
      post_trt <- post_trt[order(post_trt$id), ]
    }

    yformla <- BMisc::toformula("Y", BMisc::rhs_vars(xformula)) # nolint: object_name_linter

    # conditional rank in the treated pre-period distribution
    QR1tmin1  <- suppressWarnings(rq(yformla, data = pre_trt,   tau = u_seq, weights = .w)) # nolint: object_name_linter
    # post-period control fit: the same for all cohorts treated by this
    # period, so reused across those cells (see cached_rq())
    QR0t      <- cached_rq(yformla, post_ctrl, u_seq) # nolint: object_name_linter
    QR0tmin1  <- suppressWarnings(rq(yformla, data = pre_ctrl,  tau = u_seq, weights = .w)) # nolint: object_name_linter

    # rearranged (monotone) conditional quantiles; see qr_sorted_quantiles()
    # for why predict(..., stepfun = TRUE) isn't used
    rank1tmin1 <- qr_cdf(qr_sorted_quantiles(QR1tmin1, pre_trt), u_seq, pre_trt$Y)
    QR0tQ     <- qr_quantile(qr_sorted_quantiles(QR0t, pre_trt), u_seq, rank1tmin1)     # nolint: object_name_linter
    QR0tmin1Q <- qr_quantile(qr_sorted_quantiles(QR0tmin1, pre_trt), u_seq, rank1tmin1) # nolint: object_name_linter

    y0t <- pre_trt$Y + QR0tQ - QR0tmin1Q

    att          <- weighted.mean(Y_post_trt, w_post_trt) - weighted.mean(y0t, w_pre_trt)
    w_pre_trt_n  <- pre_trt$.w / sum(pre_trt$.w)
    F0           <- BMisc::weighted_ecdf(y0t, weights = w_pre_trt_n) # nolint: object_name_linter

    if (panel) {
      w_post_trt_n <- post_trt$.w / sum(post_trt$.w)
      Fte <- BMisc::weighted_ecdf(post_trt$Y - y0t, weights = w_post_trt_n) # nolint: object_name_linter
    }
  }

  ptetools::attgt_noif(attgt = att, extra_gt_returns = list(F0 = F0, F1 = F1, Fte = Fte))
}


# --- Unified user-facing function --------------------------------------------

#' @title Quantile Difference-in-Differences
#'
#' @description Computes Quantile Treatment effects on the Treated (QTT) and
#'   the Average Treatment Effect on the Treated (ATT) using the Quantile
#'   Difference-in-Differences identification strategy of Athey and Imbens
#'   (2006). Handles two-period data and staggered treatment adoption
#'   uniformly. Supports both panel and repeated cross sections data.
#'
#' @param yname Name of the outcome variable in \code{data}.
#' @param gname Name of the treatment group variable (first treatment period;
#'   0 for never-treated units).
#' @param tname Name of the time period variable.
#' @param idname Name of the unit id variable. Required when
#'   \code{panel = TRUE}.
#' @param data A data frame.
#' @param panel Logical; \code{TRUE} (default) for panel data,
#'   \code{FALSE} for repeated cross sections.
#' @param xformula One-sided formula for covariates used in the covariate
#'   adjustment. Default \code{~1} uses no covariates.
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
#' @param ret_quantile Passed through to \code{ptetools} for the \code{"qott"}
#'   case.
#' @param gt_type Type of group-time effect to compute. \code{"att"} (default)
#'   returns ATT(g,t). \code{"qtt"} returns the full QTT curve over \code{probs}
#'   using mixture-CDF aggregation. \code{"qott"} returns the quantile of the
#'   individual treatment effect distribution (panel only).
#' @param probs For \code{gt_type = "qtt"}, the quantile grid at which to
#'   evaluate the QTT curve. Default is \code{seq(0.05, 0.95, 0.05)}.
#'
#' @return For \code{gt_type = "att"}, a \code{pte_results} object from
#'   \code{ptetools}. For \code{gt_type = "qtt"}, a \code{pte_qtt} object
#'   with overall, group-specific, and dynamic QTT curves and bootstrap SEs.
#'
#' @references
#'   Athey, Susan and Guido Imbens. ``Identification and Inference in Nonlinear
#'   Difference-in-Differences Models.'' Econometrica 74(2), pp. 431-497, 2006.
#'
#' @examples
#' \donttest{
#' data(mpdta, package = "did")
#'
#' ## ATT aggregated across all groups and periods
#' res_att <- qdid(yname = "lemp", gname = "first.treat", tname = "year",
#'                 idname = "countyreal", data = mpdta,
#'                 gt_type = "att", biters = 20)
#' summary(res_att)
#'
#' ## Full QTT curve at selected quantiles
#' res_qtt <- qdid(yname = "lemp", gname = "first.treat", tname = "year",
#'                 idname = "countyreal", data = mpdta,
#'                 gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 20)
#' summary(res_qtt)
#' }
#'
#' @export
qdid <- function(yname,
                 gname,
                 tname,
                 idname        = NULL,
                 data,
                 panel         = TRUE,
                 xformula      = ~1,
                 weightsname   = NULL,
                 control_group = "notyettreated",
                 anticipation  = 0,
                 alp           = 0.05,
                 cband         = TRUE,
                 biters        = 100,
                 cl            = 1,
                 ret_quantile  = NULL,
                 gt_type       = "att",
                 probs         = NULL) {

  if (gt_type == "qott" && !isTRUE(panel)) {
    stop(paste0(
      "gt_type = 'qott' requires panel = TRUE: the distribution of individual ",
      "treatment effects is not identified for repeated cross sections"
    ))
  }

  if (isTRUE(panel) && is.null(idname)) {
    stop("idname must be provided when panel = TRUE")
  }

  subset_fun <- if (isTRUE(panel)) {
    ptetools::two_by_two_subset
  } else {
    ptetools::two_by_two_rcs_subset
  }

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
                                        weightsname, xformula),
    panel           = panel,
    setup_pte_fun   = ptetools::setup_pte,
    subset_fun      = subset_fun,
    attgt_fun       = qdid_gt,
    aggte_fun       = aggte_fun,
    xformula        = xformula,
    weightsname     = weightsname,
    control_group   = control_group,
    anticipation    = anticipation,
    cband           = cband,
    alp             = alp,
    boot_type       = "empirical",
    biters          = biters,
    cl              = cl,
    ret_quantile    = ret_quantile,
    gt_type         = gt_type,
    probs           = probs
  )
}
