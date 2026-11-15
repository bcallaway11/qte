# =============================================================================
# Title: Change in Changes estimator
# Description: Implements the Athey-Imbens (2006) Change in Changes estimator.
#   cic_gt is the 2x2 workhorse used directly as a ptetools attgt_fun.
#   cic is the unified user-facing function covering two-period and staggered
#   treatment adoption for both panel and repeated cross sections data.
#   (The deprecated two-period CiC() was removed in v2.1.0.)
# Author: Brant Callaway
# Last update: 2026-11-15
# Date created: 2026-05-18
# =============================================================================

# --- 2x2 workhorse -----------------------------------------------------------

#' @title Change in Changes: group-time estimator
#'
#' @description Computes the Change in Changes (CiC) ATT and counterfactual
#'   outcome distribution for a single 2x2 (pre/post x treated/control) data
#'   subset. Serves directly as the \code{attgt_fun} argument to
#'   \code{ptetools::pte}. Panel vs. repeated cross sections is detected
#'   automatically from whether the same unit ids appear in both periods.
#'
#' @param gt_data A data frame (typically a \code{gt_data_frame} from
#'   \code{ptetools}) with columns \code{name} (\code{"pre"} or
#'   \code{"post"}), \code{D} (treatment dummy), \code{Y} (outcome), \code{id}
#'   (unit identifier), and any covariate columns referenced by
#'   \code{xformula}.
#' @param xformula One-sided formula for covariates. Default \code{~1} uses no
#'   covariates. With covariates, the control group's conditional outcome
#'   distributions are estimated by the method in \code{cond_dist}, following
#'   Athey and Imbens (2006).
#' @param cond_dist How the conditional distributions are estimated when
#'   \code{xformula} has covariates: \code{"qr"} (default) for quantile
#'   regression at the percentiles 0.01, ..., 0.99, or \code{"dr"} for
#'   distribution regression (a logit of \eqn{1\{Y \le y\}} on the covariates
#'   at a grid of thresholds; Chernozhukov, Fernandez-Val and Melly 2013).
#'   Distribution regression handles outcomes with mass points (e.g. many
#'   zeros) naturally. Ignored without covariates.
#' @param ... Additional arguments passed through by \code{ptetools}; not used
#'   directly.
#'
#' @return A \code{ptetools::attgt_noif} object with the ATT estimate and, in
#'   \code{extra_gt_returns}, three objects: \code{F1} (ECDF of observed
#'   treated outcomes in the post period), \code{F0} (ECDF of counterfactual
#'   untreated outcomes for the treated group), and \code{Fte} (ECDF of
#'   individual treatment effects under rank invariance; \code{NULL} for
#'   repeated cross sections).
#'
#' @references
#'   Athey, Susan and Guido Imbens. ``Identification and Inference in Nonlinear
#'   Difference-in-Differences Models.'' Econometrica 74(2), pp. 431-497, 2006.
#'
#'   Chernozhukov, Victor, Ivan Fernandez-Val and Blaise Melly. ``Inference on
#'   Counterfactual Distributions.'' Econometrica 81(6), pp. 2205-2268, 2013.
#'
#' @export
cic_gt <- function(gt_data, xformula = ~1, cond_dist = c("qr", "dr"), ...) {
  cond_dist <- match.arg(cond_dist)
  gt_data <- droplevels(gt_data)

  # extract the four (period x treatment-status) groups
  # uppercase Y and distribution objects (F0, F1, Fte) follow econometric
  # notation from Athey-Imbens (2006) # nolint: object_name_linter
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

  # CiC transform: map each treated pre-period rank through the control
  # distributions to obtain counterfactual outcomes (A-I 2006, eq. 16).
  # Weighted CDFs are used; reduces to ecdf()/quantile(type=1) for uniform weights.

  # Weighted ECDF of ctrl pre-period: F(y) = sum(w_i * 1(Y_i <= y)),
  # evaluated at each treated pre-period outcome (handles ties and
  # out-of-range values; see wcdf_at()).
  u <- wcdf_at(Y_pre_ctrl, w_pre_ctrl, Y_pre_trt)

  # Weighted type-1 quantile of ctrl post-period: smallest y with CDF >= p.
  ord_post_c <- order(Y_post_ctrl)
  y_post_c_s <- Y_post_ctrl[ord_post_c]
  cdf_post_c <- cumsum(w_post_ctrl[ord_post_c])
  kcic <- y_post_c_s[first_geq(cdf_post_c, u)]

  att <- weighted.mean(Y_post_trt, w_post_trt) - weighted.mean(kcic, w_pre_trt)
  F0  <- BMisc::weighted_ecdf(kcic,       weights = w_pre_trt)  # nolint: object_name_linter
  F1  <- BMisc::weighted_ecdf(Y_post_trt, weights = w_post_trt) # nolint: object_name_linter

  # detect panel: same unit ids appear in both periods for the treated group
  trt_pre_ids  <- gt_data$id[gt_data$D == 1 & gt_data$name == "pre"]
  trt_post_ids <- gt_data$id[gt_data$D == 1 & gt_data$name == "post"]
  panel <- length(intersect(trt_pre_ids, trt_post_ids)) > 0

  # distribution of individual treatment effects under rank invariance
  # (panel only; not identified for repeated cross sections)
  Fte <- NULL # nolint: object_name_linter
  if (panel) {
    pre_s  <- gt_data[gt_data$name == "pre"  & gt_data$D == 1, ]
    post_s <- gt_data[gt_data$name == "post" & gt_data$D == 1, ]
    pre_s  <- pre_s[order(pre_s$id), ]
    post_s <- post_s[order(post_s$id), ]
    u_s    <- wcdf_at(Y_pre_ctrl, w_pre_ctrl, pre_s$Y)
    kcic_s <- y_post_c_s[first_geq(cdf_post_c, u_s)]
    w_post_trt_s <- post_s$.w / sum(post_s$.w)
    Fte    <- BMisc::weighted_ecdf(post_s$Y - kcic_s, weights = w_post_trt_s) # nolint: object_name_linter
  }

  # covariate adjustment via conditional distributions (Athey-Imbens 2006),
  # estimated by quantile regression (cond_dist = "qr") or distribution
  # regression ("dr"). Sampling weights enter both, so they match the weights
  # used for the ATT and CDF calculations.
  if (length(BMisc::rhs_vars(xformula)) > 0) {
    u_seq <- seq(0.01, 0.99, 0.01)

    pre_ctrl  <- gt_data[gt_data$name == "pre"  & gt_data$D == 0, ]
    post_ctrl <- gt_data[gt_data$name == "post" & gt_data$D == 0, ]
    pre_trt   <- gt_data[gt_data$name == "pre"  & gt_data$D == 1, ]
    post_trt  <- gt_data[gt_data$name == "post" & gt_data$D == 1, ]

    # for panel, sort by id so pre_trt[i] and post_trt[i] belong to the same unit
    if (panel) {
      pre_trt  <- pre_trt[order(pre_trt$id), ]
      post_trt <- post_trt[order(post_trt$id), ]
    }

    yformla  <- BMisc::toformula("Y", BMisc::rhs_vars(xformula))
    # post-period control fit: the same for all cohorts treated by this
    # period, so reused across those cells (cache = TRUE)
    fit0t    <- cond_fit(yformla, post_ctrl, cond_dist, u_seq, cache = TRUE)
    fit0tmin1 <- cond_fit(yformla, pre_ctrl, cond_dist, u_seq)

    # conditional ranks from rearranged (monotone) conditional distributions;
    # see qr_sorted_quantiles() and dr_cdf_matrix()
    F0tmin1 <- cond_cdf(fit0tmin1, pre_trt, pre_trt$Y) # nolint: object_name_linter

    if (mean(F0tmin1 >= 0.99) + mean(F0tmin1 <= 0.01) > 0.1) {
      warning(paste0(
        "many treated units have extreme ranks in the control pre-period ",
        "distribution; CiC support conditions may be violated"
      ))
    }

    y0t <- cond_quantile(fit0t, pre_trt, F0tmin1)

    att <- weighted.mean(Y_post_trt, w_post_trt) - weighted.mean(y0t, w_pre_trt)
    w_pre_trt_n  <- pre_trt$.w / sum(pre_trt$.w)
    F0  <- BMisc::weighted_ecdf(y0t, weights = w_pre_trt_n) # nolint: object_name_linter

    if (panel) {
      w_post_trt_n <- post_trt$.w / sum(post_trt$.w)
      Fte <- BMisc::weighted_ecdf(post_trt$Y - y0t, weights = w_post_trt_n) # nolint: object_name_linter
    }
  }

  ptetools::attgt_noif(attgt = att, extra_gt_returns = list(F0 = F0, F1 = F1, Fte = Fte))
}


# --- Unified user-facing function --------------------------------------------

#' @title Change in Changes
#'
#' @description Computes Quantile Treatment effects on the Treated (QTT) and
#'   the Average Treatment Effect on the Treated (ATT) using the Change in
#'   Changes identification strategy of Athey and Imbens (2006). Handles
#'   two-period data and staggered treatment adoption uniformly: a two-period,
#'   two-group dataset is the degenerate single-(g,t) case. Supports both
#'   panel and repeated cross sections data.
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
#' @inheritParams cic_gt
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
#'   case. Ignored when \code{gt_type = "qtt"} (use \code{probs} instead).
#' @param gt_type Type of group-time effect to compute. \code{"att"} (default)
#'   returns ATT(g,t). \code{"qtt"} returns the full QTT curve over \code{probs}
#'   using mixture-CDF aggregation. \code{"qott"} returns the quantile of the
#'   individual treatment effect distribution under rank invariance (panel only).
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
#' res_att <- cic(yname = "lemp", gname = "first.treat", tname = "year",
#'                idname = "countyreal", data = mpdta,
#'                gt_type = "att", biters = 20)
#' summary(res_att)
#'
#' ## Full QTT curve at selected quantiles
#' res_qtt <- cic(yname = "lemp", gname = "first.treat", tname = "year",
#'                idname = "countyreal", data = mpdta,
#'                gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 20)
#' summary(res_qtt)
#' }
#'
#' @export
cic <- function(yname,
                gname,
                tname,
                idname        = NULL,
                data,
                panel         = TRUE,
                xformula      = ~1,
                cond_dist     = c("qr", "dr"),
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
  cond_dist <- match.arg(cond_dist)

  if (gt_type == "qott" && !isTRUE(panel)) {
    stop(paste0(
      "gt_type = 'qott' requires panel = TRUE: the distribution of individual ",
      "treatment effects under rank invariance is not identified for repeated cross sections"
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
    attgt_fun       = cic_gt,
    aggte_fun       = aggte_fun,
    xformula        = xformula,
    cond_dist       = cond_dist,
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
