# =============================================================================
# Title: Difference-in-differences in CDFs (cdf_did)
# Description: QTT and ATT under parallel trends in the distribution functions
#   of untreated potential outcomes (Roth and Sant'Anna 2023). cdf_did_gt is
#   the 2x2 workhorse used as a ptetools attgt_fun: the counterfactual CDF is
#   F_pre,treated + F_post,comparison - F_pre,comparison, computed exactly on
#   the pooled support and rearranged into a valid CDF. With covariates, the
#   comparison-group change in CDFs is estimated threshold by threshold as a
#   DiD for 1{Y <= y}: by outcome regression (distribution regression), IPW,
#   or doubly robust AIPW (Sant'Anna and Zhao 2020). cdf_did is the
#   user-facing wrapper (panel or repeated cross sections, staggered
#   adoption).
# Author: Brant Callaway
# Last update: 2026-11-19
# Date created: 2026-11-19
# =============================================================================

# --- Internal: comparison-group change in CDFs, with covariates --------------

# cdf_did_delta
# E[F_post,comparison(y | X) - F_pre,comparison(y | X) | D = 1] at each point
# of `s`, the comparison group's change in conditional CDFs averaged over the
# treated group's covariates (taken from the treated pre-period rows). For
# each y this is the counterfactual-change term of a DiD for the outcome
# 1{Y <= y}, so the estimators of Sant'Anna and Zhao (2020) apply period by
# period, with p the propensity score and r = p / (1 - p):
#   reg          mean over treated of m_t(X; y), m_t from a distribution
#                regression on the comparison group in period t
#   ipw          sum over comparison rows of w r 1{Y <= y} / sum of w r
#   aipw         reg + (ipw - the same reweighted average of m_t),
#                consistent if either p or m_t is correctly specified
# and the post-period term minus the pre-period term. IPW is exact on the
# data's support points; the m_t use distribution regression's threshold grid.
#
# Args:
#   trt_pre, ctrl_pre, ctrl_post  cell data frames (Y, D, .w, covariates)
#   xformula, est_method          as in cdf_did_gt()
#   s                             sorted evaluation points
cdf_did_delta <- function(trt_pre, ctrl_pre, ctrl_post, xformula, est_method,
                          s) {
  xvars <- BMisc::rhs_vars(xformula)
  w_trt <- trt_pre$.w / sum(trt_pre$.w)

  # propensity score P(D = 1 | X) from the pre-period rows (weighted logit;
  # quasibinomial gives the binomial MLE without non-integer-weight
  # warnings), evaluated at each comparison row's own covariates
  odds <- NULL
  if (est_method != "reg") {
    pre_both <- rbind(trt_pre, ctrl_pre)
    ps_fit <- stats::glm(BMisc::toformula("D", xvars), data = pre_both,
                         family = stats::quasibinomial(), weights = .w)
    odds <- function(d) {
      p <- stats::predict(ps_fit, newdata = d, type = "response")
      p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
      p / (1 - p)
    }
  }

  # E[F_t(y | X) | D = 1] for one comparison-group period
  period_term <- function(ctrl) {
    if (est_method == "ipw") {
      r <- ctrl$.w * odds(ctrl)
      return(wcdf_at(ctrl$Y, r / sum(r), s))
    }
    fit <- dr_fit(BMisc::toformula("Y", xvars), ctrl)
    out <- dr_avg_cdf(fit, trt_pre, w_trt, s)
    if (est_method == "aipw") {
      r <- ctrl$.w * odds(ctrl)
      r <- r / sum(r)
      # reweighted residual 1{Y <= y} - m_t(X; y) over the comparison rows
      out <- out + wcdf_at(ctrl$Y, r, s) - dr_avg_cdf(fit, ctrl, r, s)
    }
    out
  }

  period_term(ctrl_post) - period_term(ctrl_pre)
}

# --- cdf_did_gt --------------------------------------------------------------

#' @title CDF Difference-in-Differences for a Single 2x2 Cell
#'
#' @description Computes the ATT and the distributions needed for the QTT for
#'   one group-time cell under parallel trends in the distribution functions
#'   of untreated potential outcomes. Used as the \code{attgt_fun} in
#'   \code{\link[ptetools]{pte}} by \code{\link{cdf_did}}.
#'
#' @details The counterfactual distribution of untreated outcomes for the
#'   treated group in the post period is
#'   \deqn{F_{0}(y) = F_{\text{pre},D=1}(y) + F_{\text{post},D=0}(y) -
#'   F_{\text{pre},D=0}(y),}
#'   evaluated exactly at every point of the pooled support of the three
#'   distributions on the right. The raw \eqn{F_0} need not be a valid
#'   distribution function (it can decrease, or leave [0, 1]); it is made one
#'   by monotone rearrangement (sorting its values; Chernozhukov,
#'   Fernandez-Val and Galichon 2010) and clipping to [0, 1]. The largest
#'   change this makes is returned as \code{cdf_violation}; large values are
#'   evidence against the identifying assumption.
#'
#'   The ATT is the mean of the treated group's observed post-period outcomes
#'   minus the mean of the \emph{raw} \eqn{F_0}, which is what the assumption
#'   identifies. Without covariates this equals the usual DiD ATT.
#'
#'   With covariates, the comparison group's change in CDFs is averaged over
#'   the treated group's covariate distribution. For each \eqn{y} this is a
#'   DiD for the outcome \eqn{1\{Y \le y\}}, estimated as in Sant'Anna and
#'   Zhao (2020): by outcome regression (distribution regression of
#'   \eqn{1\{Y \le y\}} on the covariates, in each period), by inverse
#'   probability weighting, or doubly robust (AIPW, the default). IPW is
#'   exact on the data's support points; outcome regression and AIPW use the
#'   distribution regression's threshold grid (see \code{\link{cic_gt}}'s
#'   \code{cond_dist = "dr"}). The propensity score is a weighted logit of
#'   treatment on the pre-period covariates.
#'
#'   The distribution of individual treatment effects is not identified, so
#'   \code{Fte} is \code{NULL}.
#'
#' @param gt_data A \code{gt_data_frame} for one (g,t) cell, as produced by
#'   \code{ptetools::two_by_two_subset} (panel) or
#'   \code{ptetools::two_by_two_rcs_subset} (repeated cross sections), with
#'   columns \code{Y}, \code{D}, \code{name}, \code{id}, \code{.w} and any
#'   covariates.
#' @param xformula One-sided formula of covariates. Default \code{~1}.
#' @param est_method With covariates: \code{"aipw"} (default, doubly
#'   robust), \code{"ipw"}, or \code{"reg"} (outcome regression). Ignored
#'   without
#'   covariates.
#' @param ... Additional arguments passed by \code{ptetools}; ignored.
#'
#' @return A \code{ptetools::attgt_noif} object with the ATT and, in
#'   \code{extra_gt_returns}, \code{F0} (rearranged counterfactual
#'   distribution), \code{F1} (treated post-period distribution),
#'   \code{Fte = NULL}, and \code{cdf_violation}.
#'
#' @references
#'   Roth, Jonathan, and Pedro H. C. Sant'Anna. ``When Is Parallel Trends
#'   Sensitive to Functional Form?'' Econometrica 91(2), pp. 737-747, 2023.
#'
#'   Sant'Anna, Pedro H. C., and Jun Zhao. ``Doubly Robust
#'   Difference-in-Differences Estimators.'' Journal of Econometrics 219(1),
#'   pp. 101-122, 2020.
#'
#'   Chernozhukov, Victor, Ivan Fernandez-Val, and Alfred Galichon.
#'   ``Quantile and Probability Curves Without Crossing.'' Econometrica
#'   78(3), pp. 1093-1125, 2010.
#'
#' @export
cdf_did_gt <- function(gt_data, xformula = ~1,
                       est_method = c("aipw", "ipw", "reg"), ...) {
  est_method <- match.arg(est_method)
  gt_data <- droplevels(gt_data)
  cell <- function(nm, d) {
    gt_data[gt_data$name == nm & gt_data$D == d, , drop = FALSE]
  }
  trt_pre   <- cell("pre", 1)
  trt_post  <- cell("post", 1)
  ctrl_pre  <- cell("pre", 0)
  ctrl_post <- cell("post", 0)
  if (min(nrow(trt_pre), nrow(trt_post), nrow(ctrl_pre),
          nrow(ctrl_post)) == 0L) {
    return(ptetools::attgt_noif(attgt = NA))
  }
  wnorm <- function(d) d$.w / sum(d$.w)

  # every value the counterfactual CDF can jump at
  s <- sort(unique(c(trt_pre$Y, ctrl_pre$Y, ctrl_post$Y)))

  # comparison group's change in CDFs, unconditional or averaged over the
  # treated group's covariates
  if (length(BMisc::rhs_vars(xformula)) > 0) {
    delta <- cdf_did_delta(trt_pre, ctrl_pre, ctrl_post, xformula,
                           est_method, s)
  } else {
    delta <- wcdf_at(ctrl_post$Y, wnorm(ctrl_post), s) -
      wcdf_at(ctrl_pre$Y, wnorm(ctrl_pre), s)
  }
  F0_raw <- wcdf_at(trt_pre$Y, wnorm(trt_pre), s) + delta # nolint: object_name_linter

  # ATT from the raw counterfactual, which the assumption identifies:
  # mean = max(s) - integral of F0 below it (any mass F0 leaves short of 1
  # at max(s), from the regression grids, is put at max(s)); without
  # covariates this is exactly the DiD of means
  K <- length(s) # nolint: object_name_linter
  mean0 <- s[K] - sum(F0_raw[-K] * diff(s))
  att <- stats::weighted.mean(trt_post$Y, trt_post$.w) - mean0

  # a valid CDF for the QTT: monotone rearrangement (sort the values;
  # Chernozhukov, Fernandez-Val and Galichon 2010), then clip to [0, 1],
  # ending at 1
  F0_cdf <- pmin(pmax(sort(F0_raw), 0), 1) # nolint: object_name_linter
  F0_cdf[K] <- 1 # nolint: object_name_linter
  F0 <- BMisc::weighted_ecdf(s, weights = diff(c(0, F0_cdf))) # nolint: object_name_linter
  F1 <- BMisc::weighted_ecdf(trt_post$Y, weights = wnorm(trt_post)) # nolint: object_name_linter

  ptetools::attgt_noif(
    attgt = att,
    extra_gt_returns = list(F0 = F0, F1 = F1, Fte = NULL,
                            cdf_violation = max(abs(F0_cdf - F0_raw)))
  )
}


# --- Unified user-facing function --------------------------------------------

#' @title Difference-in-Differences in Distribution Functions
#'
#' @description Computes the Quantile Treatment Effect on the Treated (QTT)
#'   and the Average Treatment Effect on the Treated (ATT) under parallel
#'   trends in the distribution functions of untreated potential outcomes,
#'   with staggered treatment adoption, for panel data or repeated cross
#'   sections.
#'
#'   \strong{Identification.} For each group \eqn{g} and period \eqn{t},
#'   \deqn{F_{Y_t(0)|G=g}(y) - F_{Y_{t-1}(0)|G=g}(y) = F_{Y_t(0)|C}(y) -
#'   F_{Y_{t-1}(0)|C}(y) \quad \text{for all } y,}
#'   with \eqn{C} the comparison group: parallel trends for the distribution
#'   function rather than the mean (DiD) or the quantile function (QDiD).
#'   Roth and Sant'Anna (2023) show this is the condition under which
#'   parallel trends holds for every transformation of the outcome. It
#'   point-identifies the counterfactual distribution, so the QTT, for
#'   discrete outcomes too, and the implied ATT is the usual DiD ATT.
#'   Covariates make the assumption conditional on \eqn{X}.
#'
#'   See \code{\link{cdf_did_gt}} for how the counterfactual distribution is
#'   estimated (including the rearrangement into a valid distribution
#'   function and the three covariate methods).
#'
#' @inheritParams cdf_did_gt
#' @param yname Name of the outcome variable in \code{data}.
#' @param gname Name of the treatment group variable (first treatment period;
#'   0 for never-treated units).
#' @param tname Name of the time period variable.
#' @param idname Name of the unit id variable. Required when
#'   \code{panel = TRUE}.
#' @param data A data frame.
#' @param panel Logical. \code{TRUE} (default) for panel data; \code{FALSE}
#'   for repeated cross sections.
#' @param xformula One-sided formula of covariates. Default \code{~1}.
#' @param est_method With covariates: \code{"aipw"} (default; doubly
#'   robust), \code{"ipw"} or \code{"reg"} (outcome regression, as in
#'   \code{did}).
#' @param weightsname Name of the column in \code{data} containing sampling
#'   weights. Default \code{NULL} uses equal weights.
#' @param clustervars Name of a column in \code{data} to cluster standard
#'   errors on, in addition to \code{idname} (as in \code{did::att_gt()}).
#'   With the empirical bootstrap, whole clusters are resampled.
#' @param control_group \code{"notyettreated"} (default) or
#'   \code{"nevertreated"}.
#' @param anticipation Number of time periods before treatment during which
#'   units may anticipate it. Default 0.
#' @param alp Significance level. Default 0.05.
#' @param cband Logical; uniform confidence bands. Default \code{TRUE}.
#' @param biters Number of bootstrap iterations. Default 100.
#' @param cl Number of cores for the bootstrap. Default 1.
#' @param gt_type \code{"att"} (default) for the ATT and event study, or
#'   \code{"qtt"} for the QTT curve.
#' @param probs For \code{gt_type = "qtt"}, the quantiles at which to report
#'   the QTT. Default \code{seq(0.05, 0.95, 0.05)}.
#'
#' @section Discrete outcomes:
#'   The counterfactual distribution is identified for discrete outcomes as
#'   well, unlike for rank- or quantile-based estimators such as
#'   \code{\link{cic}}. Quantile effects still jump at mass points and their
#'   bootstrap inference can be unreliable, so \code{gt_type = "qtt"} warns
#'   (once per call) when the outcome looks discrete; see
#'   \code{\link{cic}}.
#'
#' @return For \code{gt_type = "att"}, a \code{pte_emp_boot} object; for
#'   \code{gt_type = "qtt"}, a \code{pte_qtt} object. Both include
#'   \code{outcome_discreteness}.
#'
#' @references
#'   Roth, Jonathan, and Pedro H. C. Sant'Anna. ``When Is Parallel Trends
#'   Sensitive to Functional Form?'' Econometrica 91(2), pp. 737-747, 2023.
#'
#'   Sant'Anna, Pedro H. C., and Jun Zhao. ``Doubly Robust
#'   Difference-in-Differences Estimators.'' Journal of Econometrics 219(1),
#'   pp. 101-122, 2020.
#'
#' @examples
#' \donttest{
#' data(mpdta, package = "did")
#' res <- cdf_did(yname = "lemp", gname = "first.treat", tname = "year",
#'                idname = "countyreal", data = mpdta, gt_type = "qtt",
#'                probs = seq(0.1, 0.9, 0.1), biters = 20)
#' summary(res)
#' }
#'
#' @export
cdf_did <- function(yname,
                    gname,
                    tname,
                    idname        = NULL,
                    data,
                    panel         = TRUE,
                    xformula      = ~1,
                    est_method    = c("aipw", "ipw", "reg"),
                    weightsname   = NULL,
                    clustervars   = NULL,
                    control_group = "notyettreated",
                    anticipation  = 0,
                    alp           = 0.05,
                    cband         = TRUE,
                    biters        = 100,
                    cl            = 1,
                    gt_type       = c("att", "qtt"),
                    probs         = NULL) {
  est_method <- match.arg(est_method)
  gt_type <- match.arg(gt_type)

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

  # warn once if the outcome looks discrete (#17), for QTTs only: the
  # counterfactual distribution and the ATT are identified for discrete
  # outcomes, but quantile effects jump at mass points
  disc <- check_discrete_outcome(
    data, yname, weightsname, warn = gt_type == "qtt",
    advice = paste0("cdf_did() still identifies the counterfactual ",
                    "distribution; interpret the quantile effects with care.")
  )

  res <- ptetools::pte(
    yname           = yname,
    gname           = gname,
    tname           = tname,
    idname          = idname,
    data            = keep_used_columns(data, yname, gname, tname, idname,
                                        weightsname, xformula,
                                        clustervars = clustervars),
    panel           = panel,
    setup_pte_fun   = ptetools::setup_pte,
    subset_fun      = subset_fun,
    attgt_fun       = cdf_did_gt,
    aggte_fun       = aggte_fun,
    xformula        = xformula,
    est_method      = est_method,
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
  res$outcome_discreteness <- disc
  res
}
