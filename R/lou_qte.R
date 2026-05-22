# =============================================================================
# Title: Lagged Outcome Unconfoundedness QTT (lou_qte)
# Description: Estimates the QTT and ATT under a lagged-outcome unconfoundedness
#   assumption with staggered treatment adoption. Operates at the (g,t) level
#   via ptetools::pte() and supports IPW, outcome regression (OR), and doubly-
#   robust AIPW estimation. When lagged_outcome_cov = TRUE (the default), the
#   pre-period outcome is added as a covariate in both the propensity score
#   model and the outcome model, requiring panel data. When FALSE, the function
#   works with either panel or repeated cross sections.
#
#   The ATT version of this approach for staggered treatment adoption appears in:
#   Callaway, Brantly. "Policy Evaluation during a Pandemic."
#   Journal of Econometrics 236(2), 2023.
# Author: Brant Callaway
# Last update: 2026-05-22
# Date created: 2026-05-22
# =============================================================================

# --- (g,t)-level workhorse ---------------------------------------------------

#' @title lou_gt
#'
#' @description Estimates the ATT and counterfactual distribution F0 = F_{Y(0)|D=1}
#'   for a single (g,t) cell under lagged-outcome unconfoundedness. Called
#'   internally by \code{\link{lou_qte}} via \code{ptetools::pte()}.
#'
#' @param gt_data A \code{gt_data_frame} from \code{ptetools::two_by_two_subset}
#'   or \code{ptetools::two_by_two_rcs_subset}. Contains columns \code{Y},
#'   \code{D}, \code{name} ("pre"/"post"), \code{id}, and \code{.w}.
#' @param xformula One-sided covariate formula. Default \code{~1}.
#' @param lagged_outcome_cov Logical; if \code{TRUE} (default), the pre-period
#'   outcome is appended to \code{xformula} as an additional covariate. Requires
#'   panel data (unit ids present in both periods).
#' @param est_method One of \code{"ipw"}, \code{"or"}, or \code{"aipw"}.
#'   Default \code{"ipw"}.
#' @param ... Additional arguments (ignored).
#'
#' @keywords internal
#'
#' @return A list from \code{ptetools::attgt_noif} with \code{attgt} (ATT) and
#'   \code{extra_gt_returns} containing \code{F0} and \code{F1} (weighted ECDFs
#'   of counterfactual and observed treated outcomes).
lou_gt <- function(gt_data, xformula = ~1, lagged_outcome_cov = TRUE,
                   est_method = "ipw", ...) {
  gt_data <- droplevels(gt_data)

  pre_data  <- gt_data[gt_data$name == "pre",  ]
  post_data <- gt_data[gt_data$name == "post", ]
  post_trt  <- post_data[post_data$D == 1, ]
  post_ctrl <- post_data[post_data$D == 0, ]

  if (nrow(post_trt) == 0L || nrow(post_ctrl) == 0L) {
    return(ptetools::attgt_noif(attgt = NA))
  }

  # Attach lagged outcome from pre period; extend xformula
  if (isTRUE(lagged_outcome_cov)) {
    pre_y <- pre_data$Y[match(post_data$id, pre_data$id)]
    # Drop units with no pre-period observation (unbalanced panel)
    keep      <- !is.na(pre_y)
    post_data <- post_data[keep, ]
    pre_y     <- pre_y[keep]
    post_data$Y_lag <- pre_y
    post_trt  <- post_data[post_data$D == 1, ]
    post_ctrl <- post_data[post_data$D == 0, ]
    xformula_use <- BMisc::add_cov_to_formula("Y_lag", xformula)
  } else {
    xformula_use <- xformula
  }

  if (nrow(post_trt) == 0L || nrow(post_ctrl) == 0L) {
    return(ptetools::attgt_noif(attgt = NA))
  }

  # Normalised weights within each (period x treatment) cell
  Y_post_trt  <- post_trt$Y  # nolint: object_name_linter
  Y_post_ctrl <- post_ctrl$Y # nolint: object_name_linter
  w_post_trt  <- post_trt$.w  / sum(post_trt$.w)  # nolint: object_name_linter
  w_post_ctrl <- post_ctrl$.w / sum(post_ctrl$.w) # nolint: object_name_linter

  has_covs <- length(BMisc::rhs_vars(xformula_use)) > 0
  tau_grid <- seq(0.01, 0.99, 0.01)
  K        <- length(tau_grid)

  # Without covariates all three methods reduce to a simple distribution
  # comparison; the unconfoundedness assumption must hold unconditionally.
  if (!has_covs) {
    F1  <- BMisc::weighted_ecdf(Y_post_trt,  weights = w_post_trt)  # nolint: object_name_linter
    F0  <- BMisc::weighted_ecdf(Y_post_ctrl, weights = w_post_ctrl) # nolint: object_name_linter
    att <- weighted.mean(Y_post_trt, w_post_trt) -
           weighted.mean(Y_post_ctrl, w_post_ctrl)
    return(ptetools::attgt_noif(attgt = att,
                                extra_gt_returns = list(F0 = F0, F1 = F1,
                                                        Fte = NULL)))
  }

  # --- IPW -------------------------------------------------------------------
  if (est_method == "ipw") {
    ps_formla  <- BMisc::toformula("D", BMisc::rhs_vars(xformula_use))
    pscore_reg <- glm(ps_formla, data = post_data, family = binomial("logit"))
    pscore     <- pmax(1e-6, pmin(1 - 1e-6, fitted(pscore_reg)))

    # QTT: reweight controls by p(x)/(1-p(x)) to match treated covariate dist
    ps_ctrl <- pscore[post_data$D == 0]
    u_wts   <- post_ctrl$.w * ps_ctrl / (1 - ps_ctrl)
    u_wts   <- u_wts / sum(u_wts)

    F1  <- BMisc::weighted_ecdf(Y_post_trt,  weights = w_post_trt)  # nolint: object_name_linter
    F0  <- BMisc::weighted_ecdf(Y_post_ctrl, weights = u_wts)       # nolint: object_name_linter
    att <- weighted.mean(Y_post_trt, w_post_trt) -
           weighted.mean(Y_post_ctrl, u_wts)

  # --- OR --------------------------------------------------------------------
  } else if (est_method == "or") {
    yformla <- BMisc::toformula("Y", BMisc::rhs_vars(xformula_use))
    # Use .w by column name so rq() finds it inside the data frame;
    # avoids parent.frame() lookup failure when called inside suppressWarnings
    # Suppress "Solution may be nonunique" from quantreg LP solver — benign
    M0 <- suppressWarnings(
      quantreg::rq(yformla, data = post_ctrl, tau = tau_grid, weights = .w)
    )
    # Predict counterfactual conditional quantiles at treated units; monotonize
    pred0  <- t(apply(as.matrix(predict(M0, newdata = post_trt)), 1, sort))
    wt_rep <- rep(post_trt$.w, each = K)
    wt_rep <- wt_rep / sum(wt_rep)

    F0  <- BMisc::weighted_ecdf(c(t(pred0)), weights = wt_rep)      # nolint: object_name_linter
    F1  <- BMisc::weighted_ecdf(Y_post_trt,  weights = w_post_trt)  # nolint: object_name_linter
    att <- weighted.mean(Y_post_trt, w_post_trt) -
           weighted.mean(c(t(pred0)), wt_rep)

  # --- AIPW ------------------------------------------------------------------
  } else if (est_method == "aipw") {
    ps_formla  <- BMisc::toformula("D", BMisc::rhs_vars(xformula_use))
    pscore_reg <- glm(ps_formla, data = post_data, family = binomial("logit"))
    pscore     <- pmax(1e-6, pmin(1 - 1e-6, fitted(pscore_reg)))

    yformla <- BMisc::toformula("Y", BMisc::rhs_vars(xformula_use))
    M0 <- suppressWarnings(
      quantreg::rq(yformla, data = post_ctrl, tau = tau_grid, weights = .w)
    )
    # Predict at all post-period units for the augmentation term
    pred0_all <- t(apply(as.matrix(predict(M0, newdata = post_data)), 1, sort))

    D      <- post_data$D
    y      <- post_data$Y
    pbar   <- weighted.mean(D, post_data$.w)
    y_grid <- sort(unique(y))

    # Doubly-robust CDF for Y(0)|D=1:
    #   F(y) = E[D*F_cond(y) + (1-D)*(p/(1-p))*(1{Y<=y} - F_cond(y))] / E[D]
    F_qtt <- vapply(y_grid, function(yy) {
      F_cond    <- rowMeans(pred0_all <= yy)
      aipw_term <- D * F_cond +
                   (1 - D) * pscore / (1 - pscore) *
                   (as.numeric(y <= yy) - F_cond)
      weighted.mean(aipw_term, post_data$.w) / pbar
    }, numeric(1))
    F_qtt <- cummax(pmax(0, pmin(1, F_qtt)))

    # Represent as ecdf via probability masses at y_grid support points
    F_qtt_mass <- pmax(0, diff(c(0, F_qtt)))
    if (sum(F_qtt_mass) == 0) F_qtt_mass <- rep(1, length(y_grid))
    F0 <- BMisc::weighted_ecdf(y_grid,      # nolint: object_name_linter
                               weights = F_qtt_mass / sum(F_qtt_mass))
    F1 <- BMisc::weighted_ecdf(Y_post_trt,  # nolint: object_name_linter
                               weights = w_post_trt)

    # ATT via IPW augmentation (consistent under correct pscore)
    ps_ctrl <- pscore[D == 0]
    u_wts   <- post_ctrl$.w * ps_ctrl / (1 - ps_ctrl)
    u_wts   <- u_wts / sum(u_wts)
    att     <- weighted.mean(Y_post_trt, w_post_trt) -
               weighted.mean(Y_post_ctrl, u_wts)
  }

  ptetools::attgt_noif(attgt = att,
                       extra_gt_returns = list(F0 = F0, F1 = F1, Fte = NULL))
}


# --- User-facing function ----------------------------------------------------

#' @title Lagged Outcome Unconfoundedness QTT
#'
#' @description Estimates the Quantile Treatment Effect on the Treated (QTT)
#'   and Average Treatment Effect on the Treated (ATT) under a lagged-outcome
#'   unconfoundedness assumption with staggered treatment adoption. The key
#'   identifying assumption is
#'   \eqn{Y_{g,t}(0) \perp D \mid X, Y_{\text{pre}}}, i.e., conditional on
#'   observed covariates and the pre-treatment outcome, treatment is as good
#'   as randomly assigned within each cohort-period cell.
#'
#'   Estimation operates at the (g,t) level: for each cohort \eqn{g} and
#'   post-treatment period \eqn{t}, a cross-sectional comparison is made
#'   between the treated group and a not-yet-treated (or never-treated)
#'   comparison group, adjusting for covariates and (optionally) the
#'   pre-treatment outcome. Group-time estimates are then aggregated to
#'   overall, dynamic, and group-specific summaries.
#'
#'   Three estimation methods are available:
#'   \describe{
#'     \item{\code{"ipw"}}{Propensity-score reweighting. Control units are
#'       reweighted by \eqn{\hat p(X,Y_{\text{pre}})/(1-\hat p(X,Y_{\text{pre}}))}
#'       to approximate the covariate distribution of the treated group.}
#'     \item{\code{"or"}}{Outcome regression. A quantile regression model is
#'       fit on control units' post-period outcomes as a function of
#'       \eqn{(X, Y_{\text{pre}})}, then predicted at treated units to
#'       construct the counterfactual distribution (Melly 2006;
#'       Chernozhukov, Fernandez-Val, and Melly 2013).}
#'     \item{\code{"aipw"}}{Doubly-robust augmented IPW. Combines the
#'       propensity score and outcome models. Consistent if either model is
#'       correctly specified.}
#'   }
#'
#'   When \code{lagged_outcome_cov = FALSE} and \code{xformula = ~1}, all
#'   three methods reduce to a simple distribution comparison within each
#'   (g,t) cell, which is consistent under unconditional unconfoundedness.
#'
#'   The ATT analogue (without the QTT) for staggered adoption under
#'   lagged-outcome unconfoundedness is developed in Callaway (2023).
#'
#' @param yname Name of the outcome variable in \code{data}.
#' @param gname Name of the treatment group variable (first treatment period;
#'   0 for never-treated units).
#' @param tname Name of the time period variable.
#' @param idname Name of the unit id variable. Required when
#'   \code{panel = TRUE}.
#' @param data A data frame.
#' @param xformula One-sided formula for additional covariates used in the
#'   propensity score and/or outcome model. Default \code{~1} uses no
#'   additional covariates beyond the lagged outcome (when
#'   \code{lagged_outcome_cov = TRUE}).
#' @param lagged_outcome_cov Logical; if \code{TRUE} (default), the
#'   pre-treatment outcome is added as a covariate in both the propensity
#'   score and outcome models. Requires \code{panel = TRUE}.
#' @param est_method Estimation method: \code{"ipw"} (default),
#'   \code{"or"}, or \code{"aipw"}. See Details.
#' @param panel Logical; \code{TRUE} (default) for panel data,
#'   \code{FALSE} for repeated cross sections. \code{lagged_outcome_cov = TRUE}
#'   requires \code{panel = TRUE}.
#' @param weightsname Name of the column in \code{data} containing sampling
#'   weights. Default \code{NULL} uses equal weights.
#' @param control_group Which units to use as the comparison group:
#'   \code{"notyettreated"} (default) or \code{"nevertreated"}.
#' @param anticipation Number of periods of anticipation. Default \code{0}.
#' @param alp Significance level for confidence intervals. Default \code{0.05}.
#' @param cband Logical; if \code{TRUE} (default), report a simultaneous
#'   confidence band (uniform over all quantiles in \code{probs}) in addition
#'   to pointwise intervals.
#' @param biters Number of bootstrap iterations. Default \code{100}.
#' @param cl Number of clusters for parallel bootstrap. Default \code{1}.
#' @param gt_type Type of group-time effect to compute. \code{"att"} (default)
#'   returns the ATT curve aggregated over (g,t) cells. \code{"qtt"} returns
#'   the full QTT curve over \code{probs} using mixture-CDF aggregation.
#' @param probs For \code{gt_type = "qtt"}, the quantile grid at which to
#'   evaluate the QTT. Default \code{seq(0.05, 0.95, 0.05)}.
#'
#' @return For \code{gt_type = "att"}, a \code{pte_emp_boot} object. For
#'   \code{gt_type = "qtt"}, a \code{pte_qtt} object with overall,
#'   group-specific, and dynamic QTT curves, bootstrap standard errors,
#'   and pointwise and uniform confidence bands.
#'
#' @references
#'   Callaway, Brantly. ``Policy Evaluation during a Pandemic.''
#'   \emph{Journal of Econometrics} 236(2), 2023.
#'
#'   Melly, Blaise. ``Estimation of Counterfactual Distributions Using
#'   Quantile Regression.'' Working paper, University of St. Gallen, 2006.
#'
#'   Chernozhukov, Victor, Ivan Fernandez-Val, and Blaise Melly.
#'   ``Inference on Counterfactual Distributions.''
#'   \emph{Econometrica} 81(6), pp. 2205--2268, 2013.
#'
#' @seealso \code{\link{unc_qte}} for the cross-sectional (non-staggered)
#'   version. \code{\link{cic}}, \code{\link{qdid}}, \code{\link{mdid}} for
#'   alternative identification strategies with staggered adoption.
#'
#' @examples
#' \donttest{
#' data(mpdta, package = "did")
#'
#' ## ATT under lagged-outcome unconfoundedness (IPW with pre-period outcome)
#' res_att <- lou_qte(
#'   yname = "lemp", gname = "first.treat", tname = "year",
#'   idname = "countyreal", data = mpdta,
#'   lagged_outcome_cov = TRUE, est_method = "ipw",
#'   gt_type = "att", biters = 50
#' )
#' summary(res_att)
#'
#' ## QTT with doubly-robust estimation
#' res_qtt <- lou_qte(
#'   yname = "lemp", gname = "first.treat", tname = "year",
#'   idname = "countyreal", data = mpdta,
#'   lagged_outcome_cov = TRUE, est_method = "aipw",
#'   gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 50
#' )
#' summary(res_qtt)
#' }
#'
#' @export
lou_qte <- function(yname, gname, tname, idname = NULL, data,
                    xformula          = ~1,
                    lagged_outcome_cov = TRUE,
                    est_method        = c("ipw", "or", "aipw"),
                    panel             = TRUE,
                    weightsname       = NULL,
                    control_group     = "notyettreated",
                    anticipation      = 0,
                    alp               = 0.05,
                    cband             = TRUE,
                    biters            = 100,
                    cl                = 1,
                    gt_type           = c("att", "qtt"),
                    probs             = NULL) {
  est_method <- match.arg(est_method)
  gt_type    <- match.arg(gt_type)

  if (isTRUE(lagged_outcome_cov) && !isTRUE(panel)) {
    stop("lagged_outcome_cov = TRUE requires panel data; set panel = TRUE or lagged_outcome_cov = FALSE")
  }

  if (isTRUE(panel) && is.null(idname)) {
    stop("idname must be provided when panel = TRUE")
  }

  subset_fun <- if (isTRUE(panel)) {
    ptetools::two_by_two_subset
  } else {
    ptetools::two_by_two_rcs_subset
  }

  aggregation_fun <- if (gt_type == "qtt") {
    ptetools::qtt_pte_aggregations
  } else {
    function(al, p, eg) ptetools::attgt_pte_aggregations(al, p)
  }

  ptetools::pte(
    yname              = yname,
    gname              = gname,
    tname              = tname,
    idname             = idname,
    data               = data,
    panel              = panel,
    setup_pte_fun      = ptetools::setup_pte,
    subset_fun         = subset_fun,
    attgt_fun          = lou_gt,
    aggregation_fun    = aggregation_fun,
    xformula           = xformula,
    lagged_outcome_cov = lagged_outcome_cov,
    est_method         = est_method,
    weightsname        = weightsname,
    control_group      = control_group,
    anticipation       = anticipation,
    cband              = cband,
    alp                = alp,
    boot_type          = "empirical",
    biters             = biters,
    cl                 = cl,
    gt_type            = gt_type,
    probs              = probs
  )
}
