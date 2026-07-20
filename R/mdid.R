# =============================================================================
# Title: Mean Difference-in-Differences (MDiD)
# Description: Implements the MDiD estimator (Athey-Imbens 2006, Thuysbaert 2007).
#   mdid_gt is the 2x2 workhorse used directly as a ptetools attgt_fun.
#   mdid is the unified user-facing function covering two-period and staggered
#   treatment adoption for both panel and repeated cross sections data.
#   MDiD is a deprecated alias retained for backwards compatibility.
#   compute.MDiD is an unexported helper kept only to support the deprecated MDiD.
# Author: Brant Callaway
# Last update: 2026-05-22
# Date created: 2026-05-18
# =============================================================================

# --- 2x2 workhorse -----------------------------------------------------------

#' @title Mean Difference-in-Differences: group-time estimator
#'
#' @description Computes the MDiD ATT and counterfactual outcome distribution
#'   for a single 2x2 (pre/post x treated/control) data subset. Serves
#'   directly as the \code{attgt_fun} argument to \code{ptetools::pte}. Panel
#'   vs. repeated cross sections is detected automatically from whether the
#'   same unit ids appear in both periods.
#'
#'   \strong{Identification.} MDiD assumes that the counterfactual distribution
#'   of untreated potential outcomes for the treated group in the post period is
#'   a location shift of the treated group's pre-period distribution:
#'   \deqn{F_{Y(0),\text{post}|D=1}(y) = F_{Y,\text{pre}|D=1}(y - \Delta)}
#'   where \eqn{\Delta = E[Y_{\text{post}}|D=0] - E[Y_{\text{pre}}|D=0]} is
#'   the mean DiD (the change in mean outcomes for the untreated group). The
#'   counterfactual quantile function is therefore
#'   \eqn{Q_{Y(0),\text{post}|D=1}(\tau) = Q_{Y,\text{pre}|D=1}(\tau) + \Delta}.
#'
#'   \strong{Covariate adjustment.} With covariates, the scalar shift \eqn{\Delta}
#'   is replaced by a unit-specific conditional mean shift
#'   \eqn{\Delta(X_i) = E[Y_{\text{post}}|D=0, X_i] - E[Y_{\text{pre}}|D=0, X_i]},
#'   estimated by weighted OLS on the control group in each period. The
#'   counterfactual for treated pre-period unit \eqn{i} is
#'   \eqn{Y_{\text{pre},i} + \Delta(X_i)}. The unconditional counterfactual
#'   distribution is the empirical CDF of these shifted values, which by the
#'   law of iterated expectations consistently estimates
#'   \eqn{F_{Y(0),\text{post}|D=1}(y) = \int F_{Y,\text{pre}|D=1,X=x}(y - \Delta(x))\,dF_{X|D=1}(x)}.
#'
#' @param gt_data A data frame (typically a \code{gt_data_frame} from
#'   \code{ptetools}) with columns \code{name} (\code{"pre"} or \code{"post"}),
#'   \code{D} (treatment dummy), \code{Y} (outcome), \code{id} (unit
#'   identifier), \code{.w} (sampling weights), and any covariate columns
#'   referenced by \code{xformula}.
#' @param xformula One-sided formula for covariates. Default \code{~1} uses no
#'   covariates. With covariates, separate weighted OLS regressions are fit on
#'   the control group in each period; see Details above.
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
#'   Thuysbaert, Bram. ``Distributional Comparisons in Difference in Differences
#'   Models.'' Working Paper, 2007.
#'
#' @export
mdid_gt <- function(gt_data, xformula = ~1, ...) {
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

  # MDiD transform: shift each treated pre-period outcome by delta, the mean DiD.
  # delta is a scalar (unconditional) or a vector of unit-specific conditional
  # mean shifts (with covariates). Both cases unify in kcf = Y_pre_trt + delta.
  delta <- weighted.mean(Y_post_ctrl, w_post_ctrl) - # nolint: object_name_linter
             weighted.mean(Y_pre_ctrl, w_pre_ctrl)

  # detect panel: same unit ids appear in both periods for the treated group
  trt_pre_ids  <- gt_data$id[gt_data$D == 1 & gt_data$name == "pre"]
  trt_post_ids <- gt_data$id[gt_data$D == 1 & gt_data$name == "post"]
  panel <- length(intersect(trt_pre_ids, trt_post_ids)) > 0

  # Covariate adjustment: replace the unconditional scalar shift with a
  # unit-specific conditional mean shift delta(X_i) = E[Y_post,ctrl|X_i] -
  # E[Y_pre,ctrl|X_i], estimated by weighted OLS on each control-group cell.
  # Justified by the law of iterated expectations: the empirical CDF of
  # Y_pre_trt[i] + delta(X_i) consistently estimates F_{Y(0),post|D=1}.
  lm_pre_ctrl  <- NULL
  lm_post_ctrl <- NULL
  if (length(BMisc::rhs_vars(xformula)) > 0) {
    pre_ctrl  <- gt_data[gt_data$name == "pre"  & gt_data$D == 0, ]
    post_ctrl <- gt_data[gt_data$name == "post" & gt_data$D == 0, ]
    pre_trt   <- gt_data[gt_data$name == "pre"  & gt_data$D == 1, ]

    yformla      <- BMisc::toformula("Y", BMisc::rhs_vars(xformula)) # nolint: object_name_linter
    lm_pre_ctrl  <- lm(yformla, data = pre_ctrl,  weights = pre_ctrl$.w)  # nolint: object_name_linter
    lm_post_ctrl <- lm(yformla, data = post_ctrl, weights = post_ctrl$.w) # nolint: object_name_linter

    delta <- predict(lm_post_ctrl, newdata = pre_trt) - # nolint: object_name_linter
               predict(lm_pre_ctrl,  newdata = pre_trt)
  }

  kcf <- Y_pre_trt + delta

  att <- weighted.mean(Y_post_trt, w_post_trt) - weighted.mean(kcf, w_pre_trt)
  F0  <- BMisc::weighted_ecdf(kcf,        weights = w_pre_trt)  # nolint: object_name_linter
  F1  <- BMisc::weighted_ecdf(Y_post_trt, weights = w_post_trt) # nolint: object_name_linter

  # distribution of individual treatment effects (panel only);
  # requires unit-level matching so pre and post are sorted by id
  Fte <- NULL # nolint: object_name_linter
  if (panel) {
    pre_s  <- gt_data[gt_data$name == "pre"  & gt_data$D == 1, ]
    post_s <- gt_data[gt_data$name == "post" & gt_data$D == 1, ]
    pre_s  <- pre_s[order(pre_s$id), ]
    post_s <- post_s[order(post_s$id), ]

    if (!is.null(lm_pre_ctrl)) {
      # recompute delta on id-sorted pre_s so unit i's shift aligns with post_s[i]
      delta_s <- predict(lm_post_ctrl, newdata = pre_s) - # nolint: object_name_linter
                   predict(lm_pre_ctrl,  newdata = pre_s)
    } else {
      delta_s <- delta  # scalar; same shift for all units
    }

    kcf_s    <- pre_s$Y + delta_s
    w_post_s <- post_s$.w / sum(post_s$.w)
    Fte <- BMisc::weighted_ecdf(post_s$Y - kcf_s, weights = w_post_s) # nolint: object_name_linter
  }

  ptetools::attgt_noif(attgt = att, extra_gt_returns = list(F0 = F0, F1 = F1, Fte = Fte))
}


# --- Unified user-facing function --------------------------------------------

#' @title Mean Difference-in-Differences
#'
#' @description Computes Quantile Treatment effects on the Treated (QTT) and
#'   the Average Treatment Effect on the Treated (ATT) using the Mean
#'   Difference-in-Differences identification strategy. Handles two-period data
#'   and staggered treatment adoption uniformly. Supports both panel and
#'   repeated cross sections data.
#'
#'   \strong{Identification.} MDiD assumes the counterfactual distribution of
#'   untreated potential outcomes is a location shift of the treated group's
#'   pre-treatment distribution. The size of the shift is the mean DiD
#'   \eqn{\Delta = E[Y_{\text{post}}|D=0] - E[Y_{\text{pre}}|D=0]}, so
#'   \eqn{Q_{Y(0),\text{post}|D=1}(\tau) = Q_{Y,\text{pre}|D=1}(\tau) + \Delta}.
#'   This is stronger than parallel trends in means alone: it additionally
#'   requires that the shape of the treated group's outcome distribution is
#'   unchanged in the counterfactual. MDiD is a special case of QDiD that
#'   applies a single mean shift rather than a rank-specific distributional
#'   shift.
#'
#'   \strong{Covariate adjustment.} When \code{xformula} is specified, the
#'   scalar shift is replaced by a unit-specific conditional mean shift
#'   \eqn{\Delta(X_i) = E[Y_{\text{post}}|D=0,X_i] - E[Y_{\text{pre}}|D=0,X_i]},
#'   estimated by weighted OLS. The counterfactual for treated unit \eqn{i} is
#'   \eqn{Y_{\text{pre},i} + \Delta(X_i)}, and the unconditional counterfactual
#'   distribution is the empirical CDF of these values. By the law of iterated
#'   expectations, this consistently estimates \eqn{F_{Y(0),\text{post}|D=1}}.
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
#'   Thuysbaert, Bram. ``Distributional Comparisons in Difference in Differences
#'   Models.'' Working Paper, 2007.
#'
#' @examples
#' \donttest{
#' data(mpdta, package = "did")
#'
#' ## ATT aggregated across all groups and periods
#' res_att <- mdid(yname = "lemp", gname = "first.treat", tname = "year",
#'                 idname = "countyreal", data = mpdta,
#'                 gt_type = "att", biters = 20)
#' summary(res_att)
#'
#' ## Full QTT curve at selected quantiles
#' res_qtt <- mdid(yname = "lemp", gname = "first.treat", tname = "year",
#'                 idname = "countyreal", data = mpdta,
#'                 gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 20)
#' summary(res_qtt)
#' }
#'
#' @export
mdid <- function(yname,
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
    data            = data,
    panel           = panel,
    setup_pte_fun   = ptetools::setup_pte,
    subset_fun      = subset_fun,
    attgt_fun       = mdid_gt,
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


# --- Deprecated two-period workhorse (supports MDiD only) -------------------

# Kept only to support the deprecated MDiD() wrapper below. Do not use in new
# code; use mdid_gt() instead.
#' @keywords internal
# nolint start: object_name_linter
compute.MDiD <- function(formla, xformla = NULL, t, tmin1, tname, data,
                         panel = FALSE,
                         idname = NULL, probs = seq(0.05, 0.95, 0.05)) {
  form <- as.formula(formla)
  dta <- model.frame(terms(form, data = data), data = data)
  colnames(dta) <- c("y", "treatment")
  yname <- "y"
  treat <- "treatment"
  data <- cbind.data.frame(dta, data)

  x <- NULL
  if (!(is.null(xformla))) {
    x <- colnames(model.matrix(terms(as.formula(xformla)), data = data))[-1]
    data <- cbind(
      data[, c(yname, treat, idname, tname)],
      model.matrix(terms(as.formula(xformla)), data = data)
    )[, c(1:4, 6:(5 + length(x)))]
  }

  data <- subset(data, (data[, tname] == tmin1 | data[, tname] == t))

  if (panel) {
    data <- make_balanced_panel(data, idname, tname)
  }

  data <- droplevels(data)

  if (!(is.null(x))) {
    cov.data <- data
    cov.data$group <- as.factor(paste(cov.data[, treat], cov.data[, tname], sep = "-"))
    xmat <- cov.data[, x]
    first.stage <- lm(cov.data[, yname] ~ -1 + cov.data[, cov.data$group] +
      as.matrix(xmat))
    bet <- coef(first.stage)[5:length(coef(first.stage))]
    yfit <- cov.data[, yname] - as.matrix(xmat) %*% bet
    data[, yname] <- yfit
  }

  treated.t    <- data[data[, tname] == t    & data[, treat] == 1, ]
  treated.tmin1  <- data[data[, tname] == tmin1 & data[, treat] == 1, ]
  untreated.t    <- data[data[, tname] == t    & data[, treat] == 0, ]
  untreated.tmin1 <- data[data[, tname] == tmin1 & data[, treat] == 0, ]

  q1  <- quantile(treated.t[, yname], probs = probs)
  q0  <- quantile(treated.tmin1[, yname], probs = probs) +
           mean(untreated.t[, yname]) - mean(untreated.tmin1[, yname])
  att <- mean(treated.t[, yname]) -
           (mean(treated.tmin1[, yname]) + mean(untreated.t[, yname]) -
            mean(untreated.tmin1[, yname]))

  QTE(ate = att, qte = (q1 - q0), probs = probs)
}
# nolint end: object_name_linter


# --- Deprecated two-period wrapper -------------------------------------------

#' @title MDiD
#'
#' @description \strong{Deprecated.} Use \code{\link{mdid}} instead.
#'
#'   \code{MDiD} is a legacy two-period wrapper for the Mean
#'   Difference-in-Differences estimator. The modern replacement
#'   \code{\link{mdid}} uses a \code{yname}/\code{gname}/\code{tname}
#'   interface and supports staggered treatment adoption natively.
#'
#' @param formla formula of the form \code{y ~ treat}.
#' @param xformla optional one-sided covariate formula (unused in this
#'   estimator; accepted for consistency with other wrappers).
#' @param t post-treatment period.
#' @param tmin1 pre-treatment period.
#' @param tname name of the time column.
#' @param data data.frame.
#' @param panel logical; whether data are panel.
#' @param se logical; whether to compute bootstrap SEs.
#' @param idname name of the unit id column (required when \code{panel = TRUE}).
#' @param alp significance level.
#' @param probs quantile grid.
#' @param iters bootstrap iterations.
#' @param retEachIter logical; return per-iteration results.
#' @param biters alias for \code{iters}; takes precedence if supplied.
#'
#' @return \code{QTE} object.
#'
#' @seealso \code{\link{mdid}}
#'
#' @name MDiD-deprecated
#' @aliases MDiD
#' @keywords internal
#' @export
# nolint start: object_name_linter
MDiD <- function(formla, xformla = NULL, t, tmin1, tname, data,
                 panel = FALSE, se = TRUE,
                 idname = NULL,
                 alp = 0.05, probs = seq(0.05, 0.95, 0.05), iters = 100,
                 retEachIter = FALSE,
                 biters = NULL) {
  .Deprecated(msg = paste0(
    "MDiD() is deprecated. The replacement mdid() uses a different data model: ",
    "specify outcome, group, and time variables by name (yname, gname, tname) ",
    "and provide gname as a column holding each unit's first treatment period ",
    "(0 = never-treated), rather than a formula and explicit t/tmin1 values. ",
    "See ?mdid."
  ))
  # biters is an alias matching the did/ptetools API (MDiD has no parallel support)
  if (!is.null(biters)) iters <- biters

  form <- as.formula(formla)
  dta <- model.frame(terms(form, data = data), data = data)
  colnames(dta) <- c("y", "treatment")
  yname <- "y"
  treat <- "treatment"
  data <- cbind.data.frame(dta, data)

  data <- subset(data, (data[, tname] == tmin1 | data[, tname] == t))

  if (panel) {
    if (is.null(idname)) stop("Must provide idname when using panel option")
    data <- make_balanced_panel(data, idname, tname)
  }

  data <- droplevels(data)

  treated.t      <- data[data[, tname] == t    & data[, treat] == 1, ]
  treated.tmin1  <- data[data[, tname] == tmin1 & data[, treat] == 1, ]
  untreated.t    <- data[data[, tname] == t    & data[, treat] == 0, ]
  untreated.tmin1 <- data[data[, tname] == tmin1 & data[, treat] == 0, ]

  mdid_result <- compute.MDiD(formla, xformla, t, tmin1, tname, data, panel, idname, probs) # nolint: object_name_linter

  if (se) {
    eachIter <- list()
    if (panel) {
      treated.t      <- treated.t[order(treated.t[, idname]), ]
      treated.tmin1  <- treated.tmin1[order(treated.tmin1[, idname]), ]
      untreated.t    <- untreated.t[order(untreated.t[, idname]), ]
      untreated.tmin1 <- untreated.tmin1[order(untreated.tmin1[, idname]), ]
      nt <- nrow(treated.t)
      nu <- nrow(untreated.t)
      for (i in seq_len(iters)) {
        randy.t <- sample(seq_len(nt), nt, replace = TRUE)
        randy.u <- sample(seq_len(nu), nu, replace = TRUE)
        boot.data <- rbind(
          treated.t[randy.t, ],    untreated.t[randy.u, ],
          treated.tmin1[randy.t, ], untreated.tmin1[randy.u, ]
        )
        thisIter <- compute.MDiD(formla, xformla, t, tmin1, tname, boot.data, # nolint: object_name_linter
                                  panel = FALSE, idname, probs)
        eachIter[[i]] <- QTE(ate = thisIter$ate, qte = thisIter$qte, probs = probs) # nolint: object_name_linter
      }
    } else {
      for (i in seq_len(iters)) {
        boot.data <- rbind(
          treated.t[sample(nrow(treated.t),       replace = TRUE), ],
          treated.tmin1[sample(nrow(treated.tmin1), replace = TRUE), ],
          untreated.t[sample(nrow(untreated.t),    replace = TRUE), ],
          untreated.tmin1[sample(nrow(untreated.tmin1), replace = TRUE), ]
        )
        thisIter <- compute.MDiD(formla, xformla, t, tmin1, tname, boot.data, panel, idname, probs) # nolint: object_name_linter
        eachIter[[i]] <- QTE(ate = thisIter$ate, qte = thisIter$qte, probs = probs) # nolint: object_name_linter
      }
    }

    SEobj <- computeSE(eachIter, mdid_result, alp = alp) # nolint: object_name_linter

    QTE( # nolint: object_name_linter
      qte = mdid_result$qte, qte.upper = SEobj$qte.upper,
      qte.lower = SEobj$qte.lower, ate = mdid_result$ate,
      ate.upper = SEobj$ate.upper, ate.lower = SEobj$ate.lower,
      qte.se = SEobj$qte.se, ate.se = SEobj$ate.se,
      alp = alp,
      eachIterList = if (retEachIter) eachIter else NULL,
      probs = probs
    )
  } else {
    mdid_result
  }
}
# nolint end: object_name_linter
